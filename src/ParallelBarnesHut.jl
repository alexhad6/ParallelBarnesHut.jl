module ParallelBarnesHut

# Imports

using Plots

using Random


# Exports

export Particle

export simulate!, grav_acc

export show_particles, show_quadview, show_boxes, animate_frames

export rand_particles


# Structs

"""
Represents a particle in a Barnes Hut simulation. For a galaxy simulation, a
particle is a star or a group of stars.
"""
mutable struct Particle
    "Position coordinates of a particle (in meters)."
    pos::Array{Float64,1}
    
    "Velocity vector of a particle (in meters per second)."
    vel::Array{Float64,1}
    
    "Mass of a particle (in kilograms)."
    mass::Float64
end

"""
Used in the octree for the Barnes Hut algorithm. Represents a cubical octant of
space and stores the mass and center of mass of the contained particles.
"""
struct Node
    "Side length of the cubical region of space represented by a node."
    size::Float64
    
    "Total mass of the particles contained in the region represented by a node."
    total_mass::Float64
    
    "Center of mass of the particles contained in the region represented by a node."
    center_of_mass::Array{Float64,1}
    
    "Children of a node in the octree."
    children::Array{Union{Node, Nothing},1}
end


# Barnes Hut Tree Generation

"""
    generate_tree(particles; return_boxes = false)

Generate the Barnes Hut octree for a given array of particles. Optionally
return coordinates of the boxes represented by the nodes in the octree.
"""
function generate_tree(particles::Array{Particle,1}; return_boxes::Bool = false)
    # Calculate bounds of a box enclosing all particles
    minx, maxx = minimum(p -> p.pos[1], particles), maximum(p -> p.pos[1], particles)
    miny, maxy = minimum(p -> p.pos[2], particles), maximum(p -> p.pos[2], particles)
    minz, maxz = minimum(p -> p.pos[3], particles), maximum(p -> p.pos[3], particles)
    
    # Calculate corner and size of a cube enclosing all particles
    size = max(maxx - minx, maxy - miny, maxz - minz)
    corner = [minx, miny, minz]
    
    if return_boxes
        # Recursively compute tree of particles and return coordinates of boxes
        boxes = Tuple{Float64,Array{Float64,1}}[]
        boxes_threaded = [Tuple{Float64,Array{Float64,1}}[] for i = 1:Threads.nthreads()]
        generate_tree_helper(particles, size, corner, return_boxes, boxes_threaded)
        for thread = eachindex(boxes_threaded)
           append!(boxes, boxes_threaded[thread])
        end
        boxes
    else
        # Recursively compute tree of particles and return that tree
        generate_tree_helper(particles, size, corner, false, Array{Tuple{Float64,Array{Float64,1}},1}[])
    end
end

"""
    generate_tree_helper(particles, size, corner, save_boxes, boxes_threaded)

Generate subtree of Barnes Hut octree given a list of particles and the corner
size of a cubical region of space. Optionally save coordinates of this region in
the array `boxes_threaded`. Helper function for `generate_tree`.
"""
function generate_tree_helper(particles::Array{Particle,1}, size::Float64, corner::Array{Float64,1},
                              save_boxes::Bool, boxes_threaded::Array{Array{Tuple{Float64,Array{Float64,1}},1},1})
    # If saving boxes, save current box coordinates
    if save_boxes
        push!(boxes_threaded[Threads.threadid()], (size, corner))
    end
    
    if length(particles) == 0
        # Base case 1: if no particles in octant, return empty leaf node
        Node(size/2, 0, corner .+ size/2, [nothing for i = 1:8])
    elseif length(particles) == 1
        # Base case 2: if one particle in octant, return leaf node with that particle
        Node(size/2, particles[1].mass, particles[1].pos, [nothing for i = 1:8])
    else
        # Recursive case
        
        # Calculate total mass and center of mass
        total_mass = sum(particle.mass for particle in particles)
        center_of_mass = sum(particle.mass .* particle.pos for particle in particles) / total_mass
        
        # Sort particles into 8 octants (using threading)
        octants = [Particle[] for oct = 1:8]
        octants_threaded = [[Particle[] for oct = 1:8] for thread = 1:Threads.nthreads()]
        new_corners = vec([corner + [i, j, k] for i = (0, size/2), j = (0, size/2), k = (0, size/2)])
        @sync for particle in particles
            Threads.@spawn begin
                oct_num = 0
                for oct = 1:8
                    if all(new_corners[oct] .<= particle.pos .<= new_corners[oct] .+ size/2)
                        oct_num = oct
                    end
                end
                if oct_num != 0
                    push!(octants_threaded[Threads.threadid()][oct_num], particle)
                end
            end
        end
        for thread = eachindex(octants_threaded)
            for oct = 1:8
               append!(octants[oct], octants_threaded[thread][oct])
            end
        end
        
        # Calculate child nodes recursively (using threading)
        children::Array{Union{Node, Nothing},1} = fill(nothing, 8)
        @sync for oct = 1:8
            Threads.@spawn children[oct] = generate_tree_helper(octants[oct], size/2, new_corners[oct],
                                                                save_boxes, boxes_threaded)
        end
        
        # Return node
        Node(size, total_mass, center_of_mass, children)
    end
end


# Barnes Hut Simulation

"""
    simulate!(particles, steps; Δt = 0.1, θ = 0.5, acc_func = grav_acc)

Approximately solve N-body problem for an array of particles for number of time
steps into the future using the Barnes Hut method. Return an array of particle
arrays, where each particle array is the state of all particles at each time
step.

# Arguments
- `particles::Array{Particle,1}`: an array of particles.
- `steps::Int64`: the number of time steps to run the simulation.
- `Δt::Float64 = 0.1`: length of time step in seconds (keyword argument).
- `θ::Float64 = 0.5`: node size to particle distance threshold for Barnes Hut
  simulation (keyword argument).
- `acc_func(mass::Float64, r::Array{Float64,1}) = grav_acc`: calculates particle
  acceleration from mass and distance vector.
"""
function simulate!(particles::Array{Particle,1}, steps::Int64;
                   Δt::Float64 = 0.1, θ::Float64 = 0.5, acc_func = grav_acc)
    # Initialize array of snapshots of particles at each time step
    frames = Array{Particle,1}[]
    
    # Store initial conditions as first frame
    push!(frames, deepcopy(particles))
    
    # Advance velocity by half time step to do leapfrog integration method
    vel_step!(particles, generate_tree(particles), Δt/2, θ, acc_func)
    
    # Simulate steps and save each frame
    for i = 1:steps
        step!(particles, generate_tree(particles), Δt, θ, acc_func)
        push!(frames, deepcopy(particles))
    end
    
    # Return frames
    frames
end

"""
    vel_step!(particles, tree, Δt, θ, acc_func)

Calculate the acceleration of each particle using `net_acc` and approximate new
velocity using `Δt`.
"""
function vel_step!(particles::Array{Particle,1}, tree::Node, Δt::Float64, θ::Float64, acc_func)
    # Calculate acceleration and approximately advance velocity (using threading)
    @sync for particle in particles
        Threads.@spawn begin
            particle.vel += net_acc(particle, tree, θ, acc_func) * Δt
        end
    end
end

"""
    step!(particles, tree, Δt, θ, acc_func)

Calculate the acceleration of each particle using `net_acc` and approximate new
velocity and position using `Δt`.
"""
function step!(particles::Array{Particle,1}, tree::Node, Δt::Float64, θ::Float64, acc_func)
    # Calculate acceleration, and approximately advance velocity and position (using threading)
    @sync for particle in particles
        Threads.@spawn begin
            particle.vel += net_acc(particle, tree, θ, acc_func) * Δt
            particle.pos += particle.vel * Δt
        end
    end
end

"""
    net_acc(particle, node, θ, acc_func)

Recursively calculate the net acceleration on a particle given the Barnes Hut
octree, a threshold `θ`, and the acceleration function `acc_func`.
"""
function net_acc(particle::Particle, node::Union{Node, Nothing}, θ::Float64, acc_func)
    # Recusively calculates net acceleration on a particle given a particle tree
    if node == nothing
        # Base case 1: if tree is empty, return 0 acceleration
        zeros(Float64, 3)
    else
        # Recursive case, Barnes Hut method
        s = node.size                            # Width of the cube that this node represents
        r = node.center_of_mass .- particle.pos  # Vector from particle to node center of mass
        if s/√sum(r.^2) < θ
            # Base case 2: if node size divided by distance is less than threshold θ,
            #              calculate gravitational acceleration
            acc_func(node.total_mass, r)
        else
            # Recursive case
            sum(net_acc(particle, child, θ, acc_func) for child in node.children)
        end
    end
end

"""
    grav_acc(mass, r; ϵ = 0.02)

Calculate the gravitation acceleration on a particle given the mass of another
object and the distance vector from particle to object. The usual gravitational
force is softened by a factor `ϵ` to prevent it from blowing up when the
distance is very small.
"""
function grav_acc(mass::Float64, r::Array{Float64,1}; ϵ::Float64 = 0.02)
    # Calculate gravitational acceleration, using node center of mass
    # and softening factor ϵ (to prevent blow up at d2 = 0)
    G::Float64 = 6.67430e-11
    ((G * mass) / ((sum(r.^2)+ϵ^2)^(3/2))) .* r
end


# Display and Animate Particles

"""
    show_particles(particles, θ = 30.0, ϕ = 30.0)

Using the `Plots` package, display an array of particles in 3D with azimuthal
camera angle θ and elevation angle ϕ.
"""
function show_particles(particles::Array{Particle,1}, θ::Float64 = 30.0, ϕ::Float64 = 30.0)
    scatter([Tuple(p.pos) for p = particles],
        lims = (-1, 1),
        camera = (θ, ϕ),
        size = (500, 500),
        label = "",
        background_color = :black,
        marker = (:circle, 3, 0.8, :white, stroke(0)))
end

"""
    show_quadview(particles)

Display four views of the particles in 3D at different camera angles.
"""
function show_quadview(particles::Array{Particle,1})
    scatter([Tuple(p.pos) for p = particles],
        lims = (-1, 1),
        camera = (30, 30),
        size = (800, 800),
        label = "",
        background_color = :black,
        marker = (:circle, 3, 0.8, :white, stroke(0)),
        layout = 4,
        subplot = 1)
    scatter!([Tuple(p.pos) for p = particles],
        lims = (-1, 1),
        camera = (0, 90),
        label = "",
        marker = (:circle, 3, 0.8, :white, stroke(0)),
        subplot = 2)
    scatter!([Tuple(p.pos) for p = particles],
        lims = (-1, 1),
        camera = (90, 0),
        label = "",
        marker = (:circle, 3, 0.8, :white, stroke(0)),
        subplot = 3)
    scatter!([Tuple(p.pos) for p = particles],
        lims = (-1, 1),
        camera = (60, 30),
        label = "",
        marker = (:circle, 3, 0.8, :white, stroke(0)),
        subplot = 4)
end

"""
    show_boxes(particles)

Calculate and display in 3D the cubes represented by the Barnes Hut octree.
"""
function show_boxes(particles::Array{Particle,1})
    p = scatter([Tuple(p.pos) for p = particles],
        camera = (30.0, 30.0),
        size = (500, 500),
        label = "",
        background_color = :black,
        marker = (:circle, 3, 0.8, :white, stroke(0)))
    boxes = generate_tree(particles, return_boxes = true)
    for (s, (x, y, z)) in boxes
        plot!([(x,y,z),(x+s,y,z),(x+s,y+s,z),(x,y+s,z),(x,y,z),(x,y,z+s),
               (x+s,y,z+s),(x+s,y+s,z+s),(x,y+s,z+s),(x,y,z+s),(x+s,y,z+s),
               (x+s,y,z),(x+s,y+s,z),(x+s,y+s,z+s),(x,y+s,z+s),(x,y+s,z)],
            linecolor = :red,
            label = "")
    end
    p
end

"""
    animate_frames(frames, frame_rate = 1, file_name = "animation"; θ = 30.0, ϕ = 30.0, quadview = false)

Plot particles for every `frame_rate` frames resulting from a Barnes Hut
simulation, and generate a gif file of this animation with a custom name.
Optionally, specify azimuthal camera angle θ and elevation angle ϕ, or set
quadview to `true` in order to animate at four different angles.
"""
function animate_frames(frames::Array{Array{Particle,1},1}, frame_rate::Int64 = 1, file_name = "animation";
                        θ::Float64 = 30.0, ϕ::Float64 = 30.0, quadview::Bool = false)
    anim = Animation()
    for i = 1:frame_rate:length(frames)
        if quadview
            show_quadview(frames[i])
        else
            show_particles(frames[i], θ, ϕ)
        end
        frame(anim)
    end
    gif(anim, "$file_name.gif")
end


# Generate Random Particles

"""
    rand_particles(num_particles::Int64)

Generate `num_particles` random particles using a distribution similar to a
disk-shaped galaxy.
"""
function rand_particles(num_particles::Int64)
    particles = Particle[]
    for i = 1:num_particles
        θ = 2π*rand()
        R = randexp()/ℯ
        z = (rand() > 0.5 ? 1 : -1) * (randexp()/(10ℯ))
        v = √(R*num_particles/1e6)
        mass = 1e4
        push!(particles, Particle([R*cos(θ), R*sin(θ), z], [v*sin(θ), -v*cos(θ), 0.0], mass))
    end
    particles
end

end # module