using LinearAlgebra

"""
# TODO
"""
function add_semi_explicit_solvent!(::Type{A}, pose::Pose; selection::AbstractSelection = an"CA", probe_radius::T = 6.0, n_points::Int = 10, cut_off_radius::T = 1.85) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
    
    sync!(pose)
    
    dm = ProtoSyn.Calculators.full_distance_matrix(A, pose, selection)
    if A === ProtoSyn.CUDA_2
        dm = collect(dm)
    end
    
    atoms = selection(pose, gather = true)
    probe_radius_2 = 2 * probe_radius
    
    # Define sphere
    sphere = ProtoSyn.fibonacci_sphere(T, n_points)

    # Save pseudo atom positions to 'cloud'
    cloud = Vector{Vector{T}}()

    for i in 1:size(dm)[1]

        # define cloud i
        cloud_i  = copy(sphere)
        i_coords = pose.state[atoms[i]].t
        cloud_i  = map(x -> (x .* probe_radius) .+ i_coords, cloud_i)

        for j in 1:size(dm)[2]
            if (dm[i, j] > probe_radius_2) | (i === j)
                continue
            end

            # measure distances between atom j and each point of the cloud i
            j_coords = pose.state[atoms[j]].t
            for k in length(cloud_i):-1:1 # note the reverse loop; can't use n_points -> length(cloud_i) should be changing
                point_i = cloud_i[k]
                if norm(point_i - j_coords) < probe_radius
                    deleteat!(cloud_i, k)
                end
            end
        end

        for point in cloud_i
            push!(cloud, point)
        end
    end

    # Remove close together points in the cloud
    N = length(cloud)
    for i in N:-1:1
        point_i = cloud[i]
        for j in N:-1:i
            i == j && continue
            d = norm(point_i - cloud[j])
            if d < cut_off_radius
                N -= 1
                deleteat!(cloud, i)
                break
            end
        end
    end

    solvent = fragment(cloud)
    atoms   = collect(eachatom(solvent.graph))
    for (i, atom) in enumerate(atoms)
        atom.name = "SOL"
        atom.container.name = "SOL"
        atom.container.container.name = "SOL"

        if i > 1 # Because of re-index in the append_fragment_as_new_segment
            atom.parent = atoms[i-1]
        end
    end
    
    ProtoSyn.append_fragment_as_new_segment!(pose, solvent)

    # After re-index in the append_fragment_as_new_segment, fix the solvent root
    for atom in an"SOL"(pose, gather = true)[2:end]
        atom.parent = ProtoSyn.root(pose.graph)
        atom.container.parent = ProtoSyn.root(pose.graph).container
    end

    ProtoSyn.request_c2i!(pose.state)

    return pose
end

add_semi_explicit_solvent!(pose::Pose; selection::AbstractSelection = an"CA", probe_radius::T = 6.0, n_points::Int = 10, cut_off_radius::T = 1.85) where {T <: AbstractFloat} = begin
    add_semi_explicit_solvent!(ProtoSyn.acceleration.active, pose, selection = selection, probe_radius = probe_radius, n_points = n_points, cut_off_radius = cut_off_radius)
end


"""
# TODO
"""
add_semi_explicit_solvent(pose::Pose; selection::AbstractSelection = an"CA", probe_radius::T = 6.0, n_points::Int = 10, cut_off_radius::T = 1.85) where {T <: AbstractFloat} = begin
    pose_solvated = copy(pose)
    add_semi_explicit_solvent!(ProtoSyn.acceleration.active, pose_solvated, selection = selection, probe_radius = probe_radius, n_points = n_points, cut_off_radius = cut_off_radius)

    return pose_solvated
end


"""
# TODO
"""
function remove_semi_explicit_solvent!(pose::Pose)

    sync!(pose)
    residues = collect(eachresidue(pose.graph))
    for i in length(residues):-1:1
        if residues[i].name == "SOL"
            ProtoSyn.pop_residue!(pose, residues[i])
        end
    end

    return pose
end