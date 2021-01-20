@doc """
    module Restraints

Contains restraint energy components, such as `bond_distance_restraint`.
"""
module Restraints

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent

    function calc_bond_distance_restraint(::Any, pose::Pose; update_forces::Bool = false)

        x0 = 2.0 # max bond distance (Should be changable by the user ?)
        e  = 0.0
        f  = zeros(size(pose.state.f))
        stack = copy(origin(pose.graph).children)

        while length(stack) > 0
            i = pop!(stack)
            stack = vcat(stack, copy(i.children))
            for j in i.children
                d = ProtoSyn.distance(pose.state[i], pose.state[j])
                if d > x0
                    e += (d - x0)
                    if update_forces
                        v = collect(pose.state[j].t .- pose.state[i].t)
                        f[:, i.id] -= v
                        f[:, j.id] += v
                    end
                end
            end
        end

        return e, f
    end

    calc_bond_distance_restraint(pose::Pose; update_forces::Bool = false) = begin
        calc_bond_distance_restraint(ProtoSyn.acceleration.active, pose, update_forces = update_forces)
    end

    bond_distance_restraint = begin
        EnergyFunctionComponent("Bond_Distance_Restraint", calc_bond_distance_restraint)
    end

    # --------------------------------------------------------------------------

    function get_calc_distance_restraint(contact_map::Dict{Tuple{Int, Int}, T}) where {T <: AbstractFloat}
        return function calc_distance_restraint(pose::Pose; update_forces::Bool = false)
            
            x0 = 8.0 # in Angstrom
            e  = 0.0
            f  = zeros(size(pose.state.f))

            for ((r1, r2), α) in contact_map
                a1 = pose.graph[1][r1]["CA"]
                a2 = pose.graph[1][r2]["CA"]
                d = ProtoSyn.distance(pose.state[a1], pose.state[a2])

                if d > x0
                    e += (d - x0) * α
                    if update_forces
                        v = collect(pose.state[a2].t .- pose.state[a1].t) .* α
                        f[:, a1.id] -= v
                        f[:, a2.id] += v
                    end
                end
            end

            return e, f
        end
    end

    function load_contact_map(pose::Pose, filename::String)

        T = eltype(pose.state)
        contact_map = Dict{Tuple{Int, Int}, T}()

        open(filename, "r") do map_file
            for line in eachline(map_file)
                elems = split(line)
                length(elems[1]) > 4 && continue
                elems[1] == "END" && continue
                α = parse(T, elems[5])
                contact_map[(parse(Int, elems[1]), parse(Int, elems[2]))] = α
            end
        end

        calc_distance_restraint = get_calc_distance_restraint(contact_map)

        EnergyFunctionComponent("Distance_Restraint", calc_distance_restraint)
    end
end