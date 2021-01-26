module Restraints

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent

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
                        f[:, a1.id] += v
                        f[:, a2.id] -= v
                    end
                end
            end

            return e, f
        end
    end

    function load_contact_map(::Type{T}, filename::String) where {T <: AbstractFloat}

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

        EnergyFunctionComponent("Contact_Map_Restraint", calc_distance_restraint)
    end

    load_contact_map(filename::String) = begin
        load_contact_map(ProtoSyn.Units.defaultFloat, filename)
    end

end