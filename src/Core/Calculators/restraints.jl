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
end