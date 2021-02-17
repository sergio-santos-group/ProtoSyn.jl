@doc """
    module Restraints

Contains restraint energy components, such as `bond_distance_restraint`.
"""
module Restraints

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, distance_matrix

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

    using LinearAlgebra

    function calc_clash_energy(A::Type{M}, pose::Pose; rmin::T = 3.0, update_forces::Bool = false) where {M <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        # coords must be in AoS format
        # This version doesn't calculate forces (too slow)
        
        dm1 = distance_matrix(A, pose, an"CA")
        if dm1 === nothing # No CA atoms were found
            return 0.0, nothing
        end
        dm2 = map((rij) -> (rij >= rmin ? 0.0 : - rij + rmin), Matrix(dm1))
        dm3 = tril(dm2, -1)
        
        return sum(dm3), nothing
    end

    calc_clash_energy(pose::Pose; rmin::T = 3.0, update_forces::Bool = false) where {T <: AbstractFloat}= begin
        calc_clash_energy(ProtoSyn.acceleration.active, pose; rmin = rmin, update_forces = update_forces)
    end

    clash_restraint = begin
        EnergyFunctionComponent("Clash_Restraint", calc_clash_energy)
    end
end