@doc """
    module Restraints

Contains restraint energy components, such as `bond_distance_restraint`.
"""
module Restraints

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, distance_matrix

    """
        # TODO
    """
    function calc_bond_distance_restraint(::Any, pose::Pose; update_forces::Bool = false, x0::T = 2.0) where {T <: AbstractFloat}

        e  = 0.0
        if update_forces
            f = zeros(size(pose.state.f))
        end

        for atom in eachatom(pose.graph)
            for bond in atom.bonds
                d = ProtoSyn.distance(pose.state[atom], pose.state[bond])
                
                p = atom.symbol*bond.symbol
                if p in keys(ProtoSyn.Units.max_bond_lengths)
                    _x0 = ProtoSyn.Units.max_bond_lengths[p]
                else
                    println("Pair $p not found in max lengths list?")
                    _x0 = x0
                end

                if d > x0
                    e += (d - x0)
                    if update_forces
                        v = collect(pose.state[atom].t .- pose.state[bond].t)

                        # * Assumes atom IDs are synched with pose.state
                        f[:, atom.id] -= v
                        f[:, bond.id] += v
                    end
                end
            end
        end

        if update_forces
            return e, f
        else
            return e, nothing
        end
    end

    calc_bond_distance_restraint(pose::Pose; update_forces::Bool = false, x0::T = 2.0) where {T <: AbstractFloat} = begin
        calc_bond_distance_restraint(ProtoSyn.acceleration.active, pose, update_forces = update_forces, x0 = x0)
    end

    function get_default_bond_distance_restraint(α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Bond_Distance_Restraint",
            calc_bond_distance_restraint,
            Dict{Symbol, Any}(:x0 => 2.0),
            α,
            true)
    end

    # --------------------------------------------------------------------------

    MaskMap = Opt{Union{ProtoSyn.Mask{<: ProtoSyn.AbstractContainer}, Matrix{<: AbstractFloat}, Function}}

    function calc_flat_bottom_restraint(pose::Pose, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, selection::Opt{AbstractSelection} = nothing, mask::MaskMap = nothing) where {T <: AbstractFloat}
        fbr = ProtoSyn.Calculators.get_flat_bottom_potential(d1 = d1, d2 = d2, d3 = d3, d4 = d4)
        ProtoSyn.Calculators.apply_potential(ProtoSyn.acceleration.active, pose, fbr, mask, selection)
    end # function
end