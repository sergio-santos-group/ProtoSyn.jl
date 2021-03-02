@doc """
    module Restraints

Contains restraint energy components, such as `bond_distance_restraint`.
"""
module Restraints

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, distance_matrix

    """
        Calculators.calc_bond_distance_restraint(::Any, pose::Pose, update_forces::Bool = false; x0::T = 2.0) where {T <: AbstractFloat}
        
    Calculate the pose bond distance restraint energy according to a quadratic
    potential. This potential only defines a maximum distance for the bond,
    based on the elements involved in the bond (See
    `ProtoSyn.max_bond_lengths`). If the pair of elements is not found, a
    default `x0` value is used instead.

    # See also
    `get_default_bond_distance_restraint`

    # Examples
    ```jldoctest
    julia> Calculators.Restraints.calc_bond_distance_restraint(pose)
    (0.0, nothing)
    ```
    """
    function calc_bond_distance_restraint(::Any, pose::Pose, update_forces::Bool = false; x0::T = 2.0) where {T <: AbstractFloat}

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
                    e += (d - _x0) * (d - _x0)
                    if update_forces
                        v = collect(pose.state[atom].t .- pose.state[bond].t)
                        factor1 = 2 * (d - _x0)

                        # * Assumes atom IDs are synched with pose.state
                        f[:, atom.id] -= v * factor1
                        f[:, bond.id] += v * factor1
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

    calc_bond_distance_restraint(pose::Pose, update_forces::Bool = false; x0::T = 2.0) where {T <: AbstractFloat} = begin
        calc_bond_distance_restraint(ProtoSyn.acceleration.active, pose, update_forces, x0 = x0)
    end

    """
        get_default_bond_distance_restraint(;α::T = 1.0) where {T <: AbstractFloat}

    Return the default bond distance restraint `EnergyFunctionComponent`. `α`
    sets the component weight (on an `EnergyFunction`). This function employs
    `calc_bond_distance_restraint` and therefore only travels the bond graph in
    order to calculate each bond distance energy.

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.Restraints.get_default_bond_distance_restraint()
             Name : Bond_Distance_Restraint
       Weight (α) : 1.0
    Update forces : true
          Setings :
             :x0 => 2.0
    ```
    """
    function get_default_bond_distance_restraint(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Bond_Distance_Restraint",
            calc_bond_distance_restraint,
            Dict{Symbol, Any}(:x0 => 2.0),
            α,
            true)
    end

    # --------------------------------------------------------------------------
    # * Flat bottom restraint base function

    MaskMap = Opt{Union{ProtoSyn.Mask{<: ProtoSyn.AbstractContainer}, Matrix{<: AbstractFloat}, Function}}

    function calc_flat_bottom_restraint(pose::Pose, update_forces::Bool; d1::T = 0.0, d2::T = 0.0, d3::T = Inf, d4::T = Inf, selection::Opt{AbstractSelection} = nothing, mask::MaskMap = nothing) where {T <: AbstractFloat}
        fbr = ProtoSyn.Calculators.get_flat_bottom_potential(d1 = d1, d2 = d2, d3 = d3, d4 = d4)
        ProtoSyn.Calculators.apply_potential(ProtoSyn.acceleration.active, pose, fbr, mask, selection)
    end # function
end