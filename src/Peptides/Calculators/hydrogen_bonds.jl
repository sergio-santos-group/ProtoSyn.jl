module HydrogenBonds

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent
    using ProtoSyn.Calculators.HydrogenBonds: HydrogenBondPair,
        HydrogenBondNetwork, calc_hydrogen_bond_network

    """
    # TODO
    """
    function generate_hydrogen_bond_network(pose::Pose)
        hbn = HydrogenBondNetwork()
        for residue in eachresidue(pose.graph)
            if (residue["H"] !== nothing) & (residue["N"] !== nothing)
                push!(hbn.donors, HydrogenBondPair(residue["H"], residue["N"]))
            end

            if (residue["O"] !== nothing) & (residue["C"] !== nothing)
                push!(hbn.acceptors, HydrogenBondPair(residue["O"], residue["C"]))
            end
        end

        return hbn
    end


    """
    # TODO
    """
    function fixate_hydrogen_bond_network!(efc::EnergyFunctionComponent, pose::Pose)
        efc.settings[:hydrogen_bond_network] = efc.settings[:hydrogen_bond_network](pose)
    end


    """
    # TODO
    """
    function get_default_hydrogen_bond_network_comp(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Hydrogen_Bond_Network",
            calc_hydrogen_bond_network,
            Dict{Symbol, Any}(
                :hydrogen_bond_network => generate_hydrogen_bond_network,
                :d0 => 3.0),
            α,
            false)
    end

end