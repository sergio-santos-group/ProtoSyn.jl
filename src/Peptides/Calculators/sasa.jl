module SASA

    using ProtoSyn
    using ProtoSyn.Peptides
    using YAML

    function calculate_sasa_reference_energies(selection::Opt{AbstractSelection}, residue_list::Vector{Char}; grammar::LGrammar = Peptides.grammar, probe_radius::T = 6.0, n_points::Int = 100, hydrophobicity_map::Dict{String, T} = Peptides.doolitle_hydrophobicity) where {T <: AbstractFloat}

        ref_energy = Dict{String, T}()
        for residue in residue_list
            res = string(residue)
            pose = Peptides.build(grammar, [res])
            ref_energy[res] = ProtoSyn.Calculators.SASA.calc_sasa_energy(pose, selection, false,
                probe_radius = probe_radius,
                n_points = n_points,
                hydrophobicity_map = hydrophobicity_map,
                reference_energies = Dict{String, T}(ProtoSyn.one_2_three[residue] => 0.0))[1]
        end

        return ref_energy
    end

    const default_sidechain_sasa_reference_energy = YAML.load_file(
        joinpath(ProtoSyn.Peptides.resource_dir, "Calculators", "aa-sasa-ref.yml"))

    function get_default_sasa_energy(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "SASA_Solvation",
            calc_sasa_energy,
            nothing,
            Dict{Symbol, Any}(
                :probe_radius                 => 6.0,
                :n_points                     => 100,
                :hydrophobicity_map           => Peptides.doolitle_hydrophobicity,
                :reference_energies           => default_sidechain_sasa_reference_energy
            ),
            α,
            false
        )
    end
end