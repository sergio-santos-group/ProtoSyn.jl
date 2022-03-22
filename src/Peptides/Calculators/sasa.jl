module SASA

    using ProtoSyn
    using ProtoSyn.Peptides
    using YAML
    using Statistics

    function calculate_max_sasa(selection::AbstractSelection, residue_list::Vector{Char}; grammar::LGrammar = Peptides.grammar, probe_radius::T = 1.4, n_points::Int = 100) where {T <: AbstractFloat}

        max_sasa = Dict{String, T}()
        for residue in residue_list
            res = string(residue)
            pose = Peptides.build(grammar, [res])
            sasa = ProtoSyn.Calculators.SASA.calc_sasa(pose, selection, false,
                probe_radius = probe_radius,
                n_points = n_points)[1]
            for entry in ProtoSyn.alt_residue_names[ProtoSyn.one_2_three[residue]]
                max_sasa[entry] = sasa
            end
        end

        return max_sasa
    end

    function calculate_average_sasa(pose::Pose, selection::AbstractSelection, max_sasas::Dict{String, T}; probe_radius::T = 1.4, n_points::Int = 100) where {T <: AbstractFloat}

        s, _, sasas = ProtoSyn.Calculators.SASA.calc_sasa(pose, selection, false,
            probe_radius = probe_radius,
            n_points = n_points)

        a_resi_ids = [a.container.id for a in selection(pose, gather = true)]
        
        exposures = Vector{T}()
        for residue in eachresidue(pose.graph)
            resi_atoms_indexes = findall(==(residue.id), a_resi_ids)
            resi_sasa          = sum(sasas[resi_atoms_indexes])
            exposure           = resi_sasa / max_sasas[residue.name]
            push!(exposures, exposure)
        end

        return mean(exposures)
    end

    const default_sidechain_max_sasa = Dict{String, ProtoSyn.Units.defaultFloat}(
        YAML.load_file(
            joinpath(ProtoSyn.Peptides.resource_dir, "Calculators", "aa-max-sasa.yml")))

    function get_default_sasa_energy(;α::T = 1.0) where {T <: AbstractFloat}
        return ProtoSyn.Calculators.EnergyFunctionComponent(
            "SASA_Solvation",
            ProtoSyn.Calculators.SASA.calc_sasa_energy,
            (SidechainSelection() & !as"H") | an"CA",
            Dict{Symbol, Any}(
                :probe_radius                 => 2.8,
                :n_points                     => 100,
                :hydrophobicity_map           => Peptides.doolitle_hydrophobicity,
                :max_sasas                    => default_sidechain_max_sasa,
                :Ω                            => 0.5
            ),
            α,
            false
        )
    end
end