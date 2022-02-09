module Electrostatics

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, VerletList, MaskMap
    
    # TODO: Documentation

    """
    # TODO
    """
    function assign_acc2_eem_charges_from_file!(pose::Pose, filename::String)
        data = split(open(f->read(f, String), filename), "\n")[2]
        T = eltype(pose.state)
        charges = map((value) -> parse(T, value), split(data))
        for (i, atom) in enumerate(eachatom(pose.graph))
            pose.state[atom].δ = charges[i]
        end

        return charges
    end 
    

    """
    # TODO
    """
    function assign_default_charges!(pose::Pose, res_lib::LGrammar; supress_warn::Bool = false)

        for residue in eachresidue(pose.graph)

            # Get template residue by residue name
            name = string(ProtoSyn.three_2_one[residue.name.content])          
            template = res_lib.variables[name]
            if isa(template, Tautomer)
                # If multiple templates exist for a particular residue type,
                # find the correct tautomer matching the strucure's graph
                template = ProtoSyn.find_tautomer(template, residue)
            end
            
            # Apply the template atom's charge to each of the residue's atoms
            for atom in residue.items
                if template.graph[1][atom.name] === nothing
                    !supress_warn && @warn "Template doesn't have atom $atom"
                    continue
                end
                δ = template.state[template.graph[1][atom.name]].δ
                pose.state[atom].δ = δ
            end
        end

        # Adjust the overall charge of all atoms so that the resulting charge of
        # the structure is ≈ 0.0
        cc = sum(getproperty.(pose.state.items, :δ)) / length(pose.state.items)
        for (i, atom) in enumerate(pose.state.items)
            atom.δ    -= cc
        end

        return getproperty.(pose.state.items, :δ)
    end

    # ---

    function calc_coulomb(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing) where {A <: ProtoSyn.AbstractAccelerationType}
        cp = ProtoSyn.Calculators.coulomb_potential
        e, f = ProtoSyn.Calculators.apply_potential(A, pose, cp, update_forces, vlist, selection, mask)
        return e, f
    end # function

    calc_coulomb(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing) = begin
        calc_coulomb(ProtoSyn.acceleration.active, pose, selection, update_forces; mask = mask, vlist = vlist)
    end

    function get_default_coulomb(;α::T = 1.0) where {T <: AbstractFloat}
        # _mask = ProtoSyn.Calculators.get_bonded_mask()
        # _mask = ProtoSyn.Calculators.get_intra_residue_mask(!SidechainSelection())
        return EnergyFunctionComponent(
            "Coulomb",
            calc_coulomb,
            # !SidechainSelection(),
            nothing, 
            Dict{Symbol, Any}(:mask => nothing, :vlist => nothing),
            α,
            true)
    end

end