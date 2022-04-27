module Electrostatics

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, VerletList, MaskMap
    
    # TODO: Documentation

    """
        assign_acc2_eem_charges_from_file!(pose::Pose, filename::String, [selection::Opt{AbstractSelection}])

    Open and load ACC2 charges from a file `filename`, apply them to the given
    [`Pose`](@ref) `pose`. If an `AbstractSelection` `selection` is provided,
    the loaded charge are only applied to the selected [`Atom`](@ref) instances
    (note that both the number of loaded charge values and number of selected
    atoms must match). For more information on ACC2 charges, check
    https://acc2.ncbr.muni.cz/

    # See also
    [`assign_default_charges!`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.Electrostatics.assign_acc2_eem_charges_from_file!(pose, "charges.txt", aid"1:30")
    30-element Vector{Float64}:
      0.26429
     -0.1422
     -0.17175
     (...)
    ```
    """
    function assign_acc2_eem_charges_from_file!(pose::Pose, filename::String, selection::Opt{AbstractSelection} = nothing)

        if selection === nothing
            sele = TrueSelection{Atom}()
        else
            sele = ProtoSyn.promote(selection, Atom)
        end

        data    = split(open(f->read(f, String), filename), "\n")[2]
        T       = eltype(pose.state)
        charges = map((value) -> parse(T, value), split(data))

        atoms = sele(pose, gather = true)
        @assert length(atoms) === length(charges) "Attempted to apply $(length(charges)) to $(length(atoms)) selected atoms."
        for (i, atom) in enumerate(atoms)
            pose.state[atom].δ = charges[i]
        end

        return charges
    end 
    

    """
        assign_default_charges!(pose::Pose, res_lib::LGrammar, [selection::Opt{AbstractSelection}]; [supress_warn::Bool = false])

    Assign default charges to [`Pose`](@ref) `pose` from the given
    [`LGrammar`](@ref) `res_lib` entry, by [`Atom`](@ref) name. If an
    `AbstractSelection` `selection` is provided, only apply charges to the
    selected [`Atom`](@ref) instances. For non-canonical aminoacids and ligands
    (any [`Residue`](@ref) without an entry on `ProtoSyn.three_2_one`
    dictionary) and any [`Residue`](@ref) whose template have different
    [`Atom`](@ref) names, a warning is shown. Set `supress_warn` to `true` to
    ignore these warnings (`false`, by default).

    !!! ukw "Note:"
        Consider setting default atom names (from the same [`LGrammar`](@ref)), for example, using the [`assign_default_atom_names!`](@ref ProtoSyn.Peptides.assign_default_atom_names!) method.

    # See also
    [`assign_acc2_eem_charges_from_file!`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.Electrostatics.assign_default_charges!(pose, ProtoSyn.Peptides.grammar)
    1143-element Vector{Float64}:
     -0.025115728872692304
     -0.025115728872692304
     -0.025115728872692304
     (...)
    ```
    """
    function assign_default_charges!(pose::Pose, res_lib::LGrammar, selection::Opt{AbstractSelection} = nothing; supress_warn::Bool = false)

        if selection === nothing
            sele = TrueSelection{Residue}()
        else
            sele = ProtoSyn.promote(selection, Residue)
        end

        residues = sele(pose, gather = true)

        for residue in residues

            # Get template residue by residue name
            if !(residue.name.content in keys(ProtoSyn.three_2_one))
                @warn "Attempted to set default charge on residue $(residue.name.content) but the residue was not found on the ProtoSyn.three_2_one dictionary. Make sure the necessary LGrammar is loaded. Skipping this residue."
                continue
            end
            name = string(ProtoSyn.three_2_one[residue.name.content])
            template = res_lib.variables[name]
            if isa(template, Tautomer)
                # If multiple templates exist for a particular residue type,
                # find the correct tautomer matching the strucure's graph
                tautomer = ProtoSyn.find_tautomer(template, residue)
                if tautomer === nothing
                    # In case we don't have a template, for example, in
                    # no-sidechain models, try to apply the first tautomer
                    template = template.list[1]
                else
                    template = tautomer
                end
            end
            
            # Apply the template atom's charge to each of the residue's atoms
            for atom in residue.items
                if template.graph[1][atom.name] === nothing
                    !supress_warn && @warn "Template doesn't have atom $atom"
                    continue
                end # if
                δ = template.state[template.graph[1][atom.name]].δ
                pose.state[atom].δ = δ
            end # for
        end # for

        # Adjust the overall charge of all atoms so that the resulting charge of
        # the structure is ≈ 0.0
        cc = sum(getproperty.(pose.state.items, :δ)) / length(pose.state.items)
        for (i, atom) in enumerate(pose.state.items)
            atom.δ    -= cc
        end # for

        return getproperty.(pose.state.items, :δ)
    end # function

    # ---

    function calc_coulomb(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing, potential::Function = (x) -> 0.0) where {A <: ProtoSyn.AbstractAccelerationType}
        e, f = ProtoSyn.Calculators.apply_potential(A, pose, potential, update_forces, vlist, selection, mask)
        if e === 0.0 && all(x -> x.δ === 0.0, pose.state.items[4:end])
            @warn "The calculated Coulomb energy is 0.0 and it seems the evaluated Pose does not have any assigned charges. Consider using the `ProtoSyn.Calculators.Electrostatics.assign_default_charges!` method."
        end # if

        return e, f
    end # function

    calc_coulomb(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing, potential::Function = (x) -> 0.0) = begin
        calc_coulomb(ProtoSyn.acceleration.active, pose, selection, update_forces; mask = mask, vlist = vlist, potential = potential)
    end

    function get_default_coulomb(;α::T = 1.0)::EnergyFunctionComponent where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Coulomb",
            calc_coulomb,
            nothing, 
            Dict{Symbol, Any}(
                :mask => ProtoSyn.Calculators.get_intra_residue_mask,
                :vlist => nothing,
                :potential => ProtoSyn.Calculators.get_bump_potential_charges(c = 0.0, r = 20.0)),
            α,
            true)
    end

end