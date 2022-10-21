module Electrostatics

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, VerletList, MaskMap
    
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

        # ASSUMES UNIQUE NAMES
        data       = ProtoSyn.load(filename)
        data_atoms = collect(eachatom(data.graph))

        atoms = sele(pose, gather = true)
        @assert length(atoms) === length(data_atoms) "Attempted to apply $(length(data_atoms)) to $(length(atoms)) selected atoms."
        charges = Vector{eltype(pose.state)}()
        for atom in atoms
            i = findfirst((x) -> x.name === atom.name, data_atoms)
            δ = data.state[i].δ
            pose.state[atom].δ = δ
            push!(charges, δ)
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
    [`Atom`](@ref) names, a warning is shown. Set `supress_warn` to an
    AbstractSelection (or a boolean `true`) to ignore these warnings for the
    selected [`Atom`](@ref) instances (or all atoms, if `true`, is set to
    `false`, by default).

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
    function assign_default_charges!(pose::Pose, res_lib::LGrammar, selection::Opt{AbstractSelection} = nothing; supress_warn::Union{AbstractSelection, Bool} = false)

        if selection === nothing
            sele = TrueSelection{Residue}()
        else
            sele = ProtoSyn.promote(selection, Residue)
        end

        if typeof(supress_warn) <: AbstractSelection
            supress_warn_sele = ProtoSyn.promote(supress_warn, Atom)(pose, gather = true)
        end

        residues = sele(pose, gather = true)
        k = keys(res_lib.variables)

        for residue in residues

            # Get template residue by residue name
            if !(residue.name.content in keys(ProtoSyn.three_2_one))
                @warn "Attempted to set default charge on residue $(residue.name.content) but the residue was not found on the ProtoSyn.three_2_one dictionary. Make sure the necessary LGrammar is loaded. Skipping this residue."
                continue
            end
            name = string(ProtoSyn.three_2_one[residue.name.content])
            if !(name in k)
                @warn "Possible NCAA or ligand identified: $residue\nProtoSyn.jl will skip charge attribution in this residue.\nCheck if this is the desired behaviour."
                continue
            end
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
                    if typeof(supress_warn) === Bool && !supress_warn
                        @warn "Template doesn't have atom $atom"
                    elseif typeof(supress_warn) <: AbstractSelection && !(atom in supress_warn_sele)
                        @warn "Template doesn't have atom $atom"
                    end
                    
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

    """
        calc_coulomb([::A], pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; [mask::MaskMap = nothing], [vlist::Opt{VerletList} = nothing], [potential::Function = (x) -> 0.0]) where {A}
        
    Calculate the [`Pose`](@ref) `pose` Coulomb energy according to the given
    `potential` function, based on the cartesian coordinates (make sure the
    [`Pose`](@ref) `pose` is synched, see [`sync!`](@ref)). By default, the
    `potential` returns 0.0 for each atom-pair. This function iterates over all
    [`Atom`](@ref) instances in the provided [`Pose`](@ref) `pose` (See
    [Counters and Iterators](@ref)), unless an `AbstractSelection` `selection`
    is provided, limiting the selected [`Atom`](@ref) instances. If the
    `update_forces` flag is set to `true` (`false`, by default), also return the
    calculated forces based on this potential. Note that this function assumes
    [`Atom`](@ref)`.id` entries are synched between the
    [Graph](@ref graph-types) and [State](@ref state-types) (See
    [Indexation](@ref core-graph-methods-indexation)). An optional parameter
    `Type{<: AbstractAccelerationType}` can be provided, stating the
    acceleration type used to calculate this energetic contribution (See
    [ProtoSyn acceleration types](@ref)). Uses `ProtoSyn.acceleration.active` by
    default. This function makes use of the
    [`apply_potential!`](@ref ProtoSyn.Calculators.apply_potential!) framework. As
    such, an optional `mask` and `VerletList` `vlist` can be provided to limit
    the calculation. Make sure the [`Pose`](@ref) `pose` has charges assigned
    (see [`assign_acc2_eem_charges_from_file!`](@ref) and
    [`assign_default_charges!`](@ref)).

    # See also
    [`get_default_coulomb`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Calculators.Electrostatics.calc_coulomb(pose, nothing, false)
    (0.0, [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0])

    julia> ProtoSyn.Calculators.Electrostatics.calc_coulomb(pose, nothing, false, mask = ProtoSyn.Calculators.get_intra_residue_mask, potential = ProtoSyn.Calculators.get_bump_potential_charges(c = 0.0, r = 20.0))
    (-2.6046789109428206, nothing)
    ```
    """
    function calc_coulomb(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing, potential::Opt{Function} = nothing) where {A <: ProtoSyn.AbstractAccelerationType}
        
        if potential === nothing
            potential = function default(d::T; v::Opt{Tuple{T, T, T}} = nothing, qi::T = 0.0, qj::T = 0.0) where {T <: AbstractFloat}
                return 0.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
            end
        end

        e, f = ProtoSyn.Calculators.apply_potential(A, pose, potential, update_forces, vlist, selection, mask)
        if e === 0.0 && all(x -> x.δ === 0.0, pose.state.items[4:end])
            @warn "The calculated Coulomb energy is 0.0 and it seems the evaluated Pose does not have any assigned charges. Consider using the `ProtoSyn.Calculators.Electrostatics.assign_default_charges!` method."
        end # if

        return e, f
    end # function

    calc_coulomb(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing, potential::Opt{Function} = nothing) = begin
        calc_coulomb(ProtoSyn.acceleration.active, pose, selection, update_forces; mask = mask, vlist = vlist, potential = potential)
    end


    """
        get_default_coulomb(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default Coulomb [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). `α` sets the
    component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance). This
    component employs the [`calc_coulomb`](@ref) method, therefore defining a
    [`Pose`](@ref) energy based on a given potential. By default, this 
    [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) uses the
    [`get_bump_potential_charges`](@ref ProtoSyn.Calculators.get_bump_potential_charges)
    potential, with an intra-residue mask (see
    [`get_intra_residue_mask`](@ref ProtoSyn.Calculators.get_intra_residue_mask)).

    # Settings
    * `mask::Function` - Defines which atom-pairs to mask out of the result;
    * `vlist::VerletList` - If defined, the [`apply_potential!`](@ref ProtoSyn.Calculators.apply_potential!) method will only calculate the given atom-pairs in the `VerletList`;
    * `potential::Function` - Which potential to apply to each atom-pair;

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.Electrostatics.get_default_coulomb()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | Coulomb                          |
    | Alpha (α)      | 1.0                              |
    | Update forces  | true                             |
    | Calculator     | calc_coulomb                     |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | potential                     | bump_potential_charges                           |
     |    | vlist                         | nothing                                          |
     |    | mask                          | get_intra_residue_mask                           |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ○  Selection: nothing
    ```
    """
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