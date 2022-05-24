using ProtoSyn

# Using the ProtoSyn.eval(:(...)), PolarSelection actually becomes registered
# under the ProtoSyn module, as therefore works out of the box.
ProtoSyn.eval(:(begin

    export ProteinSelection

    """
        ProteinSelection()

    A [`ProteinSelection`](@ref) outputs a [`Mask`](@ref) of [`Residue`](@ref)
    instances where the protein aminoacids are marked as true. A known aminoacid
    is defined as having an entry both in `ProtoSyn.three_2_one` dictionary and
    in the default `Peptides` [`LGrammar`](@ref).

    # State mode
    The state mode of [`ProteinSelection`](@ref) `M` is forced to be Stateless

    # Selection type

    The selection type of [`ProteinSelection`](@ref) `T` is forced to be [`Residue`](@ref).
    
    !!! ukw "Note:"
        This selection does not have a short syntax version.

    !!! ukw "Note:"
        This selection is provided by the Peptides module but registered under
        ProtoSyn Core module, in order to directly accessible and merged with
        other `AbstractSelection` instances.

    # Examples
    ```jldoctest
    julia> ProteinSelection()
    ProteinSelection › (Residue)
    ```
    """
    mutable struct ProteinSelection{M, T} <: AbstractSelection
        ProteinSelection() = begin
            new{ProtoSyn.Stateless, Residue}()
        end
    end

    # --- Select -------------------------------------------------------------------
    function select(::ProteinSelection, container::ProtoSyn.AbstractContainer)

        n_residues = count_residues(container)
        mask = Mask{Residue}(n_residues)

        known_aminoacids = [x[1] for x in collect(keys(Peptides.grammar.variables)) if x !== "B"]
        for residue in eachresidue(container)
            if (residue.name in keys(ProtoSyn.three_2_one)) && (ProtoSyn.three_2_one[residue.name] in known_aminoacids)
                mask[residue.index] = true
            end
        end
        return mask
    end

    state_mode_type(::ProteinSelection{M, T}) where {M, T} = M
    selection_type(::ProteinSelection{M, T}) where {M, T} = T

    # --- Show -----------------------------------------------------------------
    function Base.show(io::IO, ps::ProteinSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
        lead = ProtoSyn.get_lead(level_code)
        if level_code === nothing
            level_code = LevelCode()
        end
        println(io, lead*"ProteinSelection › ($(selection_type(ps)))")
    end
end))