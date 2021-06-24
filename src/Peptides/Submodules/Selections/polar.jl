using ProtoSyn

# Using the ProtoSyn.eval(:(...)), PolarSelection actually becomes registered
# under the ProtoSyn module, as therefore works out of the box.
ProtoSyn.eval(:(begin

    export PolarSelection

    """
        PolarSelection()

    A [`PolarSelection`](@ref) outputs a [`Mask`](@ref) of [`Residue`](@ref)
    instances where the polar aminoacids are marked as true (as defined in
    `ProtoSyn.Peptides.polar_residues`).

    # State mode
    The state mode of [`PolarSelection`](@ref) `M` is forced to be Stateless

    # Selection type

    The selection type of [`PolarSelection`](@ref) `T` is forced to be [`Residue`](@ref).
    
    !!! ukw "Note:"
        This selection does not have a short syntax version.

    !!! ukw "Note:"
        This selection is provided by the Peptides module but registered under
        ProtoSyn Core module, in order to directly accessible and merged with
        other `AbstractSelection` instances.

    # Examples
    ```jldoctest
    julia> PolarSelection()
    PolarSelection › (Residue)
    ```
    """
    mutable struct PolarSelection{M, T} <: AbstractSelection
        PolarSelection() where {T <: ProtoSyn.AbstractContainer} = begin
            new{ProtoSyn.Stateless, Residue}()
        end
    end

    # --- Select -------------------------------------------------------------------
    function select(::PolarSelection, container::ProtoSyn.AbstractContainer)

        n_residues = count_residues(container)
        mask = Mask{Residue}(n_residues)

        for residue in eachresidue(container)
            if residue.name in Peptides.polar_residues
                mask[residue.index] = true
            end
        end
        return mask
    end

    state_mode_type(::PolarSelection{M, T}) where {M, T} = M
    selection_type(::PolarSelection{M, T}) where {M, T} = T

    # --- Show -----------------------------------------------------------------
    function Base.show(io::IO, ps::PolarSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
        lead = ProtoSyn.get_lead(level_code)
        if level_code === nothing
            level_code = LevelCode()
        end
        println(io, lead*"PolarSelection › ($(selection_type(ps)))")
    end
end))