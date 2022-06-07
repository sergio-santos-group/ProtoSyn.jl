using ProtoSyn

# Using the ProtoSyn.eval(:(...)), SidechainSelection actually becomes registered
# under the ProtoSyn module, as therefore works out of the box.
ProtoSyn.eval(:(begin

    export SidechainSelection

    """
        SidechainSelection()

    A [`SidechainSelection`](@ref) outputs a [`Mask`](@ref) of [`Atom`](@ref)
    instances where the sidechain atoms are marked as true (all atoms not named
    `C`, `H`, `CA`, `N` or `O` - marks caps as true).

    # State mode
    The state mode of [`SidechainSelection`](@ref) `M` is forced to be Stateless

    # Selection type

    The selection type of [`SidechainSelection`](@ref) `T` is forced to be [`Atom`](@ref).

    !!! ukw "Note:"
        This selection does not have a short syntax version.

    !!! ukw "Note:"
        This selection is provided by the Peptides module but registered under
        ProtoSyn Core module, in order to directly accessible and merged with
        other `AbstractSelection` instances. 

    # Examples
    ```jldoctest
    julia> SidechainSelection()
    SidechainSelection › (Atom)
    ```
    """
    mutable struct SidechainSelection{M, T} <: AbstractSelection
        SidechainSelection() = begin
            new{ProtoSyn.Stateless, Atom}()
        end
    end

    # --- Select -------------------------------------------------------------------
    function select(::SidechainSelection, container::ProtoSyn.AbstractContainer)

        n_atoms = count_atoms(container)
        mask = Mask{Atom}(n_atoms)

        for atom in eachatom(container)
            if !(atom.name in ["N", "H", "CA", "C", "O"])
                mask[atom.index] = true
            end
        end

        return mask
    end

    state_mode_type(::SidechainSelection{M, T}) where {M, T} = M
    selection_type(::SidechainSelection{M, T}) where {M, T} = T

    # --- Show -----------------------------------------------------------------
    function Base.show(io::IO, scs::SidechainSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
        lead = ProtoSyn.get_lead(level_code)
        if level_code === nothing
            level_code = LevelCode()
        end
        println(io, lead*"SidechainSelection › ($(selection_type(scs)))")
    end

end))