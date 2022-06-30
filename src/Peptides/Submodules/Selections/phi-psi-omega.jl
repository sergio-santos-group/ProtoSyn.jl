using ProtoSyn

# Using the ProtoSyn.eval(:(...)), PolarSelection actually becomes registered
# under the ProtoSyn module, as therefore works out of the box.
ProtoSyn.eval(:(begin

    export PhiSelection
    export PsiSelection
    export OmegaSelection

    """
        PhiSelection()

    A [`PhiSelection`](@ref) outputs a [`Mask`](@ref) of [`Atom`](@ref)
    instances where the phi-controlling atoms are marked as true.
    Phi-controlling atoms are non upstream-terminal [`Atom`](@ref) instances
    named "C". (make sure the IUPAC naming conventions are respected, consider
    using [`assign_default_atom_names!`](@ref ProtoSyn.Peptides.assign_default_atom_names!))
    for renaming peptide structures into the default atom names).

    # State mode
    The state mode of [`PhiSelection`](@ref) `M` is forced to be Stateless

    # Selection type

    The selection type of [`PhiSelection`](@ref) `T` is forced to be [`Atom`](@ref).

    !!! ukw "Note:"
        This selection is provided by the Peptides module but registered under
        ProtoSyn Core module, in order to directly accessible and merged with
        other `AbstractSelection` instances.

    # Examples
    ```jldoctest
    julia> PhiSelection()
    PhiSelection › (Atom)
    ```
    """
    mutable struct PhiSelection{M, T} <: AbstractSelection
        PhiSelection() = begin
            new{ProtoSyn.Stateless, Atom}()
        end
    end

    # --- Select ---------------------------------------------------------------
    function select(::PhiSelection, container::ProtoSyn.AbstractContainer)

        return (!UpstreamTerminalSelection{Residue}() & an"C")(container)
    end

    state_mode_type(::PhiSelection{M, T}) where {M, T} = M
    selection_type(::PhiSelection{M, T}) where {M, T} = T

    # --- Psi Selection --------------------------------------------------------

    """
        PsiSelection()

    A [`PsiSelection`](@ref) outputs a [`Mask`](@ref) of [`Atom`](@ref)
    instances where the psi-controlling atoms are marked as true.
    Psi-controlling atoms are non downstream-terminal [`Atom`](@ref) instances
    named "N". (make sure the IUPAC naming conventions are respected, consider
    using [`assign_default_atom_names!`](@ref ProtoSyn.Peptides.assign_default_atom_names!))
    for renaming peptide structures into the default atom names).

    # State mode
    The state mode of [`PsiSelection`](@ref) `M` is forced to be Stateless

    # Selection type

    The selection type of [`PsiSelection`](@ref) `T` is forced to be [`Atom`](@ref).

    !!! ukw "Note:"
        This selection is provided by the Peptides module but registered under
        ProtoSyn Core module, in order to directly accessible and merged with
        other `AbstractSelection` instances.

    # Examples
    ```jldoctest
    julia> PsiSelection()
    PsiSelection › (Atom)
    ```
    """
    mutable struct PsiSelection{M, T} <: AbstractSelection
        PsiSelection() = begin
            new{ProtoSyn.Stateless, Atom}()
        end
    end

    # --- Select ---------------------------------------------------------------
    function select(::PsiSelection, container::ProtoSyn.AbstractContainer)

        return (!DownstreamTerminalSelection{Residue}() & an"N")(container)
    end

    state_mode_type(::PsiSelection{M, T}) where {M, T} = M
    selection_type(::PsiSelection{M, T}) where {M, T} = T

    # --- Omega Selection --------------------------------------------------------

    """
        OmegaSelection()

    An [`OmegaSelection`](@ref) outputs a [`Mask`](@ref) of [`Atom`](@ref)
    instances where the omega-controlling atoms are marked as true.
    Omega-controlling atoms are non upstream-terminal [`Atom`](@ref) instances
    named "CA". (make sure the IUPAC naming conventions are respected, consider
    using [`assign_default_atom_names!`](@ref ProtoSyn.Peptides.assign_default_atom_names!))
    for renaming peptide structures into the default atom names).

    # State mode
    The state mode of [`OmegaSelection`](@ref) `M` is forced to be Stateless

    # Selection type

    The selection type of [`OmegaSelection`](@ref) `T` is forced to be [`Atom`](@ref).

    !!! ukw "Note:"
        This selection is provided by the Peptides module but registered under
        ProtoSyn Core module, in order to directly accessible and merged with
        other `AbstractSelection` instances.

    # Examples
    ```jldoctest
    julia> PsiSelection()
    PsiSelection › (Atom)
    ```
    """
    mutable struct OmegaSelection{M, T} <: AbstractSelection
        OmegaSelection() = begin
            new{ProtoSyn.Stateless, Atom}()
        end
    end

    # --- Select ---------------------------------------------------------------
    function select(::OmegaSelection, container::ProtoSyn.AbstractContainer)

        return (!UpstreamTerminalSelection{Residue}() & an"CA")(container)
    end

    state_mode_type(::OmegaSelection{M, T}) where {M, T} = M
    selection_type(::OmegaSelection{M, T}) where {M, T} = T

    # --- Show -----------------------------------------------------------------
    function Base.show(io::IO, sele::Union{PhiSelection, PsiSelection, OmegaSelection}, level_code::Opt{LevelCode} = nothing) where {M, T}
        lead = ProtoSyn.get_lead(level_code)
        if level_code === nothing
            level_code = LevelCode()
        end
        s = string(typeof(sele).name.name)
        println(io, lead*"$s › ($(selection_type(sele)))")
    end
end))