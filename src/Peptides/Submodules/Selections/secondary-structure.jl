using ProtoSyn

# Using the ProtoSyn.eval(:(...)), PolarSelection actually becomes registered
# under the ProtoSyn module, as therefore works out of the box.
ProtoSyn.eval(:(begin

    using ProtoSyn
    using ProtoSyn.Units
    export SecondaryStructureSelection

    """
        SecondaryStructureSelection(ss::Symbol, [threshold::Float64 = 0.87])

    A [`SecondaryStructureSelection`](@ref) outputs a [`Mask`](@ref) of
    [`Residue`](@ref) instances where the residues are marked as true if the
    current `phi` and `psi` dihedrals are within the `threshold` angular
    distance (in radians, 50° by default) of the reference values for the
    requested `ss`
    [`SecondaryStructure`](@ref ProtoSyn.Peptides.SecondaryStructure) type.

    # State mode
    The state mode of [`SidechainSelection`](@ref) `M` is forced to be Stateful.

    # Selection type

    The selection type of [`SidechainSelection`](@ref) `T` is forced to be [`Residue`](@ref).

    # Short syntax
    * ss"helix" = Helix secondary structure
    * ss"parallel_sheet" = Parallel beta sheet secondary structure
    * ss"antiparallel_sheet" = Antiparallel beta sheet secondary structure
    * ss"linear" = Linear secondary structure

    !!! ukw "Note:"
        This selection is provided by the Peptides module but registered under
        ProtoSyn Core module, in order to directly accessible and merged with
        other `AbstractSelection` instances. 

    # Examples
    ```jldoctest
    julia> SecondaryStructureSelection(:parallel_sheet, deg2rad(20))
    SecondaryStructureSelection › parallel_sheet (± 20.0°)

    julia> ss"helix"
    SecondaryStructureSelection › helix (± 50.0°)
    ```
    """
    mutable struct SecondaryStructureSelection{M, T} <: AbstractSelection

        ss::Symbol
        threshold::Float64 # in radians

        SecondaryStructureSelection(ss::Symbol, threshold::Float64 = 50°) = begin
            new{ProtoSyn.Stateful, Residue}(ss, threshold)
        end
    end

    # --- Select -------------------------------------------------------------------
    function select(sele::SecondaryStructureSelection, container::ProtoSyn.AbstractContainer)

        return function (state::State)

            n_residues = count_residues(container)
            mask = Mask{Residue}(n_residues)

            for residue in eachresidue(container)

                # Case residue does not have phi or psi (first of last residue)
                phi_atom = Peptides.phi(residue)
                psi_atom = Peptides.psi(residue)
                phi_atom == nothing || psi_atom == nothing && continue
                
                # Measure current phi and psi dihedral values
                phi = ProtoSyn.getdihedral(state, phi_atom)
                psi = ProtoSyn.getdihedral(state, psi_atom)

                # Retrieve reference values for the query secondary structure
                if !(sele.ss in keys(ProtoSyn.Peptides.SecondaryStructure))
                    error("':$(sele.ss)' is not a valid SecondaryStructure key.")
                end
                ss = ProtoSyn.Peptides.SecondaryStructure[sele.ss]

                # Check phi
                phi_compliant = false
                if ss.ϕ.angle - sele.threshold < phi < ss.ϕ.angle + sele.threshold
                    phi_compliant = true
                end

                # Check psi
                psi_compliant = false
                if ss.ψ.angle - sele.threshold < psi < ss.ψ.angle + sele.threshold
                    psi_compliant = true
                end

                # Set mask entry
                if phi_compliant && psi_compliant
                    mask[residue.index] = true
                end
            end

            return mask
        end
    end

    state_mode_type(::SecondaryStructureSelection{M, T}) where {M, T} = M
    selection_type(::SecondaryStructureSelection{M, T}) where {M, T} = T

    # --- Short Syntax ---------------------------------------------------------

    export @ss_str
    macro ss_str(s)
        s = Symbol(s)
        if !(s in keys(ProtoSyn.Peptides.SecondaryStructure))
            error("'$(s)' is not a valid SecondaryStructure key.")
        end
        return SecondaryStructureSelection(s)
    end

    # --- Show -----------------------------------------------------------------
    function Base.show(io::IO, sss::SecondaryStructureSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
        lead = ProtoSyn.get_lead(level_code)
        if level_code === nothing
            level_code = LevelCode()
        end
        println(io, lead*"SecondaryStructureSelection › $(sss.ss) (± $(rad2deg(sss.threshold))°)")
    end

end))