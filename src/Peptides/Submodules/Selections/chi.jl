using ProtoSyn

# Using the ProtoSyn.eval(:(...)), SidechainSelection actually becomes registered
# under the ProtoSyn module, as therefore works out of the box.
ProtoSyn.eval(:(begin

    export ChiSelection
    # Note: FieldSelection is a LEAF selection.

    """
        ChiSelection(n::Union{Int, Vector{Int}})

    A [`ChiSelection`](@ref) outputs a [`Mask`](@ref) of [`Atom`](@ref)
    instances where the chi-controlling [`Atom`](@ref) instances of sidechain
    are marked as true (all atoms defined in `Peptides.chi_dict`, by
    name - make sure the IUPAC naming conventions are respected, consider using
    [`assign_default_atom_names!`](@ref ProtoSyn.Peptides.assign_default_atom_names!))
    for renaming peptide structures into the default atom names).
    this `AbstractSelection` selects the n-chi dihedral angles (chi1, chi2, chi3
    or chi4). Optionally, a `Vector{Int}` can be provided, in wich case multiple
    chi-angles can be selected simultaneously.

    # State mode
    The state mode of [`ChiSelection`](@ref) `M` is forced to be Stateless.

    # Selection type

    The selection type of [`ChiSelection`](@ref) `T` is forced to be [`Atom`](@ref).

    # Short syntax
    * chi"..." = Chi angle

    !!! ukw "Note:"
        Using the short version syntax, it's possible to select multiple chi angles by using the 'r' flasg, as such: `chi"1|2|3|4"r`.

    !!! ukw "Note:"
        This selection is provided by the Peptides module but registered under ProtoSyn Core module, in order to directly accessible and merged with other `AbstractSelection` instances. 

    # Examples
    ```jldoctest
    julia> ChiSelection(1)
    ChiSelection (chi-1) › (Atom)


    julia> ChiSelection([1, 2])
    ChiSelection (chi-1 and chi-2) › (Atom)


    julia> chi"1"
    ChiSelection (chi-1) › (Atom)


    julia> chi"1|2|3|4"r
    ChiSelection (chi-1, chi-2, chi-3 and chi-4) › (Atom)
    ```
    """
    mutable struct ChiSelection{M, T} <: AbstractSelection

        n::Union{Int, Vector{Int}}

        ChiSelection(n::Union{Int, Vector{Int}}) = begin
            
            # ? Should this be changed to accomodate exotic NCAAs?
            @assert all(n .< 5) "ChiSelection only accepts chi values between 1 and 4."
            @assert all(n .> 0) "ChiSelection only accepts chi values between 1 and 4."

            new{ProtoSyn.Stateless, Atom}(n)
        end
    end

    # --- Select -------------------------------------------------------------------
    function select(sele::ChiSelection, container::ProtoSyn.AbstractContainer)

        @assert typeof(container) > Atom "ChiSelection cannot be applied to a single Atom instance."

        n_atoms = count_atoms(container)
        mask    = Mask{Atom}(n_atoms)

        for (index, atom) in enumerate(eachatom(container))
            if !(atom.container.name in keys(Peptides.chi_dict))
                mask[index] = false
                continue
            end

            chi_list = Peptides.chi_dict[atom.container.name]
            length(chi_list) <= 1 && begin
                mask[index] = false
                continue
            end

            N = length(chi_list)
            chi_query = [i + 1 for i in sele.n if i < N]

            if atom.name in chi_list[chi_query]
                mask[index] = true
            end
        end

        return mask
    end

    state_mode_type(::ChiSelection{M, T}) where {M, T} = M
    selection_type(::ChiSelection{M, T}) where {M, T} = T

    # --- Short Syntax -------------------------------------------------------------
    export @chi_str
    function parse_flags(n, flags)::Union{Int, Vector{Int}}
        if isempty(flags)
            return parse(Int, n)
        elseif 'r' in flags[1]
            n = [parse(Int, i) for i in split(n, "|")]
        else
            @error "Unkown flag $flags at ChiSelection()."
        end
    end
    macro chi_str(n, flags...); ChiSelection(parse_flags(n, flags)); end

    # --- Show -----------------------------------------------------------------
    function Base.show(io::IO, cs::ChiSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
        lead = ProtoSyn.get_lead(level_code)
        if level_code === nothing
            level_code = LevelCode()
        end
        s = Base.join(["chi-$i" for i in cs.n], ", ", " and ")
        println(io, lead*"ChiSelection ($s) › ($(selection_type(cs)))")
    end

end))