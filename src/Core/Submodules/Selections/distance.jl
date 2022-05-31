export DistanceSelection
# Note: DistanceSelection is a BRANCH selection.

"""
    DistanceSelection(distance::Number, sele::S) where {S <: AbstractSelection}

A [`DistanceSelection`](@ref) takes an input selection `sele` and outputs a
[`Mask`](@ref) of [`Atom`](@ref) instances within the given `distance` (in
Ansgtrom Å) of the selected atoms from `sele`.

# State mode
The state mode of [`DistanceSelection`](@ref) `M` is forced to be Stateful

# Selection type

The selection type of [`DistanceSelection`](@ref) `T` is forced to be [`Atom`](@ref).

# Examples
```jldoctest
julia> sele = DistanceSelection(2.0, rn"ALA")
DistanceSelection ❯ Within 2.0 Å (Atom)
 └── FieldSelection › Residue.name = ALA

julia> 2.0:rn"ALA"
DistanceSelection ❯ Within 2.0 Å (Atom)
 └── FieldSelection › Residue.name = ALA
```
"""
mutable struct DistanceSelection{M, T} <: AbstractSelection
    distance::Number
    sele::AbstractSelection

    DistanceSelection(distance::Number, sele::S) where {S <: AbstractSelection} = begin
        new{Stateful, Atom}(distance, sele)
    end
end

state_mode_type(::DistanceSelection{M, T}) where {M, T} = M
selection_type(::DistanceSelection{M, T}) where {M, T} = T


# --- Select -------------------------------------------------------------------
function select(sele::DistanceSelection, container::AbstractContainer)
    # Case inner selection is Stateless
    if state_mode_type(sele.sele) == Stateless
        mask      = select(sele.sele, container)
        if selection_type(mask) != Atom
            mask  = promote(mask, Atom, container)
        end

        # Case the given selection evaluates to all falses
        if !any(mask)
            return function (state::State)
                return mask
            end
        end
        
        co_sq     = sele.distance * sele.distance
        return function (state::State)
            masks = Array{Mask{Atom}, 1}()
            for (atom_i_index, atom_selected) in enumerate(mask)
                if atom_selected
                    _mask   = Mask{Atom}(falses(count_atoms(container)))
                    atom_i  = state[atom_i_index].t
                    for (atom_j_index, atom_j) in enumerate(state.items[4:end])
                        d_sq = sum(@. (atom_i - atom_j.t)^2)
                        if d_sq <= co_sq
                            _mask[atom_j_index] = true
                        end
                    end
                    push!(masks, _mask)
                end
            end
            return reduce(|, masks)
        end
    # Case inner selection is Stateful
    else
        _selector = select(sele.sele, container)
    
        co_sq     = sele.distance * sele.distance
        return function (state::State)
            mask = _selector(state)
            if selection_type(mask) != Atom
                mask  = promote(mask, Atom, container)
            end

            # Case the given selection evaluates to all falses       
            if !any(mask)
                return mask
            end

            masks = Array{Mask{Atom}, 1}()
            for (atom_i_index, atom_selected) in enumerate(mask)
                if atom_selected
                    _mask   = Mask{Atom}(falses(count_atoms(container)))
                    atom_i  = state[atom_i_index].t
                    for (atom_j_index, atom_j) in enumerate(state.items[4:end])
                        d_sq = sum(@. (atom_i - atom_j.t)^2)
                        if d_sq <= co_sq
                            _mask[atom_j_index] = true
                        end
                    end
                    push!(masks, _mask)
                end
            end
            return reduce(|, masks)
        end
    end
end

select(sele::DistanceSelection, container::AbstractContainer, state::State) = begin
    return select(sele, container)(state)
end

# --- Short Syntax -------------------------------------------------------------

Base.:(:)(l::Number, r::AbstractSelection) = DistanceSelection(l, r)

# --- Show ---------------------------------------------------------------------
function Base.show(io::IO, ds::DistanceSelection{M, T}, level_code::Opt{LevelCode} = nothing) where {M, T}
    lead = ProtoSyn.get_lead(level_code)
    if level_code === nothing
        level_code = LevelCode()
    end
    println(io, lead*"DistanceSelection ❯ Within $(ds.distance) Å ($T)")
    Base.show(io, ds.sele, vcat(level_code, 4))
end