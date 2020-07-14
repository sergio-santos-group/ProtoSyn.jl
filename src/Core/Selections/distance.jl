export DistanceSelection
mutable struct DistanceSelection{M, T} <: AbstractSelection
    is_exit_node::Bool
    distance::Number
    sele::AbstractSelection

    DistanceSelection(distance::Number, sele::S) where {S <: AbstractSelection} = begin
        new{Stateful, Atom}(true, distance, sele)
    end
end

state_mode_type(::DistanceSelection{M, T}) where {M, T} = M
selection_type(::DistanceSelection{M, T})  where {M, T} = T


# --- Select -------------------------------------------------------------------
function select(sele::DistanceSelection, container::AbstractContainer)
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
end

select(sele::DistanceSelection{Stateful, Atom}, container::AbstractContainer, state::State) = begin
    return select(sele, container)(state)
end

# --- Short Syntax -------------------------------------------------------------
Base.:(:)(l::Number, r::AbstractSelection) = DistanceSelection(l, r)


# --- Examples -----------------------------------------------------------------
# (2.1:rn"GLN")(pose)

# TODO:
# - DistanceSelection of Stateful inner selections, as such:
# DistanceSelection(1.0, DistanceSelection(1.0, rn"ALA"))
# In line 18, the selection returns a function.