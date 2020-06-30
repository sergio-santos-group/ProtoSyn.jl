abstract type AbstractSelection end

# Callable functor

function (sele::AbstractSelection)(container::AbstractContainer)
    return select(sele, container)
end

function (sele::AbstractSelection)(container::AbstractContainer, state::State)
    if typeof(sele).parameters[1] == Stateful
        return select(sele, container)(state)
    else
        return select(sele, container)
    end
end

(sele::AbstractSelection)(pose::Pose) = sele(pose.graph, pose.state)


# ----

abstract type AbstractStateMode end

struct Stateful <: AbstractStateMode end
struct Stateless <: AbstractStateMode end

const iterator = Dict(Atom => eachatom, Residue => eachresidue, Segment => eachsegment)

# ---

struct Mask{T <: AbstractContainer}
    content::BitVector
end
Mask{T}() where {T <: AbstractContainer} = Mask{T}(BitVector())

Base.getindex(m::Mask{T}, i::Int) where {T <: AbstractContainer} = m.content[i]
Base.getindex(m::Mask{T}, i::UnitRange{<: Real}) where {T <: AbstractContainer} = m.content[i]

Base.setindex!(m::Mask{T}, v::Bool, i::Real) where {T <: AbstractContainer} = m.content[i] = v
Base.setindex!(m::Mask{T}, v::Bool, i::UnitRange{<: Real}) where {T <: AbstractContainer} = m.content[i] .= v

@inline Base.length(m::Mask{T}) where {T <: AbstractContainer} = length(m.content)

Base.iterate(m::Mask{T}, (s,)=(1,)) where {T <: AbstractContainer} = begin
    if s > length(m)
        return nothing
    end
    (m.content[s], (s+1,))
end

Base.:|(m1::Mask{T}, m2::Mask{T}) where {T <: AbstractContainer} = begin
    return Mask{T}(m1.content .| m2.content)
end

Base.:&(m1::Mask{T}, m2::Mask{T}) where {T <: AbstractContainer} = begin
    return Mask{T}(m1.content .& m2.content)
end

# ---

export ParameterSelection
mutable struct ParameterSelection{S, T} <: AbstractSelection
    is_exit_node::Bool
    pattern::Regex
    field::Symbol
    
    ParameterSelection{T}(pattern::String, field::Symbol) where {T <: AbstractContainer} = begin
        new{Stateless, T}(true, Regex(pattern), field)
    end
end

export select
function select(sele::ParameterSelection{S, T}, container::AbstractContainer) where {S <: AbstractStateMode, T <: AbstractContainer}
    _iterate = iterator[T]
    mask = Mask{T}(falses(_iterate(container).size[end]))
    for item in _iterate(container)
        if occursin(sele.pattern, getproperty(item, sele.field))
            mask[item.index] = true
        end
    end
    return mask
end

# SHORT SYNTAX

export @sn_str, @rn_str, @an_str, @as_str
macro sn_str(p); ParameterSelection{Segment}(p, :name); end
macro rn_str(p); ParameterSelection{Residue}(p, :name); end
macro an_str(p); ParameterSelection{   Atom}(p, :name); end
macro as_str(p); ParameterSelection{   Atom}(p, :symbol); end

# ---
# True selections

export TrueSelection
mutable struct TrueSelection{M, T} <: AbstractSelection
    is_exit_node::Bool

    TrueSelection{T}() where {T <: AbstractContainer} = begin
        new{Stateless, T}(true)
    end
end

function select(sele::TrueSelection{M, T}, container::AbstractContainer) where {M <: AbstractStateMode, T <: AbstractContainer}
    return Mask{T}(trues(iterator[T](container).size[end]))
end

# SHORT SYNTAX ?

# ---

export DistanceSelection
mutable struct DistanceSelection{M, T} <: AbstractSelection
    is_exit_node::Bool
    distance::Number
    sele::AbstractSelection

    DistanceSelection(distance::Number, sele::S) where {S <: AbstractSelection} = begin
        new{Stateful, Atom}(true, distance, sele)
    end
end

function select(sele::DistanceSelection{M, T}, container::AbstractContainer) where {M <: AbstractStateMode, T <: AbstractContainer}
    mask      = select(sele.sele, container)
    atom_mask = demote(mask, Atom, container)
    co_sq     = sele.distance * sele.distance

    return function (state::State)
        masks = Array{Mask{Atom}, 1}()
        for (atom_i_index, atom_selected) in enumerate(atom_mask)
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

select(sele::DistanceSelection{M, T}, container::AbstractContainer, state::State) where {M <: AbstractStateMode, T <: AbstractContainer} = begin
    return select(sele, container)(state)
end

# SHORT SYNTAX
# Ex: (2.1:rn"GLN")(pose)
Base.:(:)(l::Number, r::AbstractSelection) = DistanceSelection(l, r)


# --- PROMOTION / DEMOTION OF MASKS --------------------------------------------

export promote, demote

promote(mask::Mask{Atom}, ::Type{Residue}, container::AbstractContainer, f::Function = all) = begin
    new_mask = Mask{Residue}(falses(count_residues(container)))
    return promote(mask, new_mask, eachresidue, count_atoms, container, f)
end

promote(mask::Mask{Atom}, ::Type{Segment}, container::AbstractContainer, f::Function = all) = begin
    new_mask = Mask{Residue}(falses(count_residues(container)))
    return promote(mask, new_mask, eachsegment, count_atoms, container, f)
end

promote(mask::Mask{Residue}, ::Type{Segment}, container::AbstractContainer, f::Function = all) = begin
    new_mask = Mask{Residue}(falses(count_residues(container)))
    return promote(mask, new_mask, eachsegment, count_residues, container, f)
end

function promote(mask::Mask{T1}, new_mask::Mask{T2}, iterator::Function,
    counter::Function, container::AbstractContainer,
    f::Function = all)::Mask where {T1 <: AbstractContainer, T2 <: AbstractContainer}
    
    i = 1
    for (index, item) in enumerate(iterator(container))
        if f(mask[i:(i += counter(item)) - 1])
            new_mask[index] = true
        end
    end
    
    return new_mask
end

# ---

demote(mask::Mask{Residue}, ::Type{Atom}, container::AbstractContainer) = begin
    new_mask = Mask{Atom}(falses(count_atoms(container)))
    return demote(mask, new_mask, eachresidue, count_atoms, container)
end

demote(mask::Mask{Segment}, ::Type{Atom}, container::AbstractContainer) = begin
    new_mask = Mask{Atom}(falses(count_atoms(container)))
    return demote(mask, new_mask, eachsegment, count_atoms, container)
end

demote(mask::Mask{Segment}, ::Type{Residue}, container::AbstractContainer) = begin
    new_mask = Mask{Atom}(falses(count_atoms(container)))
    return demote(mask, new_mask, eachsegment, count_residues, container)
end

function demote(mask::Mask{T1}, new_mask::Mask{T2}, iterator::Function, counter::Function,
    container::AbstractContainer)::Mask where {T1 <: AbstractContainer, T2 <: AbstractContainer}

    i = 1
    for (entry, item) in zip(mask, iterator(container))
        if entry
            new_mask[i:(i += counter(item) - 1)] = true
        else
            i += counter(item)
        end
    end

    return new_mask
end



# function design_step(pose::Pose, mask::Mask)
#     residue_list = apply(pose.graph, mask) # Return list of Residue object pointers
#     for (index, mask_i) in enumerate(mask)
#         if mask_i
#             random_residue = generate_random_residue()
#             mutate(residue_list[index], random_residue)
#         end
#     end
# end


# """
#     This function uses the same mask every step, regardless of mutations that
#     have occurred.
# """
# function design(pose::Pose, mask::Mask, steps::Int = 1000)
#     for step in range(steps)
#         design_step(pose, mask)
#     end
# end


# """
#     This function calculates the mask every step, therefore accounting for
#     previous mutations.
# """
# function design(pose::Pose, sele::AbstractSelection, steps::Int = 1000)
#     for step in range(steps)
#         f = select(sele, pose.graph)
#         mask = f(pose.state)
#         design_step(pose, mask)
#     end
# end