abstract type AbstractSelection end

abstract type AbstractStateMode end

struct Stateful <: AbstractStateMode end
struct Stateless <: AbstractStateMode end

state_mode_rule(::Type{Stateless}, ::Type{Stateless}) = Stateless
state_mode_rule(::Type{Stateful},  ::Type{Stateless}) = Stateful
state_mode_rule(::Type{Stateless}, ::Type{Stateful})  = Stateful
state_mode_rule(::Type{Stateful},  ::Type{Stateful})  = Stateful


# -------- MASKS ---------------------------------------------------------------

mutable struct Mask{T <: AbstractContainer}
    content::BitVector
end
Mask{T}() where {T <: AbstractContainer} = Mask{T}(BitVector())
Mask{T}(n::Int) where {T <: AbstractContainer} = Mask{T}(falses(n))

selection_type(m::Mask{T}) where {T} = T

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

# Cases where they are of the same type
Base.:|(m1::Mask{T}, m2::Mask{T}) where {T <: AbstractContainer} = begin
    return Mask{T}(m1.content .| m2.content)
end

Base.:&(m1::Mask{T}, m2::Mask{T}) where {T <: AbstractContainer} = begin
    return Mask{T}(m1.content .& m2.content)
end

Base.:!(m1::Mask{T}) where {T <: AbstractContainer} = begin
    m1.content = .!m1.content
    return m1
end

# ----------------- BINARY SELECTION -------------------------------------------

export BinarySelection
mutable struct BinarySelection{LM, RM} <: AbstractSelection
    is_exit_node::Bool
    op::Function
    left::AbstractSelection
    right::AbstractSelection

    BinarySelection(op::Function, l::L, r::R) where {L <: AbstractSelection, R <: AbstractSelection} = begin
        
        ML = state_mode_type(l)
        MR = state_mode_type(r)

        new{ML, MR}(true, op, r, l)
    end
end

state_mode_type(s::BinarySelection{LM, RM}) where {LM, RM} = state_mode_rule(LM, RM)

Base.:&(l::AbstractSelection, r::AbstractSelection) = BinarySelection(&, l, r)
Base.:|(l::AbstractSelection, r::AbstractSelection) = BinarySelection(|, l, r)


function promote(m1::Mask{T1}, m2::Mask{T2}, container::AbstractContainer) where {T1, T2}
    type1 = selection_type(m1)
    type2 = selection_type(m2)
    if type1 > type2
        m1 = promote(m1, type2, container)
    elseif type2 > type1
        m2 = promote(m2, type1, container)
    end

    m1, m2
end


select(sele::BinarySelection{Stateless, Stateless}, container::AbstractContainer) = begin
    
    left_mask  = select(sele.left,  container)
    right_mask = select(sele.right, container)

    l_mask, r_mask = promote(left_mask, right_mask, container)
    
    return sele.op(l_mask, r_mask)
end


select(sele::BinarySelection{Stateful, Stateless}, container::AbstractContainer) = select(sele, container, sele.right, sele.left)
select(sele::BinarySelection{Stateless, Stateful}, container::AbstractContainer) = select(sele, container, sele.left, sele.right)


select(sele::BinarySelection{T1, T2}, container::AbstractContainer, sless::AbstractSelection, sful::AbstractSelection) where {T1, T2} = begin
    
    sless_mask  = select(sless,  container)
    sful_selector = select(sful, container)

    return function (state::State)
        sful_mask = sful_selector(state)
        slessmask, sfulmask = promote(sless_mask, sful_mask, container)
        return sele.op(slessmask, sfulmask)
    end
end


select(sele::BinarySelection{Stateful, Stateful}, container::AbstractContainer) = begin
    
    left_selector  = select(sele.left,  container)
    right_selector = select(sele.right, container)

    return function (state::State)
        lmask = left_selector(state)
        rmask = right_selector(state)
        l_mask, r_mask = promote(lmask, rmask, container)
        return sele.op(l_mask, r_mask)
    end
end

# ------------------------------------------------------------------------------

export FieldSelection
mutable struct FieldSelection{S, T} <: AbstractSelection
    is_exit_node::Bool
    pattern::Regex
    field::Symbol
    
    FieldSelection{T}(pattern::String, field::Symbol) where {T <: AbstractContainer} = begin
        new{Stateless, T}(true, Regex(pattern), field)
    end
end

export select
function select(sele::FieldSelection{Stateless, T}, container::AbstractContainer) where {T <: AbstractContainer}

    n_items = counter(T)(container)
    mask = Mask{T}(n_items)

    for item in iterator(T)(container)
        if occursin(sele.pattern, getproperty(item, sele.field))
            mask[item.index] = true
        end
    end
    return mask
end

# SHORT SYNTAX

export @sn_str, @rn_str, @an_str, @as_str
macro sn_str(p); FieldSelection{Segment}(p, :name); end
macro rn_str(p); FieldSelection{Residue}(p, :name); end
macro an_str(p); FieldSelection{Atom}(p, :name); end
macro as_str(p); FieldSelection{Atom}(p, :symbol); end

state_mode_type(s::FieldSelection{M, T}) where {M, T} = M
selection_type(s::FieldSelection{M, T}) where {M, T} = T

# ------------------------------------------------------------------------------
# Unary Selection

export UnarySelection
mutable struct UnarySelection{M} <: AbstractSelection
    is_exit_node::Bool
    op::Function
    sele::AbstractSelection

    UnarySelection(op::Function, sele::AbstractSelection) where {T <: AbstractContainer} = begin
        new{Stateless}(true, op, sele)
    end
end

Base.:!(sele::AbstractSelection) = begin
    sele.is_exit_node = false
    UnarySelection(!, sele)
end

function select(sele::UnarySelection, container::AbstractContainer)
    mask = select(sele.sele, container)
    mask = (sele.op)(mask)
end

state_mode_type(s::UnarySelection{M}) where {M} = M

# ------------------------------------------------------------------------------
# True selections

export TrueSelection
mutable struct TrueSelection{M, T} <: AbstractSelection
    is_exit_node::Bool

    TrueSelection{T}() where {T <: AbstractContainer} = begin
        new{Stateless, T}(true)
    end
end

select(sele::TrueSelection{Stateless, Atom}, container::AbstractContainer) = Mask{Atom}(trues(count_atoms(container)))
select(sele::TrueSelection{Stateless, Residue}, container::AbstractContainer) = Mask{Residue}(trues(count_residues(container)))
select(sele::TrueSelection{Stateless, Segment}, container::AbstractContainer) = Mask{Segment}(trues(count_segments(container)))

# SHORT SYNTAX ?

# NOT SELECTION (TO DO)


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

function select(sele::DistanceSelection{Stateful, Atom}, container::AbstractContainer)
    mask      = select(sele.sele, container)
    if selection_type(mask) != Atom
        mask  = promote(mask, Atom, container)
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

# SHORT SYNTAX
# Ex: (2.1:rn"GLN")(pose)
Base.:(:)(l::Number, r::AbstractSelection) = DistanceSelection(l, r)


# --- PROMOTION / DEMOTION OF MASKS --------------------------------------------

function promote(mask::Mask{T1}, ::Type{T2}, container::AbstractContainer, f::Function = any)::Mask where {T1 <: AbstractContainer, T2 <: AbstractContainer}

    new_mask = Mask{T2}(counter(T2)(container))

    if T1 > T2
        count = counter(T2)
        iter = iterator(T1)
        i = 1
        for (entry, item) in zip(mask, iter(container))
            n = count(item)
            new_mask[i:(i+n-1)] = entry
            i += n
        end
    else
        count = counter(T1)
        iter = iterator(T2)
        i = 1
        for (index, item) in enumerate(iter(container))
            n = count(item)
            new_mask[index] = f(view(mask.content, i:(i+n-1)))
            i += n
        end
    end

    return new_mask
end


# export promote, demote

# promote(mask::Mask{Atom},    t::Type{Residue}, container::AbstractContainer, f::Function = all) = promote(mask, t, count_residues, eachresidue, count_atoms,    container, f)
# promote(mask::Mask{Atom},    t::Type{Segment}, container::AbstractContainer, f::Function = all) = promote(mask, t, count_segments, eachsegment, count_atoms,    container, f)
# promote(mask::Mask{Residue}, t::Type{Segment}, container::AbstractContainer, f::Function = all) = promote(mask, t, count_segments, eachsegment, count_residues, container, f)

# function promote(mask::Mask{T1}, t::Type{<: AbstractContainer}, counter1::Function, iterator::Function,
#     counter2::Function, container::AbstractContainer,
#     f::Function = all)::Mask where {T1 <: AbstractContainer, T2 <: AbstractContainer}

#     new_mask = Mask{t}(falses(counter1(container)))
    
#     i = 1
#     for (index, item) in enumerate(iterator(container))
#         if f(mask[i:(i + counter2(item)) - 1])
#             new_mask[index] = true
#         end
#         i += counter2(item)
#     end
    
#     return new_mask
# end

# # ---

# demote(mask::Mask{Residue}, t::Type{Atom},    container::AbstractContainer) = demote(mask, t, eachresidue, count_atoms,    container)
# demote(mask::Mask{Segment}, t::Type{Atom},    container::AbstractContainer) = demote(mask, t, eachsegment, count_atoms,    container)
# demote(mask::Mask{Segment}, t::Type{Residue}, container::AbstractContainer) = demote(mask, t, eachsegment, count_residues, container)

# function demote(mask::Mask{T1}, t::Type{<: AbstractContainer}, iterator::Function, counter::Function,
#     container::AbstractContainer)::Mask where {T1 <: AbstractContainer}

#     n_items  = counter(container)
#     new_mask = Mask{t}(falses(n_items))

#     i = 1
#     for (entry, item) in zip(mask, iterator(container))
#         if entry
#             new_mask[i:(i + counter(item) - 1)] = true
#         end
#         i += counter(item)
#     end

#     return new_mask
# end

# ---

# Callable functor (RESOLVE FUNCTION)

function (sele::AbstractSelection)(container::AbstractContainer)
    return select(sele, container)
end

function (sele::AbstractSelection)(container::AbstractContainer, state::State)
    if state_mode_type(sele) == Stateful
        return select(sele, container)(state)
    else
        return select(sele, container)
    end
end

(sele::AbstractSelection)(pose::Pose) = sele(pose.graph, pose.state)


# --- USE CASE EXAMPLES

# function design_step(pose::Pose, mask::Mask)
#     for (index, mask_i) in enumerate(mask)
#         if mask_i
#             random_residue = generate_random_residue()
#             mutate(pose, index, random_residue)
#         end
#     end
# end


# """
#     This function calculates the mask every step, therefore accounting for
#     previous mutations.
# """
# function design_1(pose::Pose, sele::AbstractSelection, steps::Int = 1000)
#     for step in range(steps)
#         mask = select(sele, pose) # or mask = select(sele, pose.graph, pose.state)
#         design_step(pose, mask)
#     end
# end


# """
#     This function calculates the mask every step, therefore accounting for
#     previous mutations.
#     This function uses the same mask every step, regardless of mutations that
#     have occurred.
# """
# function design_2(pose::Pose, sele::AbstractSelection, steps::Int = 1000)
#     f = select(sele, pose.graph)
#     for step in range(steps)
#         mask = f(pose.state)
#         design_step(pose, mask)
#     end
# end

export print_selection
function print_selection(io::IOStream, pose::Pose{Topology}, mask::Mask{T}) where {T <: AbstractContainer}

    if selection_type(mask) != Atom
        mask  = demote(mask, Atom, pose.graph)
    end
    
    Base.write(io, "MODEL\n")
    for (atom_index, atom) in enumerate(eachatom(pose.graph))
        sti = pose.state[atom.index] # returns and AtomState instance

        # In this file, selected atoms will be displayed in red while
        # non-selected atoms will be displayed in blue
        atom_symbol = mask[atom_index] ? "O" : "N"

        s = @sprintf("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%24s\n",
            atom.index, atom_symbol,
            atom.container.name, atom.container.container.code,
            atom.container.id,
            sti.t[1], sti.t[2], sti.t[3],
            atom_symbol)
            Base.write(io, s)
    end

    for atom in eachatom(pose.graph)
        Base.write(io, @sprintf("CONECT%5d", atom.index))
       foreach(n->Base.write(io, @sprintf("%5d",n.index)), atom.bonds)
       Base.write(io,"\n")
    end
    Base.write(io, "ENDMDL")
end