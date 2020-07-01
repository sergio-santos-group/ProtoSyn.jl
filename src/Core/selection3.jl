abstract type AbstractSelection end

abstract type AbstractStateMode end

struct Stateful <: AbstractStateMode end
struct Stateless <: AbstractStateMode end

state_mode_type(s::AbstractSelection) = typeof(s).parameters[1]
selection_type(s::AbstractSelection)  = typeof(s).parameters[2]

state_mode_rule(::Type{Stateless}, ::Type{Stateless}) = Stateless
state_mode_rule(::Type{Stateful},  ::Type{Stateless}) = Stateful
state_mode_rule(::Type{Stateless}, ::Type{Stateful})  = Stateful
state_mode_rule(::Type{Stateful},  ::Type{Stateful})  = Stateful



# ----

struct Mask{T <: AbstractContainer}
    content::BitVector
end
Mask{T}() where {T <: AbstractContainer} = Mask{T}(BitVector())

selection_type(m::Mask) = typeof(m).parameters[1]

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

# Cases where promotion and demotion need to occur
Base.:|(m1::Mask{T1}, m2::Mask{T2}) where {T1 <: AbstractContainer, T2 <: AbstractContainer} = begin
    return Mask{T}(m1.content .| m2.content)
end
# ---

# Suggestion, instead of left, right, call the Selections low_type, high_type
# since we know that relationship since the struct cosntructor

export BinarySelection
mutable struct BinarySelection{LM, HM, LT, HT} <: AbstractSelection
    is_exit_node::Bool
    op::Function
    low_sele::AbstractSelection
    high_sele::AbstractSelection

    BinarySelection(op::Function, l::L, r::R) where {L <: AbstractSelection, R <: AbstractSelection} = begin
        
        M1 = state_mode_type(l)
        M2 = state_mode_type(r)
        T1 = selection_type(l)
        T2 = selection_type(r)

        if T1 > T2
            new{M2, M1, T2, T1}(true, op, r, l)
        else
            new{M1, M2, T1, T2}(true, op, l, r)
        end
    end
end

state_mode_type(s::BinarySelection) = state_mode_rule(state_mode_type(s.high_sele), state_mode_type(s.low_sele))
selection_type(s::BinarySelection)  = min(selection_type(s.high_sele), selection_type(s.low_sele))

Base.:&(l::AbstractSelection, r::AbstractSelection) = BinarySelection(&, l, r)
Base.:|(l::AbstractSelection, r::AbstractSelection) = BinarySelection(|, l, r)


# Stateless cases

select(sele::BinarySelection{Stateless, Stateless, T, T}, container::AbstractContainer) where {T <: AbstractContainer} = begin
    sele.op(select(sele.high_sele, container), select(sele.low_sele, container))
end

select(sele::BinarySelection{Stateless, Stateless, LT, HT}, container::AbstractContainer) where {LT <: AbstractContainer, HT <: AbstractContainer} = begin
    sele.op(demote(select(sele.high_sele, container), LT, container), select(sele.low_sele, container))
end

# ---
# Stateless and Stateful cases
# This two functions are very similar, only changing wether they demote or not.
# Maybe this could be the job of another funtion?
select(sele::BinarySelection{LM, HM, T, T}, container::AbstractContainer) where {LM <: AbstractStateMode, HM <: AbstractStateMode, T <: AbstractContainer} = begin
    
    low_item  = select(sele.low_sele,  container)
    high_item = select(sele.high_sele, container)
    
    return function (state::State)
        LM === Stateful ? low_item  = low_item(state)  : nothing
        HM === Stateful ? high_item = high_item(state) : nothing
        return sele.op(high_item, low_item)
    end
end

select(sele::BinarySelection{LM, HM, LT, HT}, container::AbstractContainer) where {LM <: AbstractStateMode, HM <: AbstractStateMode, LT <: AbstractContainer, HT <: AbstractContainer} = begin
    
    low_item  = select(sele.low_sele,  container)
    high_item = select(sele.high_sele, container)
    
    return function (state::State)
        LM === Stateful ? low_item  = low_item(state)  : nothing
        HM === Stateful ? high_item = high_item(state) : nothing
        return sele.op(demote(high_item, LT, container), low_item)
    end
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
select(sele::ParameterSelection{Stateless, Atom},    container::AbstractContainer) = select(sele, container, eachatom,    count_atoms)
select(sele::ParameterSelection{Stateless, Residue}, container::AbstractContainer) = select(sele, container, eachresidue, count_residues)
select(sele::ParameterSelection{Stateless, Segment}, container::AbstractContainer) = select(sele, container, eachsegment, count_segments)

function select(sele::ParameterSelection{Stateless, T}, container::AbstractContainer, iterator::Function, counter::Function) where {T <: AbstractContainer}

    n_items = counter(container)
    mask = Mask{T}(falses(n_items))

    for item in iterator(container)
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

select(sele::TrueSelection{Stateless, Atom}, container::AbstractContainer) = Mask{Atom}(trues(count_atoms(container)))
select(sele::TrueSelection{Stateless, Residue}, container::AbstractContainer) = Mask{Residue}(trues(count_residues(container)))
select(sele::TrueSelection{Stateless, Segment}, container::AbstractContainer) = Mask{Segment}(trues(count_segments(container)))

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

function select(sele::DistanceSelection{Stateful, Atom}, container::AbstractContainer)
    mask      = select(sele.sele, container)
    if selection_type(mask) != Atom
        mask  = demote(mask, Atom, container)
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

export promote, demote

promote(mask::Mask{Atom},    t::Type{Residue}, container::AbstractContainer, f::Function = all) = promote(mask, t, count_residues, eachresidue, count_atoms,    container, f)
promote(mask::Mask{Atom},    t::Type{Segment}, container::AbstractContainer, f::Function = all) = promote(mask, t, count_segments, eachsegment, count_atoms,    container, f)
promote(mask::Mask{Residue}, t::Type{Segment}, container::AbstractContainer, f::Function = all) = promote(mask, t, count_segments, eachsegment, count_residues, container, f)

function promote(mask::Mask{T1}, t::Type{<: AbstractContainer}, counter1::Function, iterator::Function,
    counter2::Function, container::AbstractContainer,
    f::Function = all)::Mask where {T1 <: AbstractContainer, T2 <: AbstractContainer}

    new_mask = Mask{t}(falses(counter1(container)))
    
    i = 1
    for (index, item) in enumerate(iterator(container))
        if f(mask[i:(i + counter2(item)) - 1])
            new_mask[index] = true
        end
        i += counter2(item)
    end
    
    return new_mask
end

# ---

demote(mask::Mask{Residue}, t::Type{Atom},    container::AbstractContainer) = demote(mask, t, eachresidue, count_atoms,    container)
demote(mask::Mask{Segment}, t::Type{Atom},    container::AbstractContainer) = demote(mask, t, eachsegment, count_atoms,    container)
demote(mask::Mask{Segment}, t::Type{Residue}, container::AbstractContainer) = demote(mask, t, eachsegment, count_residues, container)

function demote(mask::Mask{T1}, t::Type{<: AbstractContainer}, iterator::Function, counter::Function,
    container::AbstractContainer)::Mask where {T1 <: AbstractContainer}

    new_mask = Mask{t}(falses(counter(container)))

    i = 1
    for (entry, item) in zip(mask, iterator(container))
        if entry
            new_mask[i:(i + counter(item) - 1)] = true
        end
        i += counter(item)
    end

    return new_mask
end

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