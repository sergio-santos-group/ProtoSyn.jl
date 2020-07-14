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


# --- State Mode Rules --------------------------------------------------------- 
state_mode_rule(::Type{Stateless}, ::Type{Stateless}) = Stateless
state_mode_rule(::Type{Stateful},  ::Type{Stateless}) = Stateful
state_mode_rule(::Type{Stateless}, ::Type{Stateful})  = Stateful
state_mode_rule(::Type{Stateful},  ::Type{Stateful})  = Stateful

state_mode_type(::BinarySelection{LM, RM}) where {LM, RM} = state_mode_rule(LM, RM)


# --- Unary Operations ---------------------------------------------------------
Base.:&(l::AbstractSelection, r::AbstractSelection) = BinarySelection(&, l, r)
Base.:|(l::AbstractSelection, r::AbstractSelection) = BinarySelection(|, l, r)


# --- Select -------------------------------------------------------------------
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