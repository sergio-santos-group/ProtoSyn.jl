export CastSelection
mutable struct CastSelection{M} <: AbstractSelection
    is_exit_node::Bool
    op::Function
    sele::AbstractSelection
    T::Type{<: AbstractContainer}

    CastSelection{M}(op::Function, sele::AbstractSelection, ::Type{T}) where {T <: AbstractContainer, M <: AbstractStateMode} = begin
        new{M}(true, op, sele, T)
    end
end

state_mode_type(::CastSelection{M}) where {M} = M

# --- Syntax -------------------------------------------------------------------
Base.any(sele::AbstractSelection, ::Type{T2}) where {T2 <: AbstractContainer} = begin
    T1 = selection_type(sele)
    @assert T1 < T2 "'any' function casts to higher order AbstractContainers only"
    sele.is_exit_node = false
    CastSelection{state_mode_type(sele)}(any, sele, T2)
end


Base.all(sele::AbstractSelection, ::Type{T2}) where {T2 <: AbstractContainer} = begin
    T1 = selection_type(sele)
    @assert T1 < T2 "'any' function casts to higher order AbstractContainers only"
    sele.is_exit_node = false
    CastSelection{state_mode_type(sele)}(all, sele, T2)
end


export cast
function cast(sele::AbstractSelection, ::Type{T2}) where {T2 <: AbstractContainer}
    T1 = selection_type(sele)
    @assert T1 > T2 "'cast' function casts to lower order AbstractContainers only"
    sele.is_exit_node = false
    CastSelection{state_mode_type(sele)}(all, sele, T2)
end


# --- Select -------------------------------------------------------------------
function select(sele::CastSelection{Stateless}, container::AbstractContainer)
    mask = select(sele.sele, container)
    return promote(mask, sele.T, container, sele.op)
end

function select(sele::CastSelection{Stateful}, container::AbstractContainer)

    _selector = select(sele.sele, container)

    return function (state::State)
        mask = _selector(state)
        return promote(mask, sele.T, container, sele.op)
    end
end

# --- Examples -----------------------------------------------------------------
# sele = any(an"CA", Residue)
# sele = all(rn"ALA", Segment)
# sele = cast(rn"ALA", Atom)