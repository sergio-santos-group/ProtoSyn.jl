using LinearAlgebra: I
using StaticArrays

#region NodeState --------------------------------------------------------------

export NodeState
mutable struct NodeState{T<:AbstractFloat}
    t::MVector{3,T}         # translation vector
    r::MMatrix{3,3,T,9}     # local to global rotation matrix
    
    # internal coordinates
    b::T                    # self<->parent bond length
    θ::T                    # self<->parent<->grandparent angle
    ϕ::T                    # self<->parent<->grandparent<->grand-grandparent dihedral
    
    Δϕ::T                   # change in dihedral angles (to be applied to children)
    changed::Bool           # flag
end

NodeState{T}() where {T} = begin
    NodeState{T}(
        MVector{3,T}(zeros(T, 3)),
        MMatrix{3,3,T,9}(Matrix{T}(I,3,3)),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        false
    )
end

Base.setproperty!(ns::NodeState{T}, key::Symbol, val) where T = begin
    if key == :Δϕ
        setfield!(ns, :changed, true)
    end
    setfield!(ns, key, val)
end

#endregion NodeState

#region State ------------------------------------------------------------------

export State
mutable struct State{T<:AbstractFloat}
    items::Vector{NodeState{T}}
    size::Int
    id::Int
    i2c::Bool
    c2i::Bool
    # x::Array{T,2}
    # f::Array{T,2}
    # v::Array{T,2}
    index_offset::Int
end
State{T}(n::Int) where T = begin
    items = Vector{NodeState{T}}(undef, n+3)
    for i=1:n+3
        items[i] = NodeState{T}()
    end
    items[1].t[1] = -1.0
    items[1].t[2] =  1.0
    items[2].t[1] = -1.0
    State{T}(items, n, -1, false, false, 3)
    # State{T}(items, n, -1, false, false, zeros(T,3,n))
end
State(::Type{T}, n::Int) where T = State{T}(n)
State(n::Int) = State{Float64}(n)
State(items::Vector{NodeState{T}}) where T = begin
    s = State{T}(items, length(items), -1, false, false, 0)
end
State{T}() where T = State{T}(0)

Base.getindex(s::State, i::Int) = begin
    s.items[i+s.index_offset]
end
Base.getindex(s::State, at::Atom) = begin
    s.items[at.index+s.index_offset]
end

Base.firstindex(s::State) = 1-s.index_offset
Base.lastindex(s::State) = s.size

Base.eltype(::Type{State{T}}) where T = T

Base.splice!(s::State, range::UnitRange{Int}) = begin
    i = s.index_offset + range.start
    j = s.index_offset + range.stop
    s2 = State(splice!(s.items, i:j))
    s.size -= s2.size
    s2
end

Base.append!(s1::State{T}, s2::State{T}) where T = begin
    for i = 1:s2.size
        push!(s1.items, s2[i])
    end
    s1.size += s2.size
    s1
end
Base.insert!(s1::State{T}, index::Integer, s2::State{T}) where T = begin
    index += s1.index_offset
    for i = 1:s2.size
        # push!(s1.items, s2[i])
        insert!(s1.items, index+i, s2[i])
        #println("s2 :", s2[i])
        #println(" s1:", s1[index+i])
        #println("  i:", index+i)
    end
    s1.size += s2.size
    s1
end




Base.copy(s::State{T}) where T = begin
    deepcopy(s)
end

export set!
set!(s::NodeState{T}, b::T, θ::T, ϕ::T) where T = begin
    s.b = b
    s.θ = θ
    s.ϕ = ϕ
    s
end

#endregion State