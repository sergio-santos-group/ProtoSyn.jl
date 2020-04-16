export State
mutable struct State{T <: AbstractFloat}
    size::Int
    coords::Array{T, 2}
    #forces::Opt{Array{T, 2}}
    #velocs::Opt{Array{T, 2}}
    #energy::Energy{T}
    capacity::Int
end

State{T}(size::Int) where {T<:AbstractFloat} = begin
    State(size, zeros(T, 3, size), size)
    # State(size, zeros(T, size, 3), size)
    #State(size, zeros(T, size, 3), nothing, nothing, Energy{T}())
end

Base.resize!(state::State, n::Int) = begin
    if n>state.capacity
        state.capacity = n
        coords = state.coords
        state.coords = zeros(eltype(coords), 3, n)
        copyto!(state.coords, 1, coords, 1, 3*state.size)
    end
    state
end
