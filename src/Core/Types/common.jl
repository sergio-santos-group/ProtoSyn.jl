"""
    ReadOnlyMatrix(m)

Return a read-only view into the parent array `m`. Most iteration, indexing,
abstract arrays and strided arrays interface functions defined for array `a`
are transparently defined for `ReadOnlyMatrix(a)`.
The exceptions are `setindex!` which is not allowed for `ReadOnlyMatrix`
and `similar` which uses the default definitions for `AbstractArray` from Base.
Also when used in bromdcasting `ReadOnlyMatrix` uses default bromdcast machinery.
Use `parent` function to access the parent array of `ReadOnlyMatrix`.

# Examples
```jldoctest
julia> m = [1 2; 3 4]
2×2 Array{Int64,2}:
 1  2
 3  4
julia> r = ReadOnlyMatrix(m)
2×2 ReadOnlyMatrix{Int64, 2, Array{Int64,2}}:
 1  2
 3  4
julia> r[1]
1
julia> r[1] = 10
ERROR: setindex! not defined for ReadOnlyMatrix{Int64, 2, Array{Int64,2}}
[...]
```
"""
struct ReadOnlyMatrix{T, N, P} <: AbstractArray{T, N}
    parent::P
    ReadOnlyMatrix(parent::AbstractArray{T, N}) where{T, N} =
        new{T, N, typeof(parent)}(parent)
end

Base.IteratorSize(::Type{<:ReadOnlyMatrix{T, N, P}}) where {T, N, P} =
    Base.IteratorSize(P)
Base.IteratorEltype(::Type{<:ReadOnlyMatrix{T, N, P}}) where {T, N, P} =
    Base.IteratorEltype(P)
Base.eltype(::Type{<:ReadOnlyMatrix{T, N, P}}) where {T, N, P} =
    eltype(P)
Base.size(rom::ReadOnlyMatrix, args...) = size(rom.parent, args...)
Base.elsize(rom::ReadOnlyMatrix, args...) = Base.elsize(rom.parent, args...)
Base.@propagate_inbounds Base.getindex(rom::ReadOnlyMatrix, I...) =
    getindex(rom.parent, I...)
Base.firstindex(rom::ReadOnlyMatrix) = firstindex(rom.parent)
Base.lastindex(rom::ReadOnlyMatrix) = lastindex(rom.parent)
Base.IndexStyle(::Type{<:ReadOnlyMatrix{T, N, P}}) where {T, N, P} = IndexStyle(P)
Base.iterate(rom::ReadOnlyMatrix, args...) = iterate(rom.parent, args...)
Base.length(rom::ReadOnlyMatrix) = length(rom.parent)

Base.axes(rom::ReadOnlyMatrix) = axes(rom.parent)
Base.strides(rom::ReadOnlyMatrix) = strides(rom.parent)
Base.unsafe_convert(p::Type{Ptr{T}}, rom::ReadOnlyMatrix) where {T} =
    Base.unsafe_convert(p, rom.parent)
Base.stride(rom::ReadOnlyMatrix, i::Int) = stride(rom.parent, i)
Base.parent(rom::ReadOnlyMatrix) = rom.parent