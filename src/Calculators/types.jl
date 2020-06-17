
abstract type AbstractPotential end

const KeyType = Tuple{Vararg{T,N}} where {T,N}

const ExclusionList = Dict{Int, Vector{Int}}

const BondedList = Dict{Union{Int,Symbol}, Vector{Tuple{Vararg{ProtoSyn.Atom}}}}

keylen(::Type{T}) where {T<:AbstractPotential} = length(T.types[1].types)

genkey(t::T...) where {T<:AbstractString} = join(t,":")

