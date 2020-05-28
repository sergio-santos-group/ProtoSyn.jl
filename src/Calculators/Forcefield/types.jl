
# abstract type AbstractFFComponentType{N} end
# # abstract type AbstractFFComponentType end
# abstract type AbstractFFComponent end

# const FFTypeKey=Tuple{Vararg{String}}

# parmcount(::Type{<:AbstractFFComponentType{N}}) where N = N
# struct CoulombLennardJones{T<:AbstractFloat} <: AbstractFFComponentType{5}
#     key::FFTypeKey
#     z::Int
#     m::T
#     q::T
#     σ::T
#     ϵ::T
# end

# """
# kf(x-x0)^2
# """
# struct HarmonicType{T<:AbstractFloat} <: AbstractFFComponentType{2}
#     key::FFTypeKey
#     kf::T
#     x0::T
# end

# struct CosineType{T<:AbstractFloat} <: AbstractFFComponentType{3}
#     key::FFTypeKey
#     θ::T
#     k::T
#     n::T
# end

# struct NCosineType{T<:AbstractFloat} <: AbstractFFComponentType{3}
#     key::FFTypeKey
#     # θ::NTuple{N,T}
#     # k::NTuple{N,T}
#     # n::NTuple{N,T}
#     θ::Tuple{Vararg{T}}
#     k::Tuple{Vararg{T}}
#     n::Tuple{Vararg{T}}
# end


# struct Bond{T} <: AbstractFFComponent
#     a1::Int
#     a2::Int
#     type::T
# end
# struct Angle{T} <: AbstractFFComponent
#     a1::Int
#     a2::Int
#     a3::Int
#     type::T
# end
# struct Dihedral{T} <: AbstractFFComponent
#     a1::Int
#     a2::Int
#     a3::Int
#     a4::Int
#     type::T
# end


# # const HarmonicBond{T}  where {T<:AbstractFloat} =  Bond{HarmonicType{T}}
# # const HarmonicAngle{T} where {T<:AbstractFloat} = Angle{HarmonicType{T}}
# # const CosineDihedral{T} where {T<:AbstractFloat} = Dihedral{CosineType{T}}
# # const NCosineDihedral{T} where {T<:AbstractFloat} = Dihedral{NCosineType{6,T}}
# #parmrule(::HarmonicBond) = HarmonicType

# # struct ForcefieldComponent{N,T<:AbstractFloat}
# #     i::NTuple{N,Int}
# #     type::T
# # end
# # const HarmonicBond{T}  where {T<:AbstractFloat}   = ForcefieldComponent{2,  HarmonicType{T}}
# # const HarmonicAngle{T} where {T<:AbstractFloat}   = ForcefieldComponent{3,  HarmonicType{T}}
# # const CosineDihedral{T} where {T<:AbstractFloat}  = ForcefieldComponent{4,    CosineType{T}}
# # const NCosineDihedral{T} where {T<:AbstractFloat} = ForcefieldComponent{4, NCosineType{6,T}}

# # const Forcefield = Dict{DataType, Vector{<:AbstractFFComponent}}

# # Base.push!(ff::Forcefield, item::T) where {T <: AbstractFFComponent} = begin
# #     c = get!(ff, T) do
# #         Vector{T}()
# #     end
# #     push!(c, item)
# #     ff
# # end





abstract type AbstractComponentType end
const KeyType = Tuple{Vararg{T}} where T

struct CoulombLennardJonesType{K, T<:AbstractFloat} <: AbstractComponentType
    key::K
    z::Int
    m::T
    q::T
    σ::T
    ϵ::T
end
#CoulombLennardJones{K,T}(k::Tuple{K},args::Any...) where {K,T} = CoulombLennardJones(k[1], args...)

struct HarmonicType{K, T<:AbstractFloat} <: AbstractComponentType
    key::K
    kf::T
    x0::T
end

struct CosineType{K, T<:AbstractFloat} <: AbstractComponentType
    key::K
    θ::T
    k::T
    n::T
end
struct NCosineType{K, T<:AbstractFloat} <: AbstractComponentType
    key::K
    θ::Tuple{Vararg{T}}
    k::Tuple{Vararg{T}}
    n::Tuple{Vararg{T}}
end

# 1-particle
const CoulombLennardJones{K,T} = CoulombLennardJonesType{Tuple{K},T} where {K,T}

# 2-particles
const  HarmonicBond{K,T} = HarmonicType{NTuple{2,K},T} where {K,T}

# 3-particles
const HarmonicAngle{K,T} = HarmonicType{NTuple{3,K},T} where {K,T}

# 4-particles
const  CosineDihedral{K,T} =  CosineType{NTuple{4,K},T} where {K,T}
const NCosineDihedral{K,T} = NCosineType{NTuple{4,K},T} where {K,T}

# const CoulombLennardJones{K,T} = CoulombLennardJonesType{Tuple{K},T} where {K,T}
# const HarmonicBond{K,T} = HarmonicType{Tuple{K,K},T} where {K,T}
# const HarmonicAngle{K,T} = HarmonicType{Tuple{K,K,K},T} where {K,T}
# const CosineDihedral{K,T} = CosineType{Tuple{K,K,K,K},T} where {K,T}
# const NCosineDihedral{K,T} = NCosineType{Tuple{K,K,K,K},T} where {K,T}

# keylen(::Type{CoulombLennardJones{K,T}}) where {K,T} = 1
keylen(::Type{T}) where {T<:AbstractComponentType} = length(T.parameters[1].parameters)

# const HarmonicAngle{K,T} = HarmonicType{3,K,T} where {K,T}
# const CosineDihedral{K,T} = CosineType{4,K,T} where {K,T}
# const NCosineDihedral{K,T} = NCosineType{4,K,T} where {K,T}
# keylen(::Type{HarmonicBond{K,T}}) where {K,T} = keylen(K)
# keylen(::Type{CoulombLennardJones}) = 1
# keylen(::Type{T{K,T}}) where {NK,T} = N
# keylen(::Type{HarmonicAngle}) = 3
# keylen(::Type{CosineDihedral}) = 4
# keylen(::Type{NCosineDihedral}) = 4

# const CLJAtom{T} = CoulombLennardJonesType{Int,T} where T
# const CLJAtomType{T} = CoulombLennardJonesType{String,T} where T

# const HarmonicBond{T} = HarmonicType{2,Int,T} where T
# const HarmonicBondType{T} = HarmonicType{2,String,T} where T

# const HarmonicAngle{T} = HarmonicType{3,Int,T} where T
# const HarmonicAngleType{T} = HarmonicType{3,String,T} where T

# const CosineDihedral{T} = CosineType{4,Int,T} where T
# const CosineDihedralType{T} = CosineType{4,String,T} where T

# const NCosineDihedral{T} = NCosineType{4,Int,T} where T
# const NCosineDihedralType{T} = NCosineType{4,String,T} where T

# concrete(::Type{HarmonicBondType{T}}) where T = HarmonicBond{T}
# concrete(::Type{HarmonicAngleType{T}}) where T = HarmonicAngle{T}
# concrete(::Type{CosineDihedralType{T}}) where T = CosineDihedral{T}
# concrete(::Type{NCosineDihedralType{T}}) where T = NCosineDihedral{T}
# concrete(t::T) where T = concrete(T)

# keylen(::Type{CoulombLennardJonesType{K, T}}) where {K,T} = N
# keylen(::Type{HarmonicType{N, K, T}}) where {N,K,T} = N
# keylen(::Type{CosineType{N, K, T}}) where {N,K,T} = N
# keylen(::Type{NCosineType{N, K, T}}) where {N,K,T} = N

#keylen(t::T) where T = keylen(T)
