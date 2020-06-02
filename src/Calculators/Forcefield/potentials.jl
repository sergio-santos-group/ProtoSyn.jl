
struct CoulombLennardJonesPotential{K, N, T<:AbstractFloat} <: AbstractPotential
    key::KeyType{K,N}
    z::Int
    m::T
    q::T
    σ::T
    ϵ::T
end

struct HarmonicPotential{K, N, T<:AbstractFloat} <: AbstractPotential
    key::KeyType{K,N}
    kf::T
    x0::T
end

struct CosinePotential{K, N, T<:AbstractFloat} <: AbstractPotential
    key::KeyType{K,N}
    θ::T
    k::T
    n::T
end

struct NCosinePotential{K, N, T<:AbstractFloat} <: AbstractPotential
    key::KeyType{K,N}
    θ::Tuple{Vararg{T}}
    k::Tuple{Vararg{T}}
    n::Tuple{Vararg{T}}
end

const CoulombLennardJones{K,T} = CoulombLennardJonesPotential{K,1,T} where {K,T}
const    HarmonicBond{K,T} = HarmonicPotential{K,2,T} where {K,T}
const   HarmonicAngle{K,T} = HarmonicPotential{K,3,T} where {K,T}
const  CosineDihedral{K,T} =   CosinePotential{K,4,T} where {K,T}
const NCosineDihedral{K,T} =  NCosinePotential{K,4,T} where {K,T}

concrete(x::HarmonicPotential{K,N,T}, k::KeyType{W,N}) where {K,N,T,W}= begin
    HarmonicPotential{W,N,T}(k, x.kf, x.x0)
end

concrete(x::CoulombLennardJonesPotential{K,N,T}, k::KeyType{W,N}, q::T) where {K,N,T,W}= begin
    CoulombLennardJonesPotential{W,N,T}(k, x.z, x.m, q, x.σ, x.ϵ)
end
