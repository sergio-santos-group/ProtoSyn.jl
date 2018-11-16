@doc raw"""
    HarmonicBond(a1::Int64, a2::Int64, k::Float64, b0::Float64)

Harmonic Bond of the form

```math
E(r_{ab}) = \frac{1}{2}k_{ab}(r_{ab} - b_0)^2
```

where

```math
r_{ab} = \|\vec{r_{ab}}\| = \|\vec{r_b} - \vec{r_a}\|
```

# Arguments
- `a1::Int64, a2::Int64`: *global* atom indices.
- `k::Float64`: force constant (kJ mol⁻¹ nm⁻²).
- `b0::Float64`: equilibrium bond length (nm).

# Examples
```julia-repl
julia> Forcefield.HarmonicBond(1, 2, 2500, 0.19)
Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19)
```
See algo: [`Amber.evaluate!`](@ref Forcefield)
"""
mutable struct HarmonicBond
    
    a1::Int64
    a2::Int64
    k::Float64  # kJ mol⁻¹ nm⁻²
    b0::Float64 # nm

end
Base.show(io::IO, b::HarmonicBond) = print(io, "Forcefield.HarmonicBond(a1=$(b.a1), a2=$(b.a2), k=$(b.k), b0=$(b.b0))")

# ---------------------------------------------------------------------------------------------------------

@doc raw"""
    HarmonicAngle(a1::Int64, a2::Int64, a3::Int64, k::Float64, θ::Float64)

Harmonic Angle of the form

```math
E(θ_{abc})=\frac{1}{2}k_{abc}(\theta_{abc}-\theta)^{2}
```

# Arguments
- `a1::Int64, a2::Int64, a3::Int64`: *global* atom indices.
- `k::Float64`: force constant (kJ mol⁻¹ rad⁻²).
- `θ::Float64`: equilibrium angle value (rad).

# Examples
```julia-repl
julia> Forcefield.HarmonicAngle(1, 2, 3, 670.0, 1.92)
Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92)
```
See algo: [`Amber.evaluate!`](@ref Forcefield)
"""
mutable struct HarmonicAngle
    
    a1::Int64
    a2::Int64
    a3::Int64
    k::Float64 # kJ mol⁻¹ rad⁻²
    θ::Float64 # rad
    
end
Base.show(io::IO, b::HarmonicAngle) = print(io, "Forcefield.HarmonicAngle(a1=$(b.a1), a2=$(b.a2), a3=$(b.a3), k=$(b.k), θ=$(b.θ))")

# ---------------------------------------------------------------------------------------------------------

@doc raw"""
    DihedralCos(a1::Int64, a2::Int64, a3::Int64, a4::Int64, k::Float64, θ::Float64, mult::Float64)

Periodic Dihedral of the form

```math
E(\phi_{abcd})=K_{\phi}(1+cos(n\phi_{abcd}-\phi))
```

# Arguments
- `a1::Int64, a2::Int64, a3::Int64, a4::Int64`: *global* atom indices.
- `k::Float64`: force constant (kJ mol⁻¹).
- `θ::Float64`: equilibrium angle value (rad).
- `mult::Float64`: multiplicity.

# Examples
```julia-repl
julia> Forcefield.DihedralCos(1, 2, 3, 4, 10.46, 180.0, 2.0)
Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0)
```
See algo: [`Amber.evaluate!`](@ref Forcefield)
"""
mutable struct DihedralCos

    a1::Int64
    a2::Int64
    a3::Int64
    a4::Int64
    k::Float64    # kJ mol⁻¹
    θ::Float64    # rad
    mult::Float64 # multiplicity

end
Base.show(io::IO, b::DihedralCos) = print(io, "Forcefield.DihedralCos(a1=$(b.a1), a2=$(b.a2), a3=$(b.a3), a4=$(b.a4), k=$(b.k), θ=$(b.θ), mult=$(b.mult)")

# ----------------------------------------------------------------------------------------------------------

@doc raw"""
    Atom(name::String, σ::Float64, ϵ::Float64, q::Float64, excls::Array{Int64, 1}, pairs::Array{Int64, 1})

Define an atom. 
σ, ϵ and q describe the non-bonded interactions between atoms:

The Lennard-Jones interaction is in the form:
```math
E(r_{ab}) = 4ϵ_{ab}\left ( (\frac{σ_{ab}}{r_{ab}})^{12}-(\frac{σ_{ab}}{r_{ab}})^{6}\right )
```
where the Lorentz-Berthelot rule is applied. σ is the arithmetic average and ϵ is the geometric average:
```math
σ_{ab}=\frac{σ_{a}+σ_{b}}{2}
```
```math
ϵ_{ab}=\sqrt{(ϵ_{a}ϵ_{b})}
```
For this reason, σ and ϵ are applied here in the reduced form: ``\frac{σ}{2}`` and ``\sqrt{ϵ}``.

The Coulomb interation is in the form:
```math
E(r_{ab})=k_{ϵ}\frac{q_{a}q_{b}}{r_{ab}^{2}}
```
where
```math
k_{ϵ}=\frac{1}{4πϵ_{0}}=138.935485\,kJ\,nm\,mol⁻¹\,e⁻¹
```
For this reason, q is applied here in the reduced form: ``q\times \sqrt{k_{ϵ}}``

Exclusion list contains all atom indices who are excluded from non-bonded interactions
(i.e. are at 3 or less connections from this atom - includes pairs). Pair list contains atoms that are at
3 connections from this atom, and are involved in 1-4 interactions (and have a different combination rule
as a result).

# Arguments
- `name::String`: Atom name (example: "C", "H", etc).
- `σ::Float64`: finite distance at which the inter-particle potential is zero (nm).
- `ϵ::Float64`: depth of the potential well (kJ mol⁻¹).
- `q::Float64`: atom charge (eletron).
- `excls::Array{Int64, 1}`: exclusion list (as *global* atom indices).
- `pairs::Array{Int64, 1}`: pair list containing atoms that interfere in 1-4 interations (as *global* atom indices)

# Examples
```julia-repl
julia> Forcefield.Atom("N", 0.325, 0.711, 0.0017, [0, 1, 2, 3, 4, 5], [4, 5])
Forcefield.Atom(name="N", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5])
```
See algo: [`Amber.evaluate!`](@ref Forcefield)
"""
mutable struct Atom

    name::String
    σ::Float64 #nm
    ϵ::Float64 #kJ mol⁻¹
    q::Float64 #electron
    excls::Array{Int64, 1}
    pairs::Array{Int64, 1}
    
    Atom(name::String, σ::Float64, ϵ::Float64, q::Float64, excls::Array{Int64, 1}, pairs::Array{Int64, 1}) = new(name, σ/2, sqrt(ϵ), 11.787089759563214*q, excls, pairs)
end

Base.show(io::IO, b::Atom) = print(io, "Forcefield.Atom(name=$(b.name), σ=$(b.σ), ϵ=$(b.ϵ), q=$(b.q), excls=$(b.excls), pairs=$(b.pairs))")

# ----------------------------------------------------------------------------------------------------------

@doc raw"""
    Topology(atoms::Array{Atom}, bonds::Array{HarmonicBond}, angles::Array{HarmonicAngle}, dihedralsCos::Array{DihedralCos})

Gather all topology components.

# Arguments
- `atoms::Array{Atoms}`
- `bonds::Array{HarmonicBond}`
- `angles::Array{HarmonicAngle}`
- `dihedralsCos::Array{DihedralCos}`

# Examples
```julia-repl
julia> Forcefield.Forcefield(atoms, bonds, angles, dihedrals)
Forcefield.Topology(
 atoms=ProtoSyn.Forcefield.Atom[Forcefield.Atom(name="N", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5]), ...],
 bonds=ProtoSyn.Forcefield.HarmonicBond[Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19), ...],
 angles=ProtoSyn.Forcefield.HarmonicAngle[Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92), ...],
 dihedralsCos=ProtoSyn.Forcefield.DihedralCos[Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0), ...])
```
See also: [`Amber.load_from_json`](@ref Forcefield)
"""
mutable struct Topology

    atoms::Vector{Atom}
    bonds::Vector{HarmonicBond}
    angles::Vector{HarmonicAngle}
    dihedralsCos::Vector{DihedralCos}

end
Topology() = Topology(Vector{Atom}(), Vector{HarmonicBond}(), Vector{HarmonicAngle}(), Vector{DihedralCos}())
Base.show(io::IO, b::Topology) = print(io, "Forcefield.Amber.Topology(\n atoms=$(b.atoms),\n bonds=$(b.bonds),\n angles=$(b.angles),\n dihedralsCos=$(b.dihedralsCos))")

# ----------------------------------------------------------------------------------------------------------

# @doc raw"""
#     Energy(eBond::Float64, eAngle::Float64, eDihedral::Float64, eLJ::Float64, eLJ14::Float64, eCoulomb::Float64, eCoulomb14::Float64, eTotal::Float64)

# Energy components.

# # Examples
# ```julia-repl
# julia> Forcefield.Energy()
# Forcefield.Energy(eBond=0.0, eAngle=0.0, eDihedral=0.0, eLJ=0.0, eLJ14=0.0, eCoulomb=0.0, eCoulomb14=0.0, eTotal=0.0)

# julia> Forcefield.Energy(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 2.8)
# Forcefield.Energy(eBond=0.1, eAngle=0.2, eDihedral=0.3, eLJ=0.4, eLJ14=0.5, eCoulomb=0.6, eCoulomb14=0.7, eTotal=2.8)
# ```
# See also: [`Amber.evaluate!`](@ref Forcefield)
# """
# mutable struct Energy <: Common.AbstractEnergy
#     eBond::Float64
#     eAngle::Float64
#     eDihedral::Float64
#     eLJ::Float64
#     eLJ14::Float64
#     eCoulomb::Float64
#     eCoulomb14::Float64
#     eTotal::Float64
# end

# Energy() = Energy(zeros(Float64, 8)...)
# Base.show(io::IO, b::Energy) = print(io, "Forcefield.Energy(eBond=$(b.eBond), eAngle=$(b.eAngle), eDihedral=$(b.eDihedral), eLJ=$(b.eLJ), eLJ14=$(b.eLJ14), eCoulomb=$(b.eCoulomb), eCoulomb14=$(b.eCoulomb14), eTotal=$(b.eTotal))")