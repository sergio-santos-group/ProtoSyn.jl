using Distributions: DiscreteNonParametric, probs, support
using Serialization

export SecondaryStructure
export SecondaryStructureTemplate
export DihedralTemplate

# Load Ramachandran samplers
phi_α_R = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-α-R.jls"))
psi_α_R = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-α-R.jls"))
phi_α_L = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-α-L.jls"))
psi_α_L = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-α-L.jls"))
phi_β   = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-β.jls"))
psi_β   = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-β.jls"))

function sample(dnp::DiscreteNonParametric; min_prob::T = 0.0) where {T <: AbstractFloat}
    selected = probs(dnp) .> min_prob
    x = support(dnp)[selected]
    p = probs(dnp)[selected]
    @assert length(p) > 0 "The requested minimum probability ($min_prob) didn't yield any remaining dihedral angles. Try 'min_prob' value lower than $(maximum(probs(dnp)))."
    
    return rand(DiscreteNonParametric(x, p./(sum(p))))
end

"""
# TODO
"""
struct DihedralTemplate{T <: AbstractFloat}
    angle::T
    sampler::Opt{DiscreteNonParametric}
end

DihedralTemplate{T}(angle::T) where {T <: AbstractFloat} = DihedralTemplate{T}(angle, nothing)

"""
    SecondaryStructureTemplate{T}(ϕ::T, ψ::T, ω::T) where {T <: AbstractFloat}

Return a new [`SecondaryStructureTemplate`](@ref) with the given phi `ϕ`, psi
`ψ` and omega `ω` backbone angles (in radians).

# TODO
# * Add explanation for standard deviation  

# Examples
```jldoctest
julia> t = ProtoSyn.Peptides.SecondaryStructureTemplate(-60°, -45°, 180°)
Secondary Structure Template:
 └─ Phi (ϕ):  -1.047 rad | Psi (ψ):  -0.785 rad | Omega (ω):   3.142 rad
             -60.000 deg |          -45.000 deg |            180.000 deg
```
"""
struct SecondaryStructureTemplate{T <: AbstractFloat}
    ϕ::DihedralTemplate{T}
    ψ::DihedralTemplate{T}
    ω::DihedralTemplate{T}
end

SecondaryStructureTemplate{T}(ϕ::T, ψ::T, ω::T) where {T <: AbstractFloat} = begin
    SecondaryStructureTemplate{T}(
        DihedralTemplate{T}(ϕ),
        DihedralTemplate{T}(ψ),
        DihedralTemplate{T}(ω))
end

Base.show(io::IO, sst::SecondaryStructureTemplate) = begin
    println(io, "Secondary Structure Template:")
    println(io, @sprintf " └─ Phi (ϕ): %8.3f rad | Psi (ψ): %8.3f rad | Omega (ω): %8.3f rad" sst.ϕ.angle sst.ϕ.angle sst.ψ.angle)
    println(io, @sprintf "             %8.3f deg |          %8.3f deg |            %8.3f deg" rad2deg(sst.ϕ.angle) rad2deg(sst.ϕ.angle) rad2deg(sst.ψ.angle))
end

"""
    SecondaryStructure

This constant holds default values for common
[`SecondaryStructureTemplate`](@ref) instances: `:helix`, `:linear`,
`:parallel_sheet` and `:antiparallel_sheet`.

# TODO
# * Add explanation for standard deviation 

# Examples
```jldoctest
julia> t = ProtoSyn.Peptides.SecondaryStructure[:helix]
Secondary Structure Template:
 └─ Phi (ϕ):  -1.047 rad | Psi (ψ):  -0.785 rad | Omega (ω):   3.142 rad
             -60.000 deg |          -45.000 deg |            180.000 deg
```
"""
const SecondaryStructure = Dict{Symbol, SecondaryStructureTemplate}(
    # Values taken from https://proteopedia.org/ or measured from example structures

    :antiparallel_sheet => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}(-120.0°, phi_β),   # phi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 140.0°, psi_β),   # psi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°)),         # omega

    :parallel_sheet     => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}(-119.0°, phi_β),   # phi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 113.0°, psi_β),   # psi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°)),         # omega

    :linear             => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°),          # phi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°),          # psi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°)),         # omega

    :helix              => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( -64.0°, phi_α_R), # phi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( -47.0°, psi_α_R), # psi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°)),         # omega

    :left_handed_helix  => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}(  57.0°, phi_α_L), # phi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}(  47.0°, psi_α_L), # psi
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°)))         # omega