using Distributions: DiscreteNonParametric, probs, support
using CurveFit: Polynomial
using Serialization

export SecondaryStructure
export SecondaryStructureTemplate
export DihedralTemplate

# Load Ramachandran samplers
phi_α_R_sampler  = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-α-R-sampler.jls"))
psi_α_R_sampler  = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-α-R-sampler.jls"))
phi_α_L_sampler  = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-α-L-sampler.jls"))
psi_α_L_sampler  = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-α-L-sampler.jls"))
phi_β_sampler    = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-β-sampler.jls"))
psi_β_sampler    = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-β-sampler.jls"))
phi_coil_sampler = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-coil-sampler.jls"))
psi_coil_sampler = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-coil-sampler.jls"))

# Load Ramachandran potentials
phi_α_R_potential  = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-α-R-potential.jls"))
psi_α_R_potential  = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-α-R-potential.jls"))
phi_α_L_potential  = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-α-L-potential.jls"))
psi_α_L_potential  = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-α-L-potential.jls"))
phi_β_potential    = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-β-potential.jls"))
psi_β_potential    = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-β-potential.jls"))
phi_coil_potential = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/phi-coil-potential.jls"))
psi_coil_potential = open(deserialize, joinpath(Peptides.resource_dir, "ramachandran/psi-coil-potential.jls"))


"""
    show_available_ramachandran_samplers()

Prints all available ramachandran samplers.

# See also
[`show_available_ramachandran_potentials()`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.show_available_ramachandran_samplers()
 • 8 ramachandran samplers available:
 ├── ProtoSyn.Peptides.phi_α_R_sampler
 ├── ProtoSyn.Peptides.psi_α_R_sampler
 ├── ProtoSyn.Peptides.phi_α_L_sampler
 ├── ProtoSyn.Peptides.psi_α_L_sampler
 ├── ProtoSyn.Peptides.phi_β_sampler
 ├── ProtoSyn.Peptides.psi_β_sampler
 ├── ProtoSyn.Peptides.phi_coil_sampler
 └── ProtoSyn.Peptides.psi_coil_sampler
```
"""
function show_available_ramachandran_samplers()
    samplers = ["phi_α_R_sampler", "psi_α_R_sampler", "phi_α_L_sampler", "psi_α_L_sampler", "phi_β_sampler", "psi_β_sampler", "phi_coil_sampler", "psi_coil_sampler"]
    println(" • $(length(samplers)) ramachandran samplers available:")
    for sampler in samplers[1:(end-1)]
        println(" ├── ProtoSyn.Peptides."*sampler)
    end
    println(" └── ProtoSyn.Peptides."*samplers[end])
end


"""
    show_available_ramachandran_potentials()

Prints all available ramachandran potentials.

# See also
[`show_available_ramachandran_samplers()`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.show_available_ramachandran_potentials()
 • 8 ramachandran potentials available:
 ├── ProtoSyn.Peptides.phi_α_R_potential
 ├── ProtoSyn.Peptides.psi_α_R_potential
 ├── ProtoSyn.Peptides.phi_α_L_potential
 ├── ProtoSyn.Peptides.psi_α_L_potential
 ├── ProtoSyn.Peptides.phi_β_potential
 ├── ProtoSyn.Peptides.psi_β_potential
 ├── ProtoSyn.Peptides.phi_coil_potential
 └── ProtoSyn.Peptides.psi_coil_potential
```
"""
function show_available_ramachandran_potentials()
    potentials = ["phi_α_R_potential", "psi_α_R_potential", "phi_α_L_potential", "psi_α_L_potential", "phi_β_potential", "psi_β_potential", "phi_coil_potential", "psi_coil_potential"]
    println(" • $(length(potentials)) ramachandran potentials available:")
    for potential in potentials[1:(end-1)]
        println(" ├── ProtoSyn.Peptides."*potential)
    end
    println(" └── ProtoSyn.Peptides."*potentials[end])
end


"""
    sample_ramachandran(dnp::DiscreteNonParametric; min_prob::T = 0.0) where {T <: AbstractFloat}

Sample a dihedral angle from a `DiscreteNonParametric` `dnp` distribution.
ProtoSyn makes available several ramachandran distributions (for the different
backbone dihedral angles & ramachandran zones of interest). Check
[`show_available_ramachandran_samplers`](@ref) for available options. `min_prob`
defines the minimum probability for the returned dihedral angle (higher
`min_prob` values result in more narrow distributions, closer to the ideal and
mean value). Returns the dihedral angle value in radians.

# Examples
```
julia> ProtoSyn.Peptides.sample_ramachandran(ProtoSyn.Peptides.psi_coil_sampler)
0.5934119456780721
```
"""
function sample_ramachandran(dnp::DiscreteNonParametric; min_prob::T = 0.0) where {T <: AbstractFloat}
    selected = probs(dnp) .> min_prob
    x = support(dnp)[selected]
    p = probs(dnp)[selected]
    @assert length(p) > 0 "The requested minimum probability ($min_prob) didn't yield any remaining dihedral angles. Try 'min_prob' value lower than $(maximum(probs(dnp)))."
    
    return rand(DiscreteNonParametric(x, p./(sum(p))))
end


"""
    DihedralTemplate{T}(angle::T, [sampler::Opt{DiscreteNonParametric} = nothing]) where {T <: AbstractFloat}

A [`DihedralTemplate`](@ref) offers a template for a given dihedral type, with a
mean value of `angle`. Optionally, this tempalte can also include a
`DiscreteNonParametric` distribution `sampler`, which returns angle values with
a given variation around the `value` ideal and mean value. This is used, for
example, in methods such as [`setss!`](@ref ProtoSyn.Peptides.setss!). Note that
for standardization purposes, methods in ProtoSyn expect the `angle` value to be
in radians. Optionally, this [`DihedralTemplate`](@ref) can have a name. This
serves no programatically purpose other than differentiation between
[`DihedralTemplate`](@ref) instances by the user.

# See also
[`SecondaryStructureTemplate`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.SecondaryStructure[:helix].ϕ
(Dihedral Template)    Phi (ϕ):   -1.117 rad |   -64.00 deg | Sampler: ✓ Set
```
"""
struct DihedralTemplate{T <: AbstractFloat}
    angle::T
    sampler::Opt{DiscreteNonParametric}
    name::Opt{String}
end

DihedralTemplate{T}(angle::T) where {T <: AbstractFloat} = DihedralTemplate{T}(angle, nothing, nothing)

Base.show(io::IO, dt::DihedralTemplate, level_code::Opt{LevelCode} = nothing) = begin
    init_lead       = ProtoSyn.get_lead(level_code)

    n = dt.name === nothing ? "" : dt.name
    sampler_set = dt.sampler === nothing ? "✖ Not set" : "✓ Set"
    s = @sprintf "(Dihedral Template) %10s: %8.3f rad | %8.2f deg | Sampler: %s" n dt.angle rad2deg(dt.angle) sampler_set
    println(io, init_lead*s)
end


"""
    SecondaryStructureTemplate{T}(ϕ::T, ψ::T, ω::T) where {T <: AbstractFloat}

Return a new [`SecondaryStructureTemplate`](@ref) with the given phi `ϕ`, psi
`ψ` and omega `ω` backbone angles (in radians) as [`DihedralTemplate`](@ref)
instances.

# See also
[`SecondaryStructure`](@ref ProtoSyn.Peptides.SecondaryStructure)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.SecondaryStructure[:helix]
Secondary Structure Template:
 ├── (Dihedral Template)    Phi (ϕ):   -1.117 rad |   -64.00 deg | Sampler: ✓ Set
 ├── (Dihedral Template)    Psi (ψ):   -0.820 rad |   -47.00 deg | Sampler: ✓ Set
 └── (Dihedral Template)  Omega (ω):    3.142 rad |   180.00 deg | Sampler: ✖ Not set
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
    Base.show(io, sst.ϕ, LevelCode([3]))
    Base.show(io, sst.ψ, LevelCode([3]))
    Base.show(io, sst.ω, LevelCode([4]))
end


"""
    SecondaryStructure

This constant holds default values for common
[`SecondaryStructureTemplate`](@ref) instances: `:helix`, `:linear`,
`:parallel_sheet` and `:antiparallel_sheet`.

# Examples
```jldoctest
julia> ProtoSyn.Peptides.SecondaryStructure
Dict{Symbol, ProtoSyn.Peptides.SecondaryStructureTemplate} with 5 entries:
  :antiparallel_sheet => Secondary Structure Template:…
  :linear             => Secondary Structure Template:…
  :left_handed_helix  => Secondary Structure Template:…
  :parallel_sheet     => Secondary Structure Template:…
  :helix              => Secondary Structure Template:…
```
"""
const SecondaryStructure = Dict{Symbol, SecondaryStructureTemplate}(
    # Values taken from https://proteopedia.org/ or measured from example structures

    :antiparallel_sheet => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}(-120.0°, phi_β_sampler,   "Phi (ϕ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 140.0°, psi_β_sampler,   "Psi (ψ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°,       nothing, "Omega (ω)")),

    :parallel_sheet     => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}(-119.0°, phi_β_sampler,   "Phi (ϕ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 113.0°, psi_β_sampler,   "Psi (ψ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°,       nothing, "Omega (ω)")),

    :linear             => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°, nothing,   "Phi (ϕ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°, nothing,   "Psi (ψ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°, nothing, "Omega (ω)")),

    :helix              => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( -64.0°, phi_α_R_sampler,   "Phi (ϕ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( -47.0°, psi_α_R_sampler,   "Psi (ψ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°,         nothing, "Omega (ω)")),

    :left_handed_helix  => SecondaryStructureTemplate(
        DihedralTemplate{ProtoSyn.Units.defaultFloat}(  57.0°, phi_α_L_sampler,   "Phi (ϕ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}(  47.0°, psi_α_L_sampler,   "Psi (ψ)"),
        DihedralTemplate{ProtoSyn.Units.defaultFloat}( 180.0°,         nothing, "Omega (ω)")))