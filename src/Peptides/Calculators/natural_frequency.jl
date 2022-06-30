using ProtoSyn
using ProtoSyn.Calculators: EnergyFunctionComponent

aa_frequencies_filename = joinpath(
    Peptides.resource_dir, "Calculators/aa-frequencies.yml")

default_aa_frequencies = begin
    _data = ProtoSyn.read_yml(aa_frequencies_filename)
    Dict{Char, Float64}(
        map((x) -> x[1], collect(keys(_data))) .=> collect(values(_data)))
end


"""
    calc_aa_frequency([::Type{A}], pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; [aa_frequency_map::Dict{Char, T} = default_aa_frequencies]) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}

Calculates an energy value based on the total number of [`Residue`](@ref)
instances in the given [`Pose`](@ref) of a single type. For each
[`Residue`](@ref) type, a corrresponding value should be provided in
`aa_frequency_map` (this method adds the inverse of that value as an energy
reward). By default, `aa_frequency_map` reflects the natural distribution of
aminoacids, according to Krick et al.
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4209132/). Different
`aa_frequency_map` instances can introduce bias towards certain types of
aminoacids in simulations. If provided, an `AbstractSelection` `selection`
limits the subset of considered [`Residue`](@ref) instances for aminoacid
frequency calculation (`selection` is promoted to [`Residue`](@ref) level using
the [`promote`](@ref ProtoSyn.promote) method). `update_forces` has no
effect and exists only in order to standardize calls between Calculators. An
optional parameter `Type{<: AbstractAccelerationType}` can be provided,
stating the acceleration type used to calculate this energetic contribution
(See [ProtoSyn acceleration types](@ref), if not provided defaults to
`ProtoSyn.acceleration.active`).

# See also
[`get_default_aa_frequency`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.Calculators.calc_aa_frequency(pose, nothing, false)
(-456.50000000000017, nothing)
```
"""
function calc_aa_frequency(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; aa_frequency_map::Dict{Char, T} = default_aa_frequencies) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
    if selection === nothing
        sele = TrueSelection{Residue}()
    else
        sele = ProtoSyn.promote(selection, Residue)
    end

    residues = sele(pose, gather = true)

    e = 0.0
    for residue in residues
        name = ProtoSyn.ProtoSyn.three_2_one[residue.name]
        e -= aa_frequency_map[name]
    end

    return e, nothing
end # function

calc_aa_frequency(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; aa_frequency_map::Dict{Char, T} = default_aa_frequencies) where {T <: AbstractFloat} = begin
    calc_aa_frequency(ProtoSyn.acceleration.active, pose, selection, update_forces, aa_frequency_map = aa_frequency_map)
end


"""
    get_default_aa_frequency(;[α::T = 1.0]) where {T <: AbstractFloat}

Return the default natural frequency energy [`EnergyFunctionComponent`](@ref).
`α` sets the component weight (on an
[`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance, `1.0`
by default). This function employs [`calc_aa_frequency`](@ref) as the `:calc`
function.

# Settings
* `aa_frequency_map::Dict{Char, T}` - A mapping between [`Residue`](@ref) types (in 1-letter code format) and the natural frequency of aminoacids of that type;

# See also
[`calc_aa_ss_propensity`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.Calculators.get_default_aa_frequency()
🞧  Energy Function Component:
+---------------------------------------------------+
| Name           | Natural_AA_Freq                  |
| Alpha (α)      | 1.0                              |
| Update forces  | false                            |
| Calculator     | calc_aa_frequency                |
+---------------------------------------------------+
 |    +----------------------------------------------------------------------------------+
 ├──  ● Settings                      | Value                                            |
 |    +----------------------------------------------------------------------------------+
 |    | aa_frequency_map              | Dict{Char, Float64}(20 components)               |
 |    +----------------------------------------------------------------------------------+
 |    
 └──  ○  Selection: nothing
```
"""
function get_default_aa_frequency(;α::T = 1.0) where {T <: AbstractFloat}
    return EnergyFunctionComponent(
        "Natural_AA_Freq",
        calc_aa_frequency,
        nothing,
        Dict{Symbol, Any}(:aa_frequency_map => default_aa_frequencies),
        α,
        false)
end