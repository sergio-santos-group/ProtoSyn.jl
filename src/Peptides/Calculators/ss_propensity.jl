using ProtoSyn
using ProtoSyn.Calculators: EnergyFunctionComponent

aa_ss_propensity_filename = joinpath(
    Peptides.resource_dir, "Calculators/aa-ss-propensity.yml")


"""
    load_default_aa_ss_propensity([::Type{T}]) where {T <: AbstractFloat}

Load the default aminoacid secondary structure propensity map from
`ProtoSyn.Peptides.aa_ss_propensity_filename`. Set the values to be of type `T`
(`ProtoSyn.Units.defaultFloat`, by default). ProtoSyn automatically loads this
map to `ProtoSyn.Peptides.Calculators.default_aa_ss_propensity` during loading.
The values are extracted from Contantini et al. work (See
https://www.sciencedirect.com/science/article/pii/S0006291X06002543.

# Examples
```
julia> ProtoSyn.Peptides.Calculators.load_default_aa_ss_propensity()
Dict{Char, Dict{Char, Float64}} with 3 entries:
  'H' => Dict('P'=>0.50, 'M'=>1.21, 'K'=>1.11, 'F'=>1.01…)
  'E' => Dict('P'=>0.44, 'M'=>0.99, 'K'=>0.83, 'F'=>1.43…)
  'C' => Dict('P'=>1.72, 'M'=>0.83, 'K'=>1.00, 'F'=>0.76…)
```
"""
function load_default_aa_ss_propensity(::Type{T}) where {T <: AbstractFloat}
    ss_propensities = Dict{Char, Dict{Char, T}}()
    _data = ProtoSyn.read_yml(aa_ss_propensity_filename)
    for (key, value) in _data
        _v = map((x) -> x[1], collect(keys(value))) .=> collect(values(value))
        ss_propensities[string(key)[1]] = Dict{Char, Float64}(_v)
    end
    
    return ss_propensities
end

load_default_aa_ss_propensity() = begin
    load_default_aa_ss_propensity(ProtoSyn.Units.defaultFloat)
end

default_aa_ss_propensity = load_default_aa_ss_propensity()


"""
    calc_aa_ss_propensity([::Type{A}], pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; [secondary_structure::Union{Function, String} = ProtoSyn.Peptides.categorize_ss_from_dihedral_angles], [aa_ss_propensity_map::Dict{Char, Dict{Char, T}} = default_aa_ss_propensity]) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}

Calculates the [`Pose`](@ref) `pose` energy in accordance to the
`aa_ss_propensity_map`. This map correlates aminoacid nature to secondary
structure propensity. The [`Pose`](@ref) `pose` secondary structure is given by
the `secondary_structure` argument (either a `String` using the 3-mode
categorization - "H" for helix, "E" for beta sheets and "C" for coils; or a
`Function` that takes the given [`Pose`](@ref) `pose` as a single input argument
and returns the 3-mode categorization -
[`ProtoSyn.Peptides.categorize_ss_from_dihedral_angles`](@ref), by default).
Read [`ProtoSyn.Peptides.categorize_ss_from_dihedral_angles`](@ref)
documentation for more details on [`Pose`](@ref) requirements. If an
`AbstractSelection` `selection` is given, only the selected [`Residue`](@ref)
instances are considered for the energy calculation (any given `selection` will
be promoted to [`Residue`](@ref) type, see
[`ProtoSyn.promote`](@ref ProtoSyn.promote)). This method does not calculate
forces. As such, the `update_forces` flag has no effect and serves only for
standardization purposes between `Calculators` methods. An optional parameter
`Type{<: AbstractAccelerationType}` can be provided, stating the acceleration
type used to calculate this energetic contribution (See
[ProtoSyn acceleration types](@ref), if not provided defaults to
`ProtoSyn.acceleration.active`).

!!! ukw "Note:"
    In ProtoSyn version 1.1, the [`ProtoSyn.Peptides.categorize_ss_from_dihedral_angles`](@ref) is notoriously weak. Therefore, manually assigning a secondary structure from a third-party categorization software (such as DSSP) or a ML prediction server (such as RaptorX) is preferable. This may change in future versions of ProtoSyn.

# See also
[`get_default_aa_ss_propensity`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.Calculators.calc_aa_ss_propensity(pose, nothing, false)
(-86.07000000000001, nothing)
```
"""
function calc_aa_ss_propensity(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; secondary_structure::Union{Function, String} = ProtoSyn.Peptides.categorize_ss_from_dihedral_angles, aa_ss_propensity_map::Dict{Char, Dict{Char, T}} = default_aa_ss_propensity) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
    if selection === nothing
        sele = TrueSelection{Residue}()
    else
        sele = ProtoSyn.promote(selection, Residue)
    end

    if isa(secondary_structure, Function)
        ss = secondary_structure(pose)
    else
        ss = secondary_structure
    end

    residues = sele(pose, gather = true)

    e = 0.0
    for residue in residues
        name = ProtoSyn.ProtoSyn.three_2_one[residue.name]
        e -= aa_ss_propensity_map[ss[residue.index]][name]
    end

    return e, nothing
end # function

calc_aa_ss_propensity(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; secondary_structure::Union{Function, String} = ProtoSyn.Peptides.categorize_ss_from_dihedral_angles, aa_ss_propensity_map::Dict{Char, Dict{Char, T}} = default_aa_ss_propensity) where {T <: AbstractFloat} = begin
    calc_aa_ss_propensity(ProtoSyn.acceleration.active, pose, selection, update_forces, secondary_structure = secondary_structure, aa_ss_propensity_map = aa_ss_propensity_map)
end


"""
    fixate_secondary_structure!(efc::EnergyFunctionComponent, pose::Pose)

If the given [`EnergyFunctionComponent`](@ref) `efc` is of type AA_SS_Propensity
and its current `:secondary_structure` setting is a dynamic `Function`, apply
this method to the given [`Pose`](@ref) `pose` to set a static secondary
structure `String`. This increases the [`EnergyFunctionComponent`](@ref) `efc`
performace.

# See also
[`get_default_aa_ss_propensity`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.fixate_secondary_structure!(efc, pose)
"CHHHHHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHHHHHHHHHC"
```
"""
function fixate_secondary_structure!(efc::EnergyFunctionComponent, pose::Pose)
    if :secondary_structure in keys(efc.settings) && isa(efc.settings[:secondary_structure], Function)
        efc.settings[:secondary_structure] = efc.settings[:secondary_structure](pose)
    end
end


"""
    get_default_aa_ss_propensity(;[α::T = 1.0]) where {T <: AbstractFloat}

Return the default aminoacid secondary structure propensity energy
[`EnergyFunctionComponent`](@ref). `α` sets the component weight (on an
[`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance, `1.0`
by default). This function employs
[`calc_aa_ss_propensity`](@ref ProtoSyn.Peptides.Calculators.calc_aa_ss_propensity)
as the `:calc` function.

# Settings
* `aa_ss_propensity_map::Dict{Char, Dict{Char, T}}` - The secondary structure map, correlates a 3-mode category ("H", "E" or "C") to a propensity value, for each aminoacid type (in 1-letter convention).
* `secondary_structure::Union{Function, String}` - The [`Pose`](@ref) secondary structure, either as a static `String` or a dynamic funtion that takes a [`Pose`](@ref) as the single input and returns the secondary structure `String`.

# See also
[`calc_aa_ss_propensity`](@ref)
[`fixate_secondary_structure!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.Calculators.get_default_aa_ss_propensity()
🞧  Energy Function Component:
+---------------------------------------------------+
| Name           | AA_SS_Propensity                 |
| Alpha (α)      | 1.0                              |
| Update forces  | true                             |
| Calculator     | calc_aa_ss_propensity            |
+---------------------------------------------------+
 |    +----------------------------------------------------------------------------------+
 ├──  ● Settings                      | Value                                            |
 |    +----------------------------------------------------------------------------------+
 |    | aa_ss_propensity_map          | Dict{Char, Dict{Char, Float64}}(3 components)    |
 |    | secondary_structure           | categorize_ss_from_dihedral_angles               |
 |    +----------------------------------------------------------------------------------+
 |    
 └──  ○  Selection: nothing
"""
function get_default_aa_ss_propensity(;α::T = 1.0) where {T <: AbstractFloat}
    return EnergyFunctionComponent(
        "AA_SS_Propensity",
        calc_aa_ss_propensity,
        nothing,
        Dict{Symbol, Any}(
            :aa_ss_propensity_map => default_aa_ss_propensity,
            :secondary_structure => ProtoSyn.Peptides.categorize_ss_from_dihedral_angles),
        α,
        true)
end