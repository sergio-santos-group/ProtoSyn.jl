using ProtoSyn
using ProtoSyn.Calculators: EnergyFunctionComponent

aa_ss_propensity_filename = joinpath(
    Peptides.resource_dir, "Calculators/aa-ss-propensity.yml")

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
# TODO Documentation
Make sure pose is correctly indexed
"""
function calc_aa_ss_propensity(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; secondary_structure::Union{Function, String} = ProtoSyn.Peptides.catergorize_ss_from_dihedral_angles, aa_ss_propensity_map::Dict{Char, Dict{Char, T}} = default_aa_ss_propensity) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
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
        name = ProtoSyn.Peptides.three_2_one[residue.name]
        e -= aa_ss_propensity_map[ss[residue.index]][name]
    end

    return e, nothing
end # function

calc_aa_ss_propensity(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; secondary_structure::Union{Function, String} = ProtoSyn.Peptides.catergorize_ss_from_dihedral_angles, aa_ss_propensity_map::Dict{Char, Dict{Char, T}} = default_aa_ss_propensity) where {T <: AbstractFloat} = begin
    calc_aa_ss_propensity(ProtoSyn.acceleration.active, pose, selection, update_forces, secondary_structure = secondary_structure, aa_ss_propensity_map = aa_ss_propensity_map)
end

"""
# TODO Documentation
"""
function fixate_secondary_structure!(efc::EnergyFunctionComponent, pose::Pose)
    if :secondary_structure in keys(efc.settings) && isa(efc.settings[:secondary_structure], Function)
        efc.settings[:secondary_structure] = efc.settings[:secondary_structure](pose)
    end
end

function get_default_aa_ss_propensity(;α::T = 1.0) where {T <: AbstractFloat}
    return EnergyFunctionComponent(
        "AA_Secondary_Structure_Propensity",
        calc_aa_ss_propensity,
        nothing,
        Dict{Symbol, Any}(
            :aa_ss_propensity_map => default_aa_ss_propensity,
            :secondary_structure => ProtoSyn.Peptides.catergorize_ss_from_dihedral_angles),
        α,
        true)
end