using ProtoSyn
using ProtoSyn.Calculators: EnergyFunctionComponent

aa_frequencies_filename = joinpath(
    Peptides.resource_dir, "Calculators/aa-frequencies.yml")

default_aa_frequencies = begin
    _data = ProtoSyn.read_yml(aa_frequencies_filename)
    Dict{Char, Float64}(
        map((x) -> x[1], collect(keys(_data))) .=> collect(values(_data)))
end

function calc_aa_frequency(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; aa_frequency_map::Dict{Char, T} = default_aa_frequencies) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
    if selection === nothing
        sele = TrueSelection{Residue}()
    else
        sele = ProtoSyn.promote(selection, Residue)
    end

    residues = sele(pose, gather = true)

    e = 0.0
    for residue in residues
        name = ProtoSyn.Peptides.three_2_one[residue.name]
        e -= aa_frequency_map[name]
    end

    return e, nothing
end # function

calc_aa_frequency(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; aa_frequency_map::Dict{Char, T} = default_aa_frequencies) where {T <: AbstractFloat} = begin
    calc_aa_frequency(ProtoSyn.acceleration.active, pose, selection, update_forces, aa_frequency_map = aa_frequency_map)
end

function get_default_aa_frequency(;α::T = 1.0) where {T <: AbstractFloat}
    return EnergyFunctionComponent(
        "Natural_AA_Frequency",
        calc_aa_frequency,
        nothing,
        Dict{Symbol, Any}(:aa_frequency_map => default_aa_frequencies),
        α,
        true)
end