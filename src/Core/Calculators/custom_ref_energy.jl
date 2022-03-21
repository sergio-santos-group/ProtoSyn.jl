"""
# TODO: Documentation
Careful with selection. Promotes to bigger selection type.
"""
function calc_custom_ref_energy(M::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::AbstractSelection, update_forces::Bool = false; map::Dict{AbstractSelection, T} = 0.0) where {T <: AbstractFloat}

    selection_type = ProtoSyn.selection_type(selection)

    e = 0.0
    for (map_sele, weight) in map
        map_sele_type = ProtoSyn.selection_type(map_sele)
        if map_sele_type > selection_type
            sele = ProtoSyn.promote(selection, map_sele_type) & map_sele
        elseif selection_type > map_sele_type
            sele = ProtoSyn.promote(map_sele, selection_type) & selection
        end

        e += count(sele(pose).content) * weight
    end

    return e, nothing
end

function calc_torchani_internal_energy(M::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Nothing, update_forces::Bool = false; static_ref_energy::T = 0.0, use_ensemble::Bool = false, model::Int = 3) where {T <: AbstractFloat}
    calc_torchani_internal_energy(M, pose, TrueSelection{Residue}(), update_forces, static_ref_energy = static_ref_energy, use_ensemble = use_ensemble, model = model)
end

calc_torchani_internal_energy(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; static_ref_energy::T = 0.0, use_ensemble::Bool = false, model::Int = 3) where {T <: AbstractFloat} = begin
    calc_torchani_internal_energy(ProtoSyn.acceleration.active, pose, selection, update_forces, static_ref_energy = static_ref_energy, use_ensemble = use_ensemble, model = model)
end

function get_default_custom_ref_energy(;α::T = 1.0) where {T <: AbstractFloat}
    return EnergyFunctionComponent(
        "Custom_Ref_Energy",
        calc_custom_ref_energy,
        nothing,
        Dict{Symbol, Any}(:map => Dict{AbstractSelection, T}()),
        α,
        true)
end