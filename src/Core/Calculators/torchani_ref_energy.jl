function calc_torchani_internal_energy(M::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::AbstractSelection, update_forces::Bool = false; static_ref_energy::T = 0.0, use_ensemble::Bool = false, model::Int = 3) where {T <: AbstractFloat}

    sele = ProtoSyn.promote(selection, Residue)

    e = 0.0
    for residue in sele(pose, gather = true)
        if use_ensemble
            _e, _ = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(M, pose, SerialSelection{Residue}(residue.id, :id), false)
        else
            _e, _ = ProtoSyn.Calculators.TorchANI.calc_torchani_model(M, pose, SerialSelection{Residue}(residue.id, :id), false, model = model)
        end
        e += _e
    end

    return e + static_ref_energy
end

function calc_torchani_internal_energy(M::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Nothing, update_forces::Bool = false; static_ref_energy::T = 0.0, use_ensemble::Bool = false, model::Int = 3) where {T <: AbstractFloat}
    calc_torchani_internal_energy(M, pose, TrueSelection{Residue}(), update_forces, static_ref_energy = static_ref_energy, use_ensemble = use_ensemble, model = model)
end

calc_torchani_internal_energy(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; static_ref_energy::T = 0.0, use_ensemble::Bool = false, model::Int = 3) where {T <: AbstractFloat} = begin
    calc_torchani_internal_energy(ProtoSyn.acceleration.active, pose, selection, update_forces, static_ref_energy = static_ref_energy, use_ensemble = use_ensemble, model = model)
end

function fixate_static_ref_energy!(efc::EnergyFunctionComponent, pose::Pose, selection::Opt{AbstractSelection})
    if :static_ref_energy in keys(efc.settings)
        efc.settings[:static_ref_energy] = calc_torchani_internal_energy(pose, selection, use_ensemble = efc.settings[:use_ensemble], model = efc.settings[:model])
    end
end

function get_torchani_internal_energy(;α::T = 1.0) where {T <: AbstractFloat}
    return EnergyFunctionComponent(
        "TorchANI_Ref_Energy",
        calc_torchani_internal_energy,
        nothing,
        Dict{Symbol, Any}(:static_ref_energy => 0.0, :use_ensemble => false, :model => 3),
        α,
        true)
end