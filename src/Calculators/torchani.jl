module TorchANI

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent
    using PyCall

    const torch    = PyNULL()
    const torchani = PyNULL()
    const device   = PyNULL()
    const _model   = PyNULL()

    # --- SINGLE MODEL

    """
        Calculators.calc_torchani_model([::A], pose::Pose; model_index::Int = 3) where {A}
        
    Calculate the pose energy according to a single TorchANI model neural
    network. The model can be defined using `model_index` (from model 1 to 8,
    default is 3). The optional `A` parameter defines the acceleration mode used
    (only CUDA_2 is available). If left undefined the default
    ProtoSyn.acceleration mode will be used.

    #See also:
    `calc_torchani_ensemble`

    # Examples
    ```jldoctest
    julia> Calculators.calc_torchani_model(pose)
    ```
    """
    function calc_torchani_model(::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.SIMD_1}}, pose::Pose; model_index::Int = 3)
        println("ERROR: 'calc_torchani_model' requires CUDA_2 acceleration.")
        return 0.0
    end

    function calc_torchani_model(::Type{ProtoSyn.CUDA_2}, pose::Pose; model_index::Int = 3)
        
        c           = get_cartesian_matrix(pose)
        coordinates = torch.tensor([c], requires_grad = true, device = device).float()
        
        s           = get_ani_species(pose)
        species     = torch.tensor([s], device = device)
        
        m1 = _model.species_converter((species, coordinates))
        m2 = _model.aev_computer(m1)
        return get(_model.neural_networks, model_index)(m2)[2].item()
    end

    calc_torchani_model(pose::Pose; model_index::Int = 3) = begin
        calc_torchani_model(ProtoSyn.acceleration, pose; model_index = model_index)
    end

    # --- ENSEMBLE

    """
        Calculators.calc_torchani_ensemble([::A], pose::Pose) where {A}
    
    Calculate the pose energy according to the whole TorchANI neural
    network ensemble. The optional `A` parameter defines the acceleration mode
    used (only CUDA_2 is available). If left undefined the default
    ProtoSyn.acceleration mode will be used.

    #See also:
    `calc_torchani_model`

    # Examples
    ```jldoctest
    julia> Calculators.calc_torchani_ensemble(pose)
    ```
    """
    function calc_torchani_ensemble(::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.SIMD_1}}, pose::Pose)
        println("ERROR: 'calc_torchani_ensemble' requires CUDA_2 acceleration.")
        return 0.0
    end

    function calc_torchani_ensemble(::Type{ProtoSyn.CUDA_2}, pose::Pose)
        
        c           = get_cartesian_matrix(pose)
        coordinates = torch.tensor([c], requires_grad = true, device = device).float()
        
        s           = get_ani_species(pose)
        species     = torch.tensor([s], device = device)

        m1 = _model.species_converter((species, coordinates))
        m2 = _model.aev_computer(m1)
        return _model.neural_networks(m2)[2].item()
    end

    calc_torchani_ensemble(pose::Pose) = begin
        calc_torchani_ensemble(ProtoSyn.acceleration, pose)
    end


    function __init__()
        copy!(torch, pyimport("torch"))
        copy!(torchani, pyimport("torchani"))

        if torch.cuda.is_available()
            copy!(device, torch.device("cuda"))
        else
            copy!(device, torch.device("cpu"))
        end
        
        copy!(_model, torchani.models.ANI2x(periodic_table_index = true).to(device))
    end

    torchani_model    = EnergyFunctionComponent("TorchANI_ML_Model", calc_torchani_model)
    torchani_ensemble = EnergyFunctionComponent("TorchANI_ML_Ensemble", calc_torchani_ensemble)
end
