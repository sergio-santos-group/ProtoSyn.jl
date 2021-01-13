module TorchANI

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent
    using PyCall

    const torch    = PyNULL()
    const torchani = PyNULL()
    const device   = PyNULL()
    const _model   = PyNULL()

    # --- AUX

    export get_ani_species

    """
    TO DO
    """
    function get_ani_species(container::ProtoSyn.AbstractContainer)
        
        periodic_table = Dict("H" => 1, "C" => 6, "N" => 7, "O" => 8, "S" => 16)

        species = Vector{Int64}()
        for atom in eachatom(container)
            push!(species, periodic_table[atom.symbol])
        end

        return species
    end

    get_ani_species(pose::Pose) = get_ani_species(pose.graph)

    # --- SINGLE MODEL

    """
        Calculators.calc_torchani_model([::A], pose::Pose; model_index::Int = 3) where {A}
        
    Calculate the pose energy according to a single TorchANI model neural
    network. The model can be defined using `model_index` (from model 1 to 8,
    default is 3). The optional `A` parameter defines the acceleration mode used
    (only CUDA_2 is available). If left undefined the default
    ProtoSyn.acceleration.active mode will be used.

    #See also:
    `calc_torchani_ensemble`

    # Examples
    ```jldoctest
    julia> Calculators.calc_torchani_model(pose)
    ```
    """
    function calc_torchani_model(::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.SIMD_1}}, pose::Pose; update_forces::Bool = false, model_index::Int = 3)
        println("ERROR: 'calc_torchani_model' requires CUDA_2 acceleration.")
        return 0.0
    end

    function calc_torchani_model(::Type{ProtoSyn.CUDA_2}, pose::Pose; update_forces::Bool = false, model_index::Int = 3)
        
        c           = get_cartesian_matrix(pose)
        coordinates = torch.tensor([c], requires_grad = true, device = device).float()
        
        s           = get_ani_species(pose)
        species     = torch.tensor([s], device = device)
        
        m1 = _model.species_converter((species, coordinates))
        m2 = _model.aev_computer(m1)
        m3 = get(_model.neural_networks, model_index)(m2)[2]
        if update_forces
            f = torch.autograd.grad(m3.sum(), coordinates)[1][1]
            return m3.item(), convert(Matrix{Float64}, f.cpu().numpy()')
        else
            return m3.item(), nothing
        end
    end

    calc_torchani_model(pose::Pose; update_forces::Bool = false, model_index::Int = 3) = begin
        calc_torchani_model(ProtoSyn.acceleration.active, pose, update_forces = update_forces, model_index = model_index)
    end

    # --- ENSEMBLE

    """
        Calculators.calc_torchani_ensemble([::A], pose::Pose) where {A}
    
    Calculate the pose energy according to the whole TorchANI neural
    network ensemble. The optional `A` parameter defines the acceleration mode
    used (only CUDA_2 is available). If left undefined the default
    ProtoSyn.acceleration.active mode will be used.

    #See also:
    `calc_torchani_model`

    # Examples
    ```jldoctest
    julia> Calculators.calc_torchani_ensemble(pose)
    ```
    """
    function calc_torchani_ensemble(::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.SIMD_1}}, pose::Pose, update_forces::Bool)
        println("ERROR: 'calc_torchani_ensemble' requires CUDA_2 acceleration.")
        return 0.0
    end

    function calc_torchani_ensemble(::Type{ProtoSyn.CUDA_2}, pose::Pose; update_forces::Bool = false)
        
        c           = get_cartesian_matrix(pose)
        coordinates = torch.tensor([c], requires_grad = true, device = device).float()
        
        s           = get_ani_species(pose)
        species     = torch.tensor([s], device = device)

        m1 = _model.species_converter((species, coordinates))
        m2 = _model.aev_computer(m1)
        m3 = _model.neural_networks(m2)[2]
        if update_forces
            return m3.item(), torch.autograd.grad(m3.sum(), coordinates)[1][1].cpu().numpy()
        else
            return m3.item(), nothing
        end
    end

    calc_torchani_ensemble(pose::Pose; update_forces::Bool = false) = begin
        calc_torchani_ensemble(ProtoSyn.acceleration.active, pose, update_forces = update_forces)
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
