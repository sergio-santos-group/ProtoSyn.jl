module TorchANI

    # NOTE: `server` is a reserved variable in this module.

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent
    using PyCall

    const torch    = PyNULL()
    const torchani = PyNULL()
    const device   = PyNULL()
    const _model   = PyNULL()

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

    # --- AUX

    export get_ani_species

    """
        # TODO
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
        Calculators.calc_torchani_model([::A], pose::Pose; update_forces::Bool = false, model::Int = 3) where {A}
        
    Calculate the pose energy according to a single TorchANI model neural
    network. The model can be defined using `model_index` (from model 1 to 8,
    default is 3). The optional `A` parameter defines the acceleration mode used
    (only CUDA_2 is available). If left undefined the default
    ProtoSyn.acceleration.active mode will be used. If `update_forces` is set to
    true (false, by default), return the calculated forces on each atom as well.

    # See also:
    `calc_torchani_ensemble` `calc_torchani_model_xmlrpc`

    # Examples
    ```jldoctest
    julia> Calculators.calc_torchani_model(pose)
    (...)
    ```

    # Notes:
    In ProtoSyn <=0.4, this function has a memory leak on the python side.
    Multiple calls to `calc_torchani_model` require often `GC.gc(false)` calls
    to impede the 'CUDA out of memory' error. In order to prevent/automate this
    process, consider the following options:

    (1) - Use an EnergyFunction struct (with automatic calls to `GC.gc(false)`,
    as in:

    energy_function = ProtoSyn.Calculators.EnergyFunction(
        [ProtoSyn.Calculators.TorchANI.get_default_torchani_model(α = 1.0)])

    (2) - Call `calc_torchani_model_xmlrpc` instead.

    """
    function calc_torchani_model(::Union{Type{ProtoSyn.SISD_0}, Type{ProtoSyn.SIMD_1}}, pose::Pose, update_forces::Bool = false; model::Int = 3)
        println("ERROR: 'calc_torchani_model' requires CUDA_2 acceleration.")
        return 0.0, nothing
    end

    function calc_torchani_model(::Type{ProtoSyn.CUDA_2}, pose::Pose, update_forces::Bool = false; model::Int = 3)
        
        coordinates = torch.tensor([pose.state.x.coords'], requires_grad = true, device = device).float()
        
        s           = get_ani_species(pose)
        species     = torch.tensor([s], device = device)
        
        m1 = _model.species_converter((species, coordinates))
        m2 = _model.aev_computer(m1)
        m3 = get(_model.neural_networks, model)(m2)[2]
        if update_forces
            f = torch.autograd.grad(m3.sum(), coordinates)[1][1]
            return m3.item(), convert(Matrix{Float64}, f.cpu().numpy()').*-1
        else
            return m3.item(), nothing
        end
    end

    calc_torchani_model(pose::Pose, update_forces::Bool = false; model::Int = 3) = begin
        calc_torchani_model(ProtoSyn.acceleration.active, pose, update_forces, model = model)
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

    function calc_torchani_ensemble(::Type{ProtoSyn.CUDA_2}, pose::Pose, update_forces::Bool = false)
        
        coordinates = torch.tensor([pose.state.x.coords'], requires_grad = true, device = device).float()
        
        s           = get_ani_species(pose)
        species     = torch.tensor([s], device = device)

        m1 = _model.species_converter((species, coordinates))
        m2 = _model.aev_computer(m1)
        m3 = _model.neural_networks(m2)[2]
        if update_forces
            f = torch.autograd.grad(m3.sum(), coordinates)[1][1]
            return m3.item(), convert(Matrix{Float64}, f.cpu().numpy()').*-1
        else
            return m3.item(), nothing
        end
    end

    calc_torchani_ensemble(pose::Pose, update_forces::Bool = false) = begin
        calc_torchani_ensemble(ProtoSyn.acceleration.active, pose, update_forces)
    end

    # * Default Energy Components ----------------------------------------------

    """
        get_default_torchani_model(;α::T = 1.0) where {T <: AbstractFloat}

    Return the default TorchANI model `EnergyFunctionComponent`. `α`
    sets the component weight (on an `EnergyFunction`). This component employs
    `calc_torchani_model`, therefore predicting a structure's TorchANI energy
    based on a single model.

    # TorchANI model energy settings
    - :model -> defines which model of the TorchANI ensemble to use.

    # See also
    `calc_torchani_model` | `get_default_torchani_ensemble`

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.TorchANI.get_default_torchani_model()
             Name : TorchANI_ML_Model
       Weight (α) : 1.0
    Update forces : true
          Setings :
        :model => 3
    ```
    """
    function get_default_torchani_model(;α::T = 1.0) where {T <: AbstractFloat}
        EnergyFunctionComponent("TorchANI_ML_Model", calc_torchani_model, Dict{Symbol, Any}(:model => 3), α, true)
    end
    
    """
        get_default_torchani_model(;α::T = 1.0) where {T <: AbstractFloat}

    Return the default TorchANI ensemble `EnergyFunctionComponent`. `α`
    sets the component weight (on an `EnergyFunction`). This component employs
    `calc_torchani_ensemble`, therefore predicting a structure's TorchANI energy
    based on the whole TorchANI ensemble (*Note:* This can be very slow).

    # See also
    `calc_torchani_ensemble` `get_default_torchani_model`

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.TorchANI.get_default_torchani_ensemble()
             Name : TorchANI_ML_Ensemble
       Weight (α) : 1.0
    Update forces : true
          Setings : -
    ```
    """
    function get_default_torchani_ensemble(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent("TorchANI_ML_Ensemble", calc_torchani_ensemble, Dict{Symbol, Any}(), α, true)
    end

    include("torchani_xmlrpc.jl")
end
