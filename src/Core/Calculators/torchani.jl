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
        printstyled(" | Loading TorchANI\n", color = :cyan)

        torch_is_available = false
        try
            copy!(torch, pyimport("torch"))
            torch_is_available = true
        catch LoadError
            if !("JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC" in keys(ENV))
                ENV["JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC"] = true
            end
            if ENV["JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC"] === "true"
                println()
                @warn """
                📍 ProtoSyn was not able to identify `torch` in this system.
                PyCall is currently configured to use the Python version at $(PyCall.current_python()).
                In order to use the TorchANI energy function component, make sure to:
                    - Set ENV["PYTHON"] to the path of python executable you wish to use, run Pkg.build("PyCall") and re-launch Julia and ProtoSyn.
                    - Make sure `torch` is installed in the machine trying to load ProtoSyn.

                In order to install `torch`, follow the following instructions:
                (1) Install Torch: pip install --pre torch torchvision -f https://download.pytorch.org/whl/nightly/cu100/torch_nightly.html
                (2) Re-launch Julia and ProtoSyn

                ProtoSyn will continue loading, but the `Calculators.TorchANI` module will be unavailable.
                To surpress further warnings for unavailable energy function components, set the JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC environment flag and re-launch Julia and ProtoSyn. 
                \$ export JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC=false
                Optionally, add the above line to ~/.bashrc to persistently supress warnings in further sessions.

                """
            end
        end

        torchani_is_available = false
        try
            copy!(torchani, pyimport("torchani"))
            torchani_is_available = true
        catch LoadError
            if !("JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC" in keys(ENV))
                ENV["JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC"] = true
            end
            if ENV["JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC"] === "true"
                println()
                @warn """
                📍 ProtoSyn was not able to identify `torchani` in this system.
                PyCall is currently configured to use the Python version at $(PyCall.current_python()).
                In order to use the TorchANI energy function component, make sure to:
                    - Set ENV["PYTHON"] to the path of python executable you wish to use, run Pkg.build("PyCall") and re-launch Julia and ProtoSyn.
                    - Make sure `torchani` is installed in the machine trying to load ProtoSyn.

                In order to install `torchani`, follow the following instructions:
                (1) Install TorchANI: pip install torchani
                (2) Re-launch Julia and ProtoSyn

                ProtoSyn will continue loading, but the `Calculators.TorchANI` module will be unavailable.
                To surpress further warnings for unavailable energy function components, set the JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC environment flag and re-launch Julia and ProtoSyn. 
                \$ export JULIA_PROTOSYN_WARN_NON_AVALIABLE_EFC=false
                Optionally, add the above line to ~/.bashrc to persistently supress warnings in further sessions.

                """
            end
        end

        if torch_is_available
            if torch.cuda.is_available()
                copy!(device, torch.device("cuda"))
            else
                copy!(device, torch.device("cpu"))
            end
        end

        if torchani_is_available
            copy!(_model, torchani.models.ANI2x(periodic_table_index = true).to(device))
        end

        if torch_is_available && torchani_is_available
            # @info "TorchANI is using:"
            # @info " torch version $(torch.__version__)"
            # @info " cuda-toolkit version $(torch.version.cuda)"
            # @info " torchani version $(torchani.__version__)"
        end
    end

    # --- AUX

    export get_ani_species

    """
        get_ani_species(container::ProtoSyn.AbstractContainer)

    Return a `Vector{Int}` with the atomic number of each [`Atom`](@ref)
    instance in the given `AbstractContainer` `container`, according to a
    periodic table.

    # See also
    [`calc_torchani_model`](@ref) [`calc_torchani_ensemble`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.TorchANI.get_ani_species(pose.graph[1][1])
    11-element Vector{Int64}:
     7
     1
     6
     1
     6
     1
     1
     8
     1
     6
     8
    ```
    """
    function get_ani_species(pose::Pose, selection::Opt{AbstractSelection})
        
        periodic_table = Dict("H" => 1, "C" => 6, "N" => 7, "O" => 8, "S" => 16, "Se" => 34)

        if selection === nothing
            sele = TrueSelection{Atom}()
        else
            sele = ProtoSyn.promote(selection, Atom)
        end

        species = Vector{Int64}()
        for atom in sele(pose, gather = true)
            push!(species, periodic_table[atom.symbol])
        end

        return species
    end

    # get_ani_species(pose::Pose, selection::Opt{AbstractSelection}) = get_ani_species(pose, selection)

    # --- SINGLE MODEL
    
    """
        calc_torchani_model([::A], pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; model::Int = 3) where {A}
        
    Calculate and return the [`Pose`](@ref) `pose` energy according to a single
    TorchANI model neural network. The model can be defined using `model_index`
    (from model 1 to 8, default is 3).The optional `A` parameter defines the
    acceleration type used. If left undefined the default
    `ProtoSyn.acceleration.active` mode will be used. By setting the
    `update_forces` flag to `true` (`false` by default), this function will also
    calculate and return the forces acting on each atom based on a single
    TorchANI model neural network.

    !!! ukw "Note:"
        Each `model` will return a slightly different value for the energy of the molecular system. Use [`calc_torchani_ensemble`](@ref) for a more accurate (and slow) energy prediction. However, as in most cases the energy value is used in comparison with multiple states/frames, [`calc_torchani_ensemble`](@ref) prediction may be suficient.

    # See also:
    [`calc_torchani_ensemble`](@ref) [`calc_torchani_model_xmlrpc`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose)
    (-0.12573561072349548, nothing)

    julia> ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose, true)
    (-0.12573561072349548, [ ... ])
    ```

    !!! ukw "Note:"
        In ProtoSyn >= 1.0, this function has a memory leak on the Python call. Multiple calls to [`calc_torchani_model`](@ref) require often `GC.gc(false)` calls to impede the 'CUDA out of memory' error. In order to prevent/automate this process, consider the following options:

        (1) - Use an [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) struct (with automatic calls to `GC.gc(false)`, as in:

        `EnergyFunction([ProtoSyn.Calculators.TorchANI.get_default_torchani_model()])`

        (2) - Use [`calc_torchani_model_xmlrpc`](@ref) instead.
    """
    function calc_torchani_model(::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; model::Int = 3)
        
        if selection === nothing
            sele = TrueSelection{Atom}()
        else
            sele = ProtoSyn.promote(selection, Atom)
        end

        coords = pose.state.x.coords[:, sele(pose).content]'
        coordinates = torch.tensor([coords], requires_grad = update_forces, device = device).float()
        s           = get_ani_species(pose, selection)
        species     = torch.tensor([s], device = device)
        m1 = _model.species_converter((species, coordinates))
        m2 = _model.aev_computer(m1)
        m3 = get(_model.neural_networks, model)(m2)[2]
        if update_forces
            f = torch.autograd.grad(m3.sum(), coordinates)[1][1]
            _f = convert(Matrix{Float64}, f.cpu().numpy()' .* -1)
            e = m3.item()

            return e, _f
        else
            e = m3.item()
            # Note: When not consuming the m3 gradient (if requires_grad is set
            # to true), as is the case when update_forces is set to false (this
            # function returns e, nothing), a memory leak on the python site
            # occurs. Two workarounds have been identified:
            # (1) - Set requires_grad = update_forces.
            # (2) - Call GC.gc(false)
            # Fix 1 causes slightly higher memory allocation that option 2.
            # Instead of calling GC.gc(false) every step, a counter can be used
            # to call it every N steps. This can be done either before returning
            # from this function or outside. For ProtoSyn >= 1.0, GC.gc(false)
            # is being called outside, in the encompassing EnergyFunction call.
            # This seems to result in a slightly lower mameory allocation. 
            # However, if for some reason this function is being called solo, 
            # CUDA out of memory erros might occur is users are not aware of
            # this issue. Future iterations of ProtoSyn may bring GC.gc(false)
            # inside the scope of this function.
            return e, nothing
        end
    end

    calc_torchani_model(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; model::Int = 3) = begin
        calc_torchani_model(ProtoSyn.acceleration.active, pose, selection, update_forces, model = model)
    end

    # --- ENSEMBLE

    """
        Calculators.calc_torchani_ensemble([::A], pose::Pose, update_forces::Bool = false) where {A <: ProtoSyn.AbstractAccelerationType}
    
    Calculate and return the [`Pose`](@ref) `pose` energy according to the whole
    TorchANI neural network ensemble. The optional `A` parameter defines the
    acceleration type used. If left undefined the default
    `ProtoSyn.acceleration.active` mode will be used. By setting the
    `update_forces` flag to `true` (`false` by default), this function will also
    calculate and return the forces acting on each atom based on the whole
    TorchANI neural network ensemble.

    # See also:
    [`calc_torchani_model`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(pose)
    (-0.12801790237426758, nothing)

    julia> ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(pose, true)
    (-0.12801788747310638, [ ... ])
    ```
    """
    function calc_torchani_ensemble(::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false)
        
        if selection === nothing
            sele = TrueSelection{Atom}()
        else
            sele = ProtoSyn.promote(selection, Atom)
        end

        coords = pose.state.x.coords[:, sele(pose).content]'
        coordinates = torch.tensor([coords], requires_grad = true, device = device).float()
        
        s           = get_ani_species(pose, selection)
        species     = torch.tensor([s], device = device)

        m1 = _model.species_converter((species, coordinates))
        m2 = _model.aev_computer(m1)
        # m3 = _model.neural_networks(m2)[2]
        m3 = _model.neural_networks((m2[1].float(), m2[2].float()))[2]
        if update_forces
            f = torch.autograd.grad(m3.sum(), coordinates)[1][1]
            return m3.item(), convert(Matrix{Float64}, f.cpu().numpy()').*-1
        else
            return m3.item(), nothing
        end
    end

    calc_torchani_ensemble(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false) = begin
        calc_torchani_ensemble(ProtoSyn.acceleration.active, pose, selection, update_forces)
    end

    # * Default Energy Components ----------------------------------------------

    """
        get_default_torchani_model(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default TorchANI model [`EnergyFunctionComponent`](@ref). `α`
    sets the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance). This
    component employs the [`calc_torchani_model`](@ref) method, therefore
    predicting a structure's TorchANI energy based on a single model.

    # Settings
    * `model::Int` - Defines which model of the TorchANI ensemble to use.

    # See also
    [`calc_torchani_model`](@ref) [`get_default_torchani_ensemble`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.TorchANI.get_default_torchani_model()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | TorchANI_ML_Model                |
    | Alpha (α)      | 1.0                              |
    | Update forces  | true                             |
    | Calculator     | calc_torchani_model              |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | model                         | 3                                                |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ○  Selection: nothing
    ```
    """
    function get_default_torchani_model(;α::T = 1.0) where {T <: AbstractFloat}
        EnergyFunctionComponent(
            "TorchANI_ML_Model",
            calc_torchani_model,
            nothing,
            Dict{Symbol, Any}(:model => 3),
            α,
            true)
    end
    
    """
        get_default_torchani_ensemble(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default TorchANI ensemble [`EnergyFunctionComponent`](@ref). `α`
    sets the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance). This
    component employs the [`calc_torchani_ensemble`](@ref) method, therefore
    predicting a structure's TorchANI energy based on the whole TorchANI
    ensemble (Note: This can be very slow).

    # See also
    [`calc_torchani_ensemble`](@ref) [`get_default_torchani_model`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.TorchANI.get_default_torchani_ensemble()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | TorchANI_ML_Ensemble             |
    | Alpha (α)      | 1.0                              |
    | Update forces  | true                             |
    | Calculator     | calc_torchani_ensemble           |
    +---------------------------------------------------+
     └──  ○  Selection: nothing
    ```
    """
    function get_default_torchani_ensemble(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "TorchANI_ML_Ensemble",
            calc_torchani_ensemble,
            nothing,
            Dict{Symbol, Any}(),
            α,
            true)
    end

    include("torchani_xmlrpc.jl")
    include("torchani_ref_energy.jl")
end