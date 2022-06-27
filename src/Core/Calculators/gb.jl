module GB

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent
    using ONNX
    using Ghost

    mutable struct GBModels
        C::Opt{Ghost.Tape}
        N::Opt{Ghost.Tape}
        H::Opt{Ghost.Tape}
        O::Opt{Ghost.Tape}
        S::Opt{Ghost.Tape}
    end

    models_onnx = nothing

    function __init__()
        printstyled(" | Loading ONNX models\n", color = :cyan)
        
        @eval(GB, models_onnx = GBModels(
            ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_C.onnx"), rand(Float64, 400, 1)),
            ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_N.onnx"), rand(Float64, 400, 1)),
            ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_H.onnx"), rand(Float64, 400, 1)),
            ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_O.onnx"), rand(Float64, 400, 1)),
            ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_S.onnx"), rand(Float64, 400, 1))))
    end

    Base.show(io::IO, models::GBModels) = begin
        s = ""
        for elem in fieldnames(typeof(models))
            s *= string(elem) * ": "
            s *= getproperty(models, elem) === nothing ? "⨯ | " : "✓ | "
        end
        print(io, s[1:end-2])
    end


    """
        predict_igbr_nn_born_radii(pose::Pose, selection::Opt{AbstractSelection} = nothing; dm::Opt{Matrix{T}} = nothing, models::GBModels = models_onnx) where {T <: AbstractFloat}
    
        Returns the Born Radii for each [`Atom`](@ref) instance in the given
    [`Pose`](@ref) `pose` (selected by the `AbstractSelection` `selection` -
    [`TrueSelection`](@ref), by default), according to the IGBR neural network
    model (See Fogolari et al. work - https://pubmed.ncbi.nlm.nih.gov/31693089/).
    By default, uses the pre-trained `ProtoSyn.Calculators.GB.models_onnx`
    models, but can be set to use other models with the `models` argument.
    Optionally, a pre-calculated full distance matrix `dm` can be provided,
    otherwise will calculate one using the `ProtoSyn.acceleration.active` mode.

    # See also
    [`calc_gb`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.GB.predict_igbr_nn_born_radii(pose)
    1140-element Vector{Float64}:
     0.4011280834674835
     0.3859025537967682
     0.5258529782295227
     (...)
    ```
    """
    function predict_igbr_nn_born_radii(pose::Pose, selection::Opt{AbstractSelection} = nothing; dm::Opt{Matrix{T}} = nothing, models::GBModels = models_onnx) where {T <: AbstractFloat}
        
        if selection === nothing
            sele = TrueSelection{Atom}()
        else
            sele = ProtoSyn.promote(selection, Atom)
        end
        
        # 1. Calculate distance matrix, if required
        if dm === nothing
            dm = collect(ProtoSyn.Calculators.full_distance_matrix(pose, sele))
        end
        
        # 2. Pre-calculate histogram for born radii prediction
        hist = ProtoSyn.hist_by_distance_by_elem(pose, sele, dm = dm)
        
        # 3. Predict born-radii
        n_selected_atoms = size(hist)[1]
        born_radii = zeros(eltype(pose.state), n_selected_atoms)
        atoms = sele(pose, gather = true)
        for elem in fieldnames(typeof(models))
            model = getproperty(models, elem)
            elem  = string(elem)
            sele  = FieldSelection{Atom}(string(elem), :symbol) & selection
            radii = Ghost.play!(model, hist')

            # 3.1) For this specific case, the mask needs to consider only the
            # selected atoms (not all the atoms in the pose). Currently (in
            # ProtoSyn 1.01), only the whole pose or individual graph containers
            # (Residue & Segment) have resolving functions for the selections.
            # In future iterations, selection resolving methods for vectors of
            # Atoms (& Residues, Segments) can be implemented. For now, the 
            # following custom implementation resolves the selection only on the
            # selected atoms.
            mask = Mask{Atom}(n_selected_atoms)
            for (index, atom) in enumerate(atoms)
                if atom.symbol == elem
                    mask[index] = true
                end
            end

            born_radii[mask.content] .= radii[mask.content]
        end

        return born_radii
    end


    """
        calc_gb([::A], pose::Pose, selection::Opt{AbstractSelection}, [update_forces::Bool = false]; [born_radii::Union{Vector{T}, Function} = predict_igbr_nn_born_radii], [ϵ_protein::T = 1.0], [ϵ_solvent::T = 80.0], [models::GBModels = models_onnx]) where {A}
        
    Calculate the [`Pose`](@ref) `pose` solvation energy according to the
    Generalized Born function (make sure the [`Pose`](@ref) `pose` is synched,
    see [`sync!`](@ref)). ProtoSyn uses the Generalzied Born function as
    described in Simmerling et al. work (See
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4361090/). The protein and
    solvent dieletric constants can be defined using the `ϵ_protein` and
    `ϵ_solvent` arguments (default is 1.0 and 80.0, respectively). The Born
    Radii can be statically provided in the `born_radii` argument. Optionally,
    this argument can be a `Function` instance, in which case the Born Radii
    will be re-calculated, using the provided `models`
    (uses `ProtoSyn.Calculators.GB.models_onnx`, by default). Custom Born Radii
    predictors can be defined, using the following signature:

    ```my_born_radii_predictor(pose::Pose, selection::Opt{AbstractSelection} = nothing; dm::Opt{Matrix{T}} = nothing, models::GBModels = models_onnx)```
    
    Note the presence of pre-calculated full-distance matrix `dm` in the
    `Function` arguments. Even if no `GBModels` are used in the Born Radii
    prediction, `calc_gb` still provides them and as such the custom `Function`
    signature should expect them (or use kwargs). An optional parameter `A`
    (`Type{<: AbstractAccelerationType}`) can be provided, stating the
    acceleration type used to calculate this energetic contribution (See
    [ProtoSyn acceleration types](@ref)). Uses `ProtoSyn.acceleration.active` by
    default.
    # See also
    [`get_default_gb`](@ref)
    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.GB.calc_gb(pose, nothing)
    (124.4289232784785, nothing)
    ```
    """
    function calc_gb(::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; born_radii::Union{Vector{T}, Function} = predict_igbr_nn_born_radii, ϵ_protein::T = 1.0, ϵ_solvent::T = 80.0, models::GBModels = models_onnx) where {T <: AbstractFloat}

        if selection !== nothing
            sele = ProtoSyn.promote(selection, Atom)
        else
            sele = TrueSelection{Atom}()
        end

        # Pre-calculate atomic distances
        atoms  = sele(pose, gather = true)
        natoms = length(atoms)
        dm     = collect(ProtoSyn.Calculators.full_distance_matrix(pose, sele))

        # Predict born radii if necessary
        if isa(born_radii, Function)
            born_radii = born_radii(pose, sele, dm = dm, models = models)
        end

        env = - (1/2) * ((1/ϵ_protein) - (1/ϵ_solvent)) # Dieletric term
        int = T(0.0)

        for i in 1:natoms
            atomi = atoms[i]
            qi = pose.state[atomi].δ
            αi = born_radii[i]

            for j in 1:natoms
                i === j && continue

                atomj = atoms[j]
                qj    = pose.state[atomj].δ
                αj    = born_radii[j]

                d_sqr = dm[i, j] * dm[i, j]
                @debug "$i <-> $j"
                @debug " Distance: $(dm[i, j])\n qi: $qi | qj: $qj"
                @debug " αi: $αi | αj: $αj"
                αi = αi < 0 ? 0.0 : αi
                αj = αj < 0 ? 0.0 : αj
                f = sqrt(d_sqr + (αi * αj * exp((-d_sqr) / (4 * αi * αj))))
                int += (qi * qj) / f
            end
        end

        @debug "Env: $env\nInt: $int"
        e = T(env * int)

        return e, nothing
    end

    calc_gb(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; born_radii::Union{Vector{T}, Function} = predict_igbr_nn_born_radii, ϵ_protein::T = 1.0, ϵ_solvent::T = 80.0, models::GBModels = models_onnx) where {T <: AbstractFloat} = begin
        calc_gb(ProtoSyn.acceleration.active, pose, selection, update_forces, born_radii = born_radii, ϵ_protein = ϵ_protein, ϵ_solvent = ϵ_solvent, models = models)
    end


    """
        get_default_gb(;[α::T = 1.0]) where {T <: AbstractFloat}
    Return the default Generalized Born [`EnergyFunctionComponent`](@ref). `α`
    sets the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance). This
    component employs the [`calc_gb`](@ref) method,
    therefore defining a [`Pose`](@ref) energy based on the Generalized Born
    function. By default, this [`EnergyFunctionComponent`](@ref) uses the
    [`predict_igbr_nn_born_radii`](@ref ProtoSyn.Calculators.GB.predict_igbr_nn_born_radii)
    function to predict Born Radii every call. Define
    `efc.settings[:born_radii]` as a `Vector{Float64}` to use static born radii.
    # Settings
    * `born_radii::Union{Function, Vector{T}}` - Defines either the born radii predictor function or static list of born radii (where T <: AbstractFloat);
    * `ϵ_protein::T` - Define the protein dieletric constant (where T <: AbstractFloat);
    * `ϵ_solvent::T` - Define the solvent dieletric constant (where T <: AbstractFloat);
    * `models::GBModels` - Define the `GBModels` to use if `:born_radii` is a predictor function.
    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.GB.get_default_gb()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | GB_Solvation                     |
    | Alpha (α)      | 1.0                              |
    | Update forces  | false                            |
    | Calculator     | calc_gb                          |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | ϵ_solvent                     | 80.0                                             |
     |    | models                        | C: ✓ | N: ✓ | H: ✓ | O: ✓ | S: ✓                 |
     |    | ϵ_protein                     | 4.0                                              |
     |    | born_radii                    | predict_igbr_nn_born_radii                       |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ○  Selection: nothing
    ```
    """
    function get_default_gb(;α::T = 1.0) where {T <: AbstractFloat}
        EnergyFunctionComponent(
            "GB_Solvation",
            calc_gb,
            nothing,
            Dict{Symbol, Any}(
                :born_radii => predict_igbr_nn_born_radii,
                :ϵ_protein  => 4.0,
                :ϵ_solvent  => 80.0,
                :models     => models_onnx
            ),
            α,
            false)
    end

end