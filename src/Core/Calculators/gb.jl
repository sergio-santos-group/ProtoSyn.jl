module GB

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent
    # using PyCall

    # const tf     = PyNULL()
    # const models_keras = Dict{String, PyObject}(
    #     "C" => PyNULL(),
    #     "N" => PyNULL(),
    #     "H" => PyNULL(),
    #     "O" => PyNULL(),
    #     "S" => PyNULL(),
    # )

    # function __init__()
    #     copy!(tf, pyimport("tensorflow"))
    #     filedir = joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN")
        
    #     # Setting debug level on tensorflow
    #     if ProtoSyn.verbose.mode === false
    #         ENV["KMP_WARNINGS"] = "0"
    #         tf.get_logger().setLevel("ERROR")
    #         ENV["TF_CPP_MIN_LOG_LEVEL"] = "3"
    #     end
        
    #     # Load model from the corresponding json and h5 file
    #     function model_from_json!(t::PyObject, json_m::String, weight_m::String)
    #         json_file = open(joinpath(filedir, json_m), "r")
    #         loaded_model_json = read(json_file, String)
    #         close(json_file)
    #         copy!(t, tf.keras.models.model_from_json(loaded_model_json))
    #         t.load_weights(weight_m)
    #     end

    #     # Populate the models
    #     for key in keys(models_keras)
    #         model_from_json!(models_keras[key], "model_$key.json", "model_$key.h5")
    #     end
    # end

    using ONNX
    using Ghost
    const models_onnx = Dict{String, Ghost.Tape}(
        "C" => ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_C.onnx"), rand(Float64, 400, 1)),
        "N" => ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_N.onnx"), rand(Float64, 400, 1)),
        "H" => ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_H.onnx"), rand(Float64, 400, 1)),
        "O" => ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_O.onnx"), rand(Float64, 400, 1)),
        "S" => ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_S.onnx"), rand(Float64, 400, 1)),
    )


    function predict_igbr_nn_born_radii(pose::Pose, selection::Opt{AbstractSelection} = nothing; dm::Opt{Matrix{T}} = nothing, models::Dict{String, <: Any} = models_keras, batch_mode::Bool = true) where {T <: AbstractFloat}
        
        # 1. Calculate distance matrix, if required
        if dm === nothing
            dm = collect(ProtoSyn.Calculators.distance_matrix(pose, selection))
        end
        
        # 2. Pre-calculate histogram for born radii prediction
        hist = ProtoSyn.hist_by_distance_by_elem(pose, selection, dm = dm)
                        
        
        # 3. Predict born-radii
        if batch_mode
            # Batch processing (≈ 2x slower)
            born_radii = zeros(eltype(pose.state), size(hist)[1])
            for (elem, model) in models
                sele  = FieldSelection{Atom}(elem, :symbol) & selection
                # radii = model.predict(hist)
                radii = Ghost.play!(model, hist')
                mask  = sele(pose).content
                born_radii[mask] .= radii[mask]
            end
        else
            # Element-wise processing
            born_radii = Vector{eltype(pose.state)}()
            for (i, atom) in enumerate(eachatom(pose.graph))
                lane = hist[i, :]
                # push!(born_radii, models[atom.symbol].predict(lane)[1])
                push!(born_radii, Ghost.play!(models[atom.symbol], lane)[1])
            end
        end

        return born_radii
    end

    """
    # TODO DOCUMENTATION

    """
    function calc_gb(::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; born_radii::Opt{Vector{T}} = nothing, ϵ_protein::T = 1.0, ϵ_solvent::T = 80.0, models::Dict{String, <: Any} = models_keras, batch_mode::Bool = true) where {T <: AbstractFloat}

        # Pre-calculate atomic distances
        atoms  = collect(eachatom(pose.graph))
        natoms = length(atoms)
        dm     = collect(ProtoSyn.Calculators.distance_matrix(pose, selection))

        if born_radii === nothing
            born_radii = predict_igbr_nn_born_radii(pose, selection, dm = dm, models = models, batch_mode = batch_mode)
        end

        env = (1/(8*π*ϵ_protein)) * (1 - (1/ϵ_solvent)) # Dieletric term
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
                if ProtoSyn.verbose.mode
                    println("$i <-> $j")
                    println(" Distance: $(dm[i, j])\n qi: $qi | qj: $qj")
                    println(" αi: $αi | αj: $αj")
                end
                f = sqrt(d_sqr + (αi * αj * exp((-d_sqr) / (4 * αi * αj))))
                int += (qi * qj) / f
            end
        end

        ProtoSyn.verbose.mode && println("Env: $env\nInt: $int")
        e = T(env * int)

        return e, nothing
    end

    calc_gb(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; born_radii::Opt{Vector{T}} = nothing, ϵ_protein::T = 1.0, ϵ_solvent::T = 80.0, models::Dict{String, <: Any} = models_keras, batch_mode::Bool = true) where {T <: AbstractFloat} = begin
        calc_gb(ProtoSyn.acceleration.active, pose, selection, update_forces, born_radii = born_radii, ϵ_protein = ϵ_protein, ϵ_solvent = ϵ_solvent, models = models, batch_mode = batch_mode)
    end

    function get_default_gb(;α::T = 1.0) where {T <: AbstractFloat}
        EnergyFunctionComponent(
            "GB_Solvation",
            calc_gb,
            nothing,
            Dict{Symbol, Any}(
                born_radii => nothing,
                ϵ_protein  => 1.0,
                ϵ_solvent  => 80.0,
                models     => models_onnx,
                batch_mode => true),
            α,
            false)
    end

end