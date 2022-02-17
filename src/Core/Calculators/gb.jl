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

    const models_onnx = GBModels(
        ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_C.onnx"), rand(Float64, 400, 1)),
        ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_N.onnx"), rand(Float64, 400, 1)),
        ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_H.onnx"), rand(Float64, 400, 1)),
        ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_O.onnx"), rand(Float64, 400, 1)),
        ONNX.load(joinpath(ProtoSyn.resource_dir, "Calculators/iGBR-NN/model_S.onnx"), rand(Float64, 400, 1)))

    Base.show(io::IO, models::GBModels) = begin
        s = ""
        for elem in fieldnames(typeof(models))
            s *= string(elem) * ": "
            s *= getproperty(models, elem) === nothing ? "⨯ | " : "✓ | "
        end
        print(io, s[1:end-2])
    end


    function predict_igbr_nn_born_radii(pose::Pose, selection::Opt{AbstractSelection} = nothing; dm::Opt{Matrix{T}} = nothing, models::GBModels = models_onnx) where {T <: AbstractFloat}
        
        # 1. Calculate distance matrix, if required
        if dm === nothing
            dm = collect(ProtoSyn.Calculators.distance_matrix(pose, selection))
        end
        
        # 2. Pre-calculate histogram for born radii prediction
        hist = ProtoSyn.hist_by_distance_by_elem(pose, selection, dm = dm)
        
        # 3. Predict born-radii
        born_radii = zeros(eltype(pose.state), size(hist)[1])
        for elem in fieldnames(typeof(models))
            model = getproperty(models, elem)
            sele  = FieldSelection{Atom}(string(elem), :symbol) & selection
            radii = Ghost.play!(model, hist')
            mask  = sele(pose).content
            born_radii[mask] .= radii[mask]
        end

        return born_radii
    end

    """
    # TODO DOCUMENTATION

    """
    function calc_gb(::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; born_radii::Union{Vector{T}, Function} = predict_igbr_nn_born_radii, ϵ_protein::T = 1.0, ϵ_solvent::T = 80.0, models::GBModels = models_onnx) where {T <: AbstractFloat}

        # Pre-calculate atomic distances
        atoms  = collect(eachatom(pose.graph))
        natoms = length(atoms)
        dm     = collect(ProtoSyn.Calculators.distance_matrix(pose, selection))

        # Predict born radii if necessary
        if isa(born_radii, Function)
            born_radii = born_radii(pose, selection, dm = dm, models = models)
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

    calc_gb(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; born_radii::Union{Vector{T}, Function} = predict_igbr_nn_born_radii, ϵ_protein::T = 1.0, ϵ_solvent::T = 80.0, models::GBModels = models_onnx) where {T <: AbstractFloat} = begin
        calc_gb(ProtoSyn.acceleration.active, pose, selection, update_forces, born_radii = born_radii, ϵ_protein = ϵ_protein, ϵ_solvent = ϵ_solvent, models = models)
    end

    function get_default_gb(;α::T = 1.0) where {T <: AbstractFloat}
        EnergyFunctionComponent(
            "GB_Solvation",
            calc_gb,
            nothing,
            Dict{Symbol, Any}(
                :born_radii => predict_igbr_nn_born_radii,
                :ϵ_protein  => 1.0,
                :ϵ_solvent  => 80.0,
                :models     => models_onnx
            ),
            α,
            false)
    end

end