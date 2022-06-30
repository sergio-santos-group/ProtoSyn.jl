using StaticArrays

@testset verbose = true "Drivers $(repeat("-", 49))" begin

    @testset verbose = true "$(@sprintf "%-54s" "Thermostats")" begin
        pose = copy(backup)
        t1 = ProtoSyn.Drivers.get_constant_temperature(1.0)
        @test t1(0)   === 1.0
        @test t1(50)  === 1.0
        @test t1(100) === 1.0

        t2 = ProtoSyn.Drivers.get_linear_quench(1.0, 100)
        @test t2(0)   === 1.0
        @test t2(50)  === 0.5
        @test t2(100) === 0.0

        t3 = ProtoSyn.Drivers.get_quadratic_quench(1.0, 100, 0.0)
        @test t3(0)   === 1.0
        @test t3(50)  === 0.25
        @test t3(100) === 0.0
    end

    @testset verbose = true "$(@sprintf "%-54s" "Callbacks")" begin
        pose = copy(backup)

        event = (pose::Pose, driver_state::ProtoSyn.Drivers.DriverState) -> return "OK"
        cb = ProtoSyn.Drivers.Callback(event, 2)

        driver_state = ProtoSyn.Drivers.MonteCarloState{Float64}()
        @test cb(pose, driver_state) === "OK"

        driver_state.step += 1
        @test cb(pose, driver_state) === nothing
    end

    if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2
        @testset verbose = true "$(@sprintf "%-54s" "Monte Carlo")" begin
            pose = copy(backup)
            energy_function = ProtoSyn.Calculators.EnergyFunction([
                ProtoSyn.Calculators.TorchANI.get_default_torchani_ensemble()
            ])
            m = ProtoSyn.Mutators.DihedralMutator(() -> pi, 1.0, 1.0, an"C" & rid"2")
            t1 = ProtoSyn.Drivers.get_constant_temperature(1e10)

            monte_carlo = ProtoSyn.Drivers.MonteCarlo(energy_function, m, nothing, 10, t1)

            @test pose.state.e[:Total] === Inf
            @test pose.state[pose.graph[1][2]["CA"]].Δϕ === 0.0
            pose = monte_carlo(pose)
            @test pose.state.e[:Total] ≈ -0.12801788747310638 atol = 1e-5
            @test pose.state[pose.graph[1][2]["CA"]].Δϕ ≈ 31.41592653589793 atol = 1e-5
            @test pose.state.i2c === false
            @test pose.state.c2i === false
        end
    else
        @warn "Skipping Monte Carlo test on this system: No CUDA found."
    end

    if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2
        @testset verbose = true "$(@sprintf "%-54s" "Steepest Descent")" begin
            pose = copy(backup)
            energy_function = ProtoSyn.Calculators.EnergyFunction([
                ProtoSyn.Calculators.TorchANI.get_default_torchani_ensemble()
            ])
            energy_function.update_forces = true

            sd = ProtoSyn.Drivers.SteepestDescent(energy_function, nothing, 100, 0.001, 0.1)

            T = eltype(pose.state)
            @test pose.state.e[:Total] === Inf
            @test pose.state[pose.graph[1][2]["CA"]].Δϕ === 0.0
            @test pose.state[pose.graph[1][2]["CA"]].ϕ === 3.141592653589793
            @test pose.state[pose.graph[1][2]["N"]].t ≈ [3.762615121185313, -2.6664117066805666, -8.203956756092598e-7]
            pose = sd(pose)
            @test pose.state[pose.graph[1][2]["N"]].t[1] ≈ 3.6981331962060873 atol = 1e-5
            @test pose.state[pose.graph[1][2]["N"]].t[2] ≈ -2.7050704999395263 atol = 1e-5
            @test pose.state[pose.graph[1][2]["N"]].t[3] ≈ -0.04025743850412392 atol = 1e-5
            @test pose.state.e[:Total] ≈ -0.26939642429351807 atol = 1e-5
            @test pose.state[pose.graph[1][2]["CA"]].Δϕ === 0.0
            @test pose.state[pose.graph[1][2]["CA"]].ϕ ≈ -3.115530980993034 atol = 1e-5
            @test pose.state.i2c === false
            @test pose.state.c2i === false
        end
    else
        @warn "Skipping Steepest Descent test on this system: No CUDA found."
    end

    if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2
        @testset verbose = true "$(@sprintf "%-54s" "Iterated Local Search")" begin
            pose = copy(backup)
            energy_function = ProtoSyn.Calculators.EnergyFunction([
                ProtoSyn.Calculators.TorchANI.get_default_torchani_ensemble()
            ])

            m = ProtoSyn.Mutators.DihedralMutator(() -> pi, 1.0, 1.0, an"C" & rid"2")
            t1 = ProtoSyn.Drivers.get_constant_temperature(1e10)
            monte_carlo = ProtoSyn.Drivers.MonteCarlo(energy_function, m, nothing, 10, t1)
            ils = ProtoSyn.Drivers.ILS(energy_function, m, monte_carlo, nothing, 2, t1)

            @test pose.state.e[:Total] === Inf
            @test pose.state[pose.graph[1][2]["CA"]].Δϕ === 0.0
            pose = ils(pose)
            @test pose.state.e[:Total] ≈ 0.3060266077518463 atol = 1e-5
            @test pose.state[pose.graph[1][2]["CA"]].Δϕ ≈ 65.97344572538569 atol = 1e-5
            @test pose.state.i2c === false
            @test pose.state.c2i === false
        end
    else
        @warn "Skipping Iterated Local Search test on this system: No CUDA found."
    end

    if ProtoSyn.acceleration.active === ProtoSyn.CUDA_2
        @testset verbose = true "$(@sprintf "%-54s" "Compound")" begin
            pose = copy(backup)
            energy_function = ProtoSyn.Calculators.EnergyFunction([
                ProtoSyn.Calculators.TorchANI.get_default_torchani_ensemble()
            ])
            energy_function.update_forces = true

            m = ProtoSyn.Mutators.DihedralMutator(() -> pi, 1.0, 1.0, an"C" & rid"2")
            t1 = ProtoSyn.Drivers.get_constant_temperature(1e10)
            monte_carlo = ProtoSyn.Drivers.MonteCarlo(energy_function, m, nothing, 10, t1)
            ils = ProtoSyn.Drivers.ILS(energy_function, m, monte_carlo, nothing, 2, t1)
            cdriver = ProtoSyn.Drivers.CompoundDriver([ils, m, monte_carlo]) 

            @test pose.state.e[:Total] === Inf
            @test pose.state[pose.graph[1][2]["CA"]].Δϕ === 0.0
            pose = cdriver(pose)
            @test pose.state.e[:Total] ≈ -0.12801788747310638 atol = 1e-5
            @test pose.state[pose.graph[1][2]["CA"]].Δϕ ≈ 100.53096491487345 atol = 1e-5
            @test pose.state.i2c === false
            @test pose.state.c2i === false
        end
    else
        @warn "Skipping Compound test on this system: No CUDA found."
    end
end