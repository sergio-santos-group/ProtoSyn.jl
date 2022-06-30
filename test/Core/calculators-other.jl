using PyCall

@testset verbose = true "$(@sprintf "%-54s" "Custom Reference Energy")" begin
    pose = copy(backup)
    ref = ProtoSyn.Calculators.get_default_custom_ref_energy()
    ref.settings[:map] = Dict{AbstractSelection, Float64}(rn"GLU" => 1.0)
    ef = ProtoSyn.Calculators.EnergyFunction([ref])
    e1 = ef(pose)

    @test e1 ≈ 1.0 atol = 1e-5
end 


@testset verbose = true "$(@sprintf "%-54s" "Electrostatics")" begin
    pose = copy(backup)
    coulomb = ProtoSyn.Calculators.Electrostatics.get_default_coulomb()
    ef = ProtoSyn.Calculators.EnergyFunction([coulomb])
    e1 = ef(pose)

    @test e1 ≈ 0.7892871860858698 atol = 1e-5
end 


@testset verbose = true "$(@sprintf "%-54s" "Generalized Born")" begin
    pose = copy(backup)
    gb = ProtoSyn.Calculators.GB.get_default_gb()
    ef = ProtoSyn.Calculators.EnergyFunction([gb])
    e1 = ef(pose)

    @test e1 ≈ 1.12156798870317 atol = 1e-5
end 


@testset verbose = true "$(@sprintf "%-54s" "Hydrogen Bonds")" begin
    pose = copy(backup)
    hb = ProtoSyn.Calculators.HydrogenBonds.get_default_hydrogen_bond_network()
    ef = ProtoSyn.Calculators.EnergyFunction([hb])
    e1 = ef(pose)

    @test e1 ≈ 0.0 atol = 1e-5
end 


@testset verbose = true "$(@sprintf "%-54s" "Radius Gyration")" begin
    pose = copy(backup)
    rg = ProtoSyn.Calculators.RG.get_default_rg()
    ef = ProtoSyn.Calculators.EnergyFunction([rg])
    e1 = ef(pose)

    @test e1 ≈ 6.606306008191303 atol = 1e-5
end 


@testset verbose = true "$(@sprintf "%-54s" "REF-15")" begin
    REF_15_available = false
    pyrosetta = PyNULL()
    try
        copy!(pyrosetta, pyimport("pyrosetta"))
        REF_15_available = true
    catch LoadError
        nothing
    end

    if REF_15_available
        pose = copy(backup)
        ref15 = ProtoSyn.Calculators.REF15.get_default_ref15()
        ef = ProtoSyn.Calculators.EnergyFunction([ref15])
        e1 = ef(pose)

        @test e1 ≈ 13.57779308416654 atol = 1e-5
    end
end


@testset verbose = true "$(@sprintf "%-54s" "SASA")" begin
    pose = copy(backup)
    sasa = ProtoSyn.Calculators.SASA.get_default_sasa()
    ef = ProtoSyn.Calculators.EnergyFunction([sasa])
    e1 = ef(pose)

    @test e1 ≈ 1501.0 atol = 1e-5

    sasa = ProtoSyn.Calculators.SASA.get_default_sasa_energy()
    ef = ProtoSyn.Calculators.EnergyFunction([sasa])
    e1 = ef(pose)

    @test e1 ≈ -901.5000000000002 atol = 1e-5
end 