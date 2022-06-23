using ProtoSyn.Units

@testset verbose = true "Units $(repeat("-", 51))" begin

    @testset verbose = true "Angles" begin
        @test 0° == 0 
        @test 1° ≈ deg2rad(1) atol = 1e-10
        @test 180° ≈ pi atol = 1e-10
    end

    @testset verbose = true "Distances" begin
        @test 0Å == 0
        @test 0nm == 0
        @test 0m == 0
        @test 10Å == 1nm
        @test 1e10Å == 1m
    end

    @testset verbose = true "Conversion" begin
        @test tonumber(1) == 1
        @test tonumber("1") == 1
        @test tonumber("1+1") == 2
        @test tonumber("0°") == 0
        @test tonumber("1°") ≈ deg2rad(1) atol = 1e-10
        @test tonumber("1°+pi") ≈ deg2rad(1) + pi atol = 1e-10
        @test tonumber("1+pi") ≈ 1 + pi
    end

    @testset verbose = true "Composite" begin
        @test tonumber("kJ/mol/rad^2") == 1
        @test 1*kJ/mol/rad^2 == 1
    end
    
end