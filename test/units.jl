using ProtoSyn.Units

@testset "angles" begin
    @test 0° == 0 
    @test 1°≈deg2rad(1) atol=1e-10
    @test 180°≈pi atol=1e-10
end

@testset "distances" begin
    @test 0Å == 0
    @test 0nm == 0
    @test 0m == 0
    @test 10Å == 1nm
    @test 1e10Å == 1m
end

@testset "conversion" begin
    @test tonumber(1)==1
    @test tonumber("1")==1
    @test tonumber("1+1")==2
    @test tonumber("0°")==0
    @test tonumber("1°")≈deg2rad(1) atol=1e-10
    @test tonumber("1°+pi")≈deg2rad(1)+pi atol=1e-10
    @test tonumber("1+pi")≈1+pi
end
