function three_atom_angle(r::Float64, k::Float64)
    type = ProtoSyn.Forcefields.HarmonicAngleType("A", r, k)
    angle = ProtoSyn.Forcefields.HarmonicAngle(1, 2, 3, type)
    
    st = ProtoSyn.State(3)
    st.forces = zeros(st.size, 3)
    st.energy = ProtoSyn.Energy()

    [angle],st
end

@testset "harmonic angle: equilibrium" begin
    # 1----2----3
    angles, state = three_atom_angle(1π, 1.0)
    state.coords[1,1] = -1.0
    state.coords[3,1] =  1.0

    @test ProtoSyn.eval!(state, angles, false) ≈ 0.0
    @test ProtoSyn.eval!(state, angles, true) ≈ 0.0
    @test state.forces ≈ [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]

end

@testset "harmonic angle: off-equilibrium" begin
    #      3
    #      |
    # 1----2
    angles, state = three_atom_angle(0.5π, 1.0)
    
    θ = 0.5π
    θ = 0.95π
    # θ = 1π
    state.coords[1,1] = -1.0
    state.coords[3,1] = cos(π-θ)
    state.coords[3,2] = sin(π-θ)
    
    @test ProtoSyn.eval!(state, angles, false) ≈ 0.5*(θ-angles[1].type.θ)^2
    @test ProtoSyn.eval!(state, angles, true) ≈ 0.5*(θ-angles[1].type.θ)^2
    # @test state.forces ≈ [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]

end

