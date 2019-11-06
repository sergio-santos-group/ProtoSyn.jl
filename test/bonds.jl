function two_atom_bond(r::Float64, k::Float64)
    type = ProtoSyn.Forcefields.HarmonicBondType("A", r, k)
    bond = ProtoSyn.Forcefields.HarmonicBond(1, 2, type)

    st = ProtoSyn.State(2)
    st.forces = zeros(st.size, 3)
    st.energy = ProtoSyn.Energy()

    [bond], st
end

@testset "2-atom bond: equilibrium" begin
    bonds, state = two_atom_bond(1.0, 1.0)
    state.coords[2,1] = 1.0

    @test ProtoSyn.eval!(state, bonds, false) ≈ 0.0
    @test ProtoSyn.eval!(state, bonds, true) ≈ 0.0
    @test state.forces ≈ [0.0 0.0 0.0; 0.0 0.0 0.0]

end

@testset "2-atom bond: off-equilibrium" begin
    bonds, state = two_atom_bond(1.0, 1.0)
    state.coords[2,1] = 2.0

    @test ProtoSyn.eval!(state, bonds, false) ≈ 0.5
    @test ProtoSyn.eval!(state, bonds, true) ≈ 0.5
    @test state.forces ≈ [1.0 0.0 0.0; -1.0 0.0 0.0]

end