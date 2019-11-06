@testset "2-atom LJ: r=rmin" begin
    rmin = 2.0
    ϵ = 1.0
    q = 0.0
    m = 1.0
    z = 1
    natoms = 2

    atomtype = ProtoSyn.Forcefields.AtomType("X", z, m, rmin / (2.0^(1/6)), ϵ)
    atoms = [ProtoSyn.Forcefields.Atom(q, atomtype) for i=1:natoms]
    
    state = ProtoSyn.State(natoms)
    state.forces = zeros(natoms, 3)
    state.energy = ProtoSyn.Energy()

    state.coords[2,1] = rmin

    @test ProtoSyn.eval!(state, atoms, false) ≈ -ϵ
    @test ProtoSyn.eval!(state, atoms, true) ≈ -ϵ
    @test state.forces ≈ [0.0 0.0 0.0; 0.0 0.0 0.0] atol=1e-14

end


@testset "2-atom LJ: r=sigma" begin
    rmin = 2.0
    ϵ = 1.0
    q = 0.0
    m = 1.0
    z = 1
    natoms = 2

    atomtype = ProtoSyn.Forcefields.AtomType("X", z, m, rmin / (2.0^(1/6)), ϵ)
    atoms = [ProtoSyn.Forcefields.Atom(q, atomtype) for i=1:natoms]
    
    state = ProtoSyn.State(natoms)
    state.forces = zeros(natoms, 3)
    state.energy = ProtoSyn.Energy()

    state.coords[2,1] = atomtype.σ

    @test ProtoSyn.eval!(state, atoms, false) ≈ 0.0
    @test ProtoSyn.eval!(state, atoms, true) ≈ 0.0
    @test state.forces[1,1] < 0.0 && state.forces[2,1] > 0.0
end



@testset "2-atom LJ: r=10sigma" begin
    rmin = 2.0
    ϵ = 1.0
    q = 0.0
    m = 1.0
    z = 1
    natoms = 2

    atomtype = ProtoSyn.Forcefields.AtomType("X", z, m, rmin / (2.0^(1/6)), ϵ)
    atoms = [ProtoSyn.Forcefields.Atom(q, atomtype) for i=1:natoms]
    
    state = ProtoSyn.State(natoms)
    state.forces = zeros(natoms, 3)
    state.energy = ProtoSyn.Energy()

    state.coords[2,1] = 10*atomtype.σ

    @test ProtoSyn.eval!(state, atoms, false) > -ϵ/60
    @test ProtoSyn.eval!(state, atoms, true) > -ϵ/60
    @test state.forces[1,1] > 0.0 && state.forces[2,1] < 0.0
    
end





@testset "2-atom Coulomb: r=1, q1=q2=-1" begin
    natoms = 2

    atomtype = ProtoSyn.Forcefields.AtomType("X", -1, -1, 1.0, 0.0)
    atoms = [
        ProtoSyn.Forcefields.Atom(-1.0, atomtype),
        ProtoSyn.Forcefields.Atom(-1.0, atomtype)
    ]
    
    state = ProtoSyn.State(natoms)
    state.forces = zeros(natoms, 3)
    state.energy = ProtoSyn.Energy()

    state.coords[2,1] = 1.0

    @test ProtoSyn.eval!(state, atoms, false) ≈ 1.0
    @test ProtoSyn.eval!(state, atoms, true) ≈ 1.0
    @test state.forces[1,1] < 0.0 && state.forces[2,1] > 0.0

end

@testset "2-atom Coulomb: r=2, q1=q2=-1" begin
    natoms = 2

    atomtype = ProtoSyn.Forcefields.AtomType("X", -1, -1, 1.0, 0.0)
    atoms = [
        ProtoSyn.Forcefields.Atom(-1.0, atomtype),
        ProtoSyn.Forcefields.Atom(-1.0, atomtype)
    ]
    
    state = ProtoSyn.State(natoms)
    state.forces = zeros(natoms, 3)
    state.energy = ProtoSyn.Energy()

    state.coords[2,1] = 2

    @test ProtoSyn.eval!(state, atoms, false) ≈ 0.5
    @test ProtoSyn.eval!(state, atoms, true) ≈ 0.5
    @test state.forces[1,1] < 0.0 && state.forces[2,1] > 0.0

end

@testset "2-atom Coulomb: r=1, q1=1, q2=-1" begin
    natoms = 2

    atomtype = ProtoSyn.Forcefields.AtomType("X", -1, -1, 1.0, 0.0)
    atoms = [
        ProtoSyn.Forcefields.Atom(-1.0, atomtype),
        ProtoSyn.Forcefields.Atom( 1.0, atomtype)
    ]
    
    state = ProtoSyn.State(natoms)
    state.forces = zeros(natoms, 3)
    state.energy = ProtoSyn.Energy()

    state.coords[2,1] = 1.0

    @test ProtoSyn.eval!(state, atoms, false) ≈ -1.0
    @test ProtoSyn.eval!(state, atoms, true) ≈ -1.0
    @test state.forces[1,1] > 0.0 && state.forces[2,1] < 0.0

end