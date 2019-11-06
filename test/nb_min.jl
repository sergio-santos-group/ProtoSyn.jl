module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Forcefields
using .ProtoSyn.Drivers.SteepestDescent
const SD = ProtoSyn.Drivers.SteepestDescent

const natoms = 4
rmin = 2.0
ϵ = 1.0
q = 0.0
m = 1.0
z = 1
atomtype = ProtoSyn.Forcefields.AtomType("X", z, m, rmin / (2.0^(1/6)), ϵ)

atoms = ProtoSyn.Forcefields.AtomContainer([
    ProtoSyn.Forcefields.Atom(q, atomtype)
    for i = 1:natoms
])

state = ProtoSyn.State(natoms)
state.coords =  rand(natoms, 3)
state.forces = zeros(natoms, 3)
state.energy = ProtoSyn.Energy()

ev = ProtoSyn.Forcefields.aggregate(atoms)

dconf = SD.DriverConfig(
    eval! = ev,
    max_steps=20000,
    nblist_freq=0,
    force_tolerance=0.001,
    max_displacement=0.1,
)

fout = open("../tmp.pdb", "w")

write(fout, state)

dstat = ProtoSyn.run!(state, dconf) do st, ds
    #if ds.step%20 == 0
        energy = state.energy
        println("$(ds.step): E=$(energy.total) Fmax=$(ds.max_force), stepsize=$(ds.stepsize)")
        println("   $(energy.components)")
        write(fout, st)
    #end
end

close(fout)


end

