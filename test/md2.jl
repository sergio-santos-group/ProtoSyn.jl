module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Drivers
using .ProtoSyn.Calculators.Forcefield

# construct model
aa = include(joinpath(ProtoSyn.Peptides.resource_dir, "aminoacids.jl"))
mol, state = ProtoSyn.Peptides.build("A"^20, aa)


# Assign forcefield
ff = include(joinpath(ProtoSyn.Calculators.Forcefield.resource_dir, "amber99/forcefield.jl"))
mp = ProtoSyn.Calculators.Forcefield.loadmap("amber99/aminoacids.yml");
top = ProtoSyn.Calculators.Forcefield.gentop(mol, ff, mp)

# Define custom evaluator
ffeval = ProtoSyn.Calculators.Forcefield.aggregate(top)

function custom_eval(state::State, do_forces::Bool)
    rij = state.xyz[10,:] - state.xyz[20,:]
    dij = sqrt(dot(rij,rij))

    penalty = 100*(12.3 - dij)^2
    return penalty + ffeval(state, do_forces)

end

thermostat = ProtoSyn.Drivers.VRescaleThermostat(
    300.0,  # temperature
)

md = ProtoSyn.Drivers.MolecularDynamics(
    eval! = ProtoSyn.Calculators.Forcefield.aggregate(top),
    thermostat! = thermostat,
    max_steps = 100,
    masses = [at.type.m for at in top.components[:Atom].items],
    timestep = 0.002,

    remove_com_mode = :angular,
    remove_com_freq = 100,

) 




fout = open("traj.pdb", "w")

dstat = md(state) do st,ds
    if ds.step%100 == 0
        energy = state.energy
        println("$(ds.step): E=$(energy.total) temp=$(ds.temperature) t0=$(thermostat.temperature)")
        write(fout, mol, st)
    end
end

close(fout)

end