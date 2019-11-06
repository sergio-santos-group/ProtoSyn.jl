module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Drivers
using .ProtoSyn.Calculators.Forcefield

using DelimitedFiles: readdlm

# top_soluto
# top_solvent

# top_total = top_soluto + 1000*top_solvent

ffresources = ProtoSyn.Calculators.Forcefield.resource_dir

tip3p0  = include(joinpath(ffresources, "amber03/tip3p.jl"))
spc216 = readdlm(joinpath(ffresources, "amber03/spc216.xyz"); comments=true)

box = spc216[1]
coords = spc216[2:end, :]
coords = [coords; [0.0 box 0.0] .+ coords]

natoms = size(coords, 1)

state = State(natoms)
state.coords = copy(coords)

cutoff = 1.0
state.pairlist = ProtoSyn.VerletList(natoms)
state.pairlist.cutoff = cutoff

tip3p = (2 * 216) * tip3p0

tip3p.components[:Atom].config = ProtoSyn.Calculators.Forcefield.PlainCutoff(cutoff)

const Δt = 0.002

thermostat = ProtoSyn.Drivers.VRescaleThermostat(300.0)

ev = ProtoSyn.Calculators.Forcefield.aggregate(tip3p)
#     tip3p.components[:Atom],
#     tip3p.components[:HarmonicBond],
#     tip3p.components[:HarmonicAngle],
# )

md = ProtoSyn.Drivers.MolecularDynamics(
    eval! = ev, #ProtoSyn.Calculators.Forcefield.aggregate(tip3p),
    max_steps = 5_000,
    masses = [at.type.m for at in tip3p.components[:Atom].items],
    timestep = Δt,
    thermostat! = thermostat,
    pairlist_freq = 20,
    remove_com_mode = :angular,
    remove_com_freq = 100,

) 

# md.max_steps = 10
# @time md(state)

# md.max_steps = 500
# @time md(state)




# tip3p.components[:Atom].config = nothing
# md.max_steps = 10
# @time md(state)

# md.max_steps = 500
# @time md(state)


fout = open("../tmp.pdb", "w")
write(fout, tip3p, state)

dstat = md(state) do st,ds
    if ds.step%50 == 0
        energy = state.energy
        println("$(ds.step): E=$(energy.total) temp=$(ds.temperature) t0=$(thermostat.temperature)")
        write(fout, tip3p, st)
        # write(fout, st)
    end
end

close(fout)

#md(state)

end