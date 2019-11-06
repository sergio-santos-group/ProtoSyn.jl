module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Calculators.Forcefield
using .ProtoSyn.Drivers




# construct model
aa = include(joinpath(ProtoSyn.Peptides.resource_dir, "aminoacids.jl"))
mol, state = ProtoSyn.Peptides.build("A"^20, aa)



# Assign forcefield
ff = include(joinpath(ProtoSyn.Calculators.Forcefield.resource_dir, "amber03/forcefield2.jl"))
mp = ProtoSyn.Calculators.Forcefield.loadmap("amber03/aminoacids.yml");
top = ProtoSyn.Calculators.Forcefield.assign(mol, ff, mp)


# configure nonbonded evaluation
top[:Atom].config = ProtoSyn.Calculators.Forcefield.PlainCutoff(2.0)


# specify which components are to be calculated
ev = ProtoSyn.Calculators.Forcefield.aggregate(
    top[:Atom],
    top[:HarmonicBond],
    top[:HarmonicAngle]
)


sd = ProtoSyn.Drivers.SteepestDescent(
    eval! = ProtoSyn.Calculators.Forcefield.aggregate(top),
    max_steps = 100,
    pairlist_freq = 0,
    force_tolerance = 0.1,
    max_displacement = 0.1
)


fout = open("../tmp.pdb", "w")
# write(fout, mol, state)

# Version 1 - specify callback using "do" syntax
dstat2 = sd(state) do st, ds
    if ds.step%10 == 0
        energy = state.energy
        println("$(ds.step): E=$(energy.total) Fmax=$(ds.max_force), stepsize=$(ds.stepsize)")
        # write(fout, mol, st)
    end
end


# # Version 2 - specify callback as first argument
# function cb(st, ds)
#     if ds.step%10 == 0
#         energy = state.energy
#         println("$(ds.step): E=$(energy.total) Fmax=$(ds.max_force), stepsize=$(ds.stepsize)")
#         # write(fout, mol, st)
#     end
# end

# dstat = sd(cb, state)

# # Version 3 - no callback
# dstat = sd(state)


write(fout, mol, state)
close(fout)





end

