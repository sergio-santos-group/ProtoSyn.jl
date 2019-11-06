module Tst

include("../src/ProtoSyn.jl")

using .ProtoSyn
#using .ProtoSyn.Drivers
#using .ProtoSyn.Calculators.Forcefield




# construct model
aa = include(joinpath(ProtoSyn.Peptides.resource_dir, "aminoacids.jl"))

mol, state = ProtoSyn.Peptides.build("A"^5, aa)
#tip3p = include(joinpath(ProtoSyn.Calculators.Forcefield.resource_dir, "amber03/tip3p.jl"))

# identify phi,psi pairs
dihedrals = Vector{NTuple{2,Dihedral}}()
let
    ar = zeros(Int, 5)
    for residue in ProtoSyn.iterbyresidue(mol)
        ProtoSyn.cproduct(residue, ("-C","N","CA","C","+N"), 0, ar) do indices
            d1 = Dihedral(indices[1:4]..., ProtoSyn.Peptides.DihedralTypes.phi)
            d2 = Dihedral(indices[2:5]..., ProtoSyn.Peptides.DihedralTypes.psi)
            push!(dihedrals, (d1,d2))
        end
    end
end

# # Assign forcefield
# ff = include(joinpath(ProtoSyn.Calculators.Forcefield.resource_dir, "amber03/forcefield2.jl"))
# mp = ProtoSyn.Calculators.Forcefield.loadmap("amber03/aminoacids.yml");
# top = ProtoSyn.Calculators.Forcefield.gentop(mol, ff, mp)


# top.components[:Atom].config = ProtoSyn.Calculators.Forcefield.PlainCutoff(1.2)

# const Δt = 0.002

# # masses = top.
# thermostat = ProtoSyn.Drivers.VRescaleThermostat(
#     300.0,  # temperature
#     #1.0,    # coupling
#     #Δt      # integration timestep
# )

# md = ProtoSyn.Drivers.MolecularDynamics(
#     eval! = ProtoSyn.Calculators.Forcefield.aggregate(top),
#     max_steps = 100,
#     masses = [at.type.m for at in top.components[:Atom].items],
#     timestep = Δt,
#     thermostat! = thermostat,

#     remove_com_mode = :angular,
#     remove_com_freq = 100,

# ) 




# fout = open("../tmp.pdb", "w")

# dstat = md(state) do st,ds
#     if ds.step%100 == 0
#         energy = state.energy
#         println("$(ds.step): E=$(energy.total) temp=$(ds.temperature) t0=$(thermostat.temperature)")
#         write(fout, mol, st)
#     end
# end

# close(fout)

end