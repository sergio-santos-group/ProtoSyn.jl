# Load resources
include("../src/ProtoSyn.jl")

using .ProtoSyn
using .ProtoSyn.Peptides
using .ProtoSyn
using .ProtoSyn.Units
println(" | ProtoSyn loaded successfully.")

using Printf
using Random
using PyCall

seed = convert(Int, floor(rand() * 100000))
# seed = 46747
Random.seed!(seed)

ase = pyimport("ase")
np  = pyimport("numpy")
aimnet = pyimport("aimnet")
model_smd = aimnet.load_AIMNetSMD_ens().cuda()
calc_smd = aimnet.AIMNetCalculator(model_smd)
println(" | AIMNet model loaded successfully.")

# Start algorithm
res_lib = grammar();
pose = Peptides.build(res_lib, seq"MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAY");
# Peptides.setss!(pose, SecondaryStructure[:helix], rid"1:20");
# Peptides.setss!(pose, SecondaryStructure[:helix], rid"27:44");

s = get_ani_species(pose)
c = get_cartesian_matrix(pose)
atoms = ase.Atoms(positions=c, numbers=s)
atoms.set_calculator(calc_smd)
aimnet_energy = atoms.get_potential_energy()[1]
println(e)