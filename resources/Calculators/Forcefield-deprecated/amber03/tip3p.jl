let

using DelimitedFiles: readdlm
using ProtoSyn.Calculators.Forcefield

# instantiate topology
tip3p = Topology(3)

# define atom types for oxygen and hydrogen
HW = AtomType("HW", 1,  1.008, 0.00000e+00, 0.00000e+00)
OW = AtomType("OW", 8, 16.000, 3.15061e-01, 6.36386e-01)

# bond and angles types
bOH = HarmonicBondType(genkey("OW", "HW"), 0.09572, 502416.0)
aHOH = HarmonicAngleType(genkey("HW", "OW", "HW"), deg2rad(104.52), 628.02)

# fill topology
push!(tip3p, ProtoSyn.Calculators.Forcefield.Atom(-0.834, OW, [2, 3, -1]))
push!(tip3p, ProtoSyn.Calculators.Forcefield.Atom( 0.417, HW,    [3, -1]))
push!(tip3p, ProtoSyn.Calculators.Forcefield.Atom( 0.417, HW,       [-1]))
push!(tip3p, HarmonicAngle(2, 1, 3, aHOH))
push!(tip3p, HarmonicBond(1, 2, bOH))
push!(tip3p, HarmonicBond(1, 3, bOH))


# 216H2O,WATJP01,SPC216,SPC-MODEL,300K,BOX(M)=1.86206NM,WFVG,MAR. 1984
# OW1, HW1, HW2
#1.86206 1.86206 1.86206

# data = readdlm("spc216.xyz"; comments=true) #, ' ', Float64)#, '\n', comments=true)

# box = data[1,:]
# xyz = data[2:end,:]
# box, xyz

tip3p

end
