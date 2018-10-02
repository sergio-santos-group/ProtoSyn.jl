var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "(Image: header)"
},

{
    "location": "index.html#O-PAI-NATAL-MORREU-ENGASGADO-1",
    "page": "Home",
    "title": "O PAI NATAL MORREU ENGASGADO",
    "category": "section",
    "text": "PAI NATAL MORRE ENGASGADO A COMER COELHO DA PÁSCOA À CAÇADOR!!(Image: Santa)CurrentModule = Forcefield<!– ```@docs evalbond!javascript function fancyAlert(arg) {   if(arg) {     $.facebox({div:\'#foo\'})   } } ```"
},

{
    "location": "forcefield.html#",
    "page": "Forcefield",
    "title": "Forcefield",
    "category": "page",
    "text": ""
},

{
    "location": "forcefield.html#Force-field-1",
    "page": "Forcefield",
    "title": "Force-field",
    "category": "section",
    "text": "CurrentModule = Forcefield"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.HarmonicBond",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.HarmonicBond",
    "category": "type",
    "text": "HarmonicBond(a1::Int64, a2::Int64, k::Float64, b0::Float64)\n\nHarmonic Bond of the form\n\nE(r_ab) = frac12k_ab(r_ab - b_0)^2\n\nwhere\n\nr_ab = vecr_ab = vecr_b - vecr_a\n\nArguments\n\na1::Int64, a2::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹ nm⁻²).\nb0::Float64: equilibrium bond length (nm).\n\nExamples\n\njulia> Forcefield.HarmonicBond(1, 2, 2500, 0.19)\nForcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19)\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.HarmonicAngle",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.HarmonicAngle",
    "category": "type",
    "text": "HarmonicAngle(a1::Int64, a2::Int64, a3::Int64, k::Float64, θ::Float64)\n\nHarmonic Angle of the form\n\nE(θ_abc)=frac12k_abc(theta_abc-theta)^2\n\nArguments\n\na1::Int64, a2::Int64, a3::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹ rad⁻²).\nθ::Float64: equilibrium angle value (rad).\n\nExamples\n\njulia> Forcefield.HarmonicAngle(1, 2, 3, 670.0, 1.92)\nForcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92)\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.DihedralCos",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.DihedralCos",
    "category": "type",
    "text": "DihedralCos(a1::Int64, a2::Int64, a3::Int64, a4::Int64, k::Float64, θ::Float64, mult::Float64)\n\nPeriodic Dihedral of the form\n\nE(phi_abcd)=K_phi(1+cos(nphi_abcd-phi))\n\nArguments\n\na1::Int64, a2::Int64, a3::Int64, a4::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹).\nθ::Float64: equilibrium angle value (deg).\nmult::Float64: multiplicity.\n\nExamples\n\njulia> Forcefield.DihedralCos(1, 2, 3, 4, 10.46, 180.0, 2.0)\nForcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0)\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Atom",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Atom",
    "category": "type",
    "text": "Atom(name::String, σ::Float64, ϵ::Float64, q::Float64, excls::Array{Int64, 1}, pairs::Array{Int64, 1})\n\nDefine an atom.  σ, ϵ and q describe the non-bonded interactions between atoms:\n\nThe Lennard-Jones interaction is in the form:\n\nE(r_ab) = 4ϵ_ableft ( (fracσ_abr_ab)^12-(fracσ_abr_ab)^6right )\n\nwhere the Lorentz-Berthelot rule is applied. σ is the arithmetic average and ϵ is the geometric average:\n\nσ_ab=fracσ_a+σ_b2\n\nϵ_ab=sqrt(ϵ_aϵ_b)\n\nFor this reason, σ and ϵ are applied here in the reduced form: fracσ2 and sqrtϵ.\n\nThe Coulomb interation is in the form:\n\nE(r_ab)=k_ϵfracq_aq_br_ab^2\n\nwhere\n\nk_ϵ=frac14πϵ_0=138935485kJnmmol¹e¹\n\nFor this reason, q is applied here in the reduced form: qtimes sqrtk_ϵ\n\nExclusion list contains all atom indices who are excluded from non-bonded interactions (i.e. are at 3 or less connections from this atom - includes pairs). Pair list contains atoms that are at 3 connections from this atom, and are involved in 1-4 interactions (and have a different combination rule as a result).\n\nArguments\n\nname::String: Atom name (example: \"C\", \"H\", etc).\nσ::Float64: half the finite distance at which the inter-particle potential is zero (nm).\nϵ::Float64: square root of the depth of the potential well (kJ mol⁻¹).\nq::Float64: adjusted atom charge (eletron).\nexcls::Array{Int64, 1}: exclusion list (as global atom indices).\npairs::Array{Int64, 1}: pair list containing atoms that interfere in 1-4 interations (as global atom indices)\n\nExamples\n\njulia> Forcefield.Atom(\"N\", 0.325, 0.711, 0.0017, [0, 1, 2, 3, 4, 5], [4, 5])\nForcefield.Atom(name=\"N\", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5])\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Topology",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Topology",
    "category": "type",
    "text": "Topology(atoms::Array{Atom}, bonds::Array{HarmonicBond}, angles::Array{HarmonicAngle}, dihedralsCos::Array{DihedralCos})\n\nGather all topology components.\n\nArguments\n\natoms::Array{Atoms}\nbonds::Array{HarmonicBond}\nangles::Array{HarmonicAngle}\ndihedralsCos::Array{DihedralCos}\n\nExamples\n\njulia> Forcefield.Forcefield(atoms, bonds, angles, dihedrals)\nForcefield.Topology(\n atoms=ProtoSyn.Forcefield.Atom[Forcefield.Atom(name=\"N\", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5]), ...],\n bonds=ProtoSyn.Forcefield.HarmonicBond[Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19), ...],\n angles=ProtoSyn.Forcefield.HarmonicAngle[Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92), ...],\n dihedralsCos=ProtoSyn.Forcefield.DihedralCos[Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0), ...])\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Energy",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Energy",
    "category": "type",
    "text": "DihedralCos(eBond::Float64, eAngle::Float64, eDihedral::Float64, eLJ::Float64, eLJ14::Float64, eCoulomb::Float64, eCoulomb14::Float64, eTotal::Float64)\n\nEnergy components.\n\nArguments\n\na1::Int64, a2::Int64, a3::Int64, a4::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹).\nθ::Float64: equilibrium angle value (deg).\nmult::Float64: multiplicity.\n\nExamples\n\njulia> Forcefield.Energy()\nForcefield.Energy(\n eBond=0.0,\n eAngle=0.0,\n eDihedral=0.0,\n eLJ=0.0,\n eLJ14=0.0,\n eCoulomb=0.0,\n eCoulomb14=0.0,\n eTotal=0.0\n)\n\njulia> Forcefield.Energy(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)\nForcefield.Energy(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#Components-1",
    "page": "Forcefield",
    "title": "Components",
    "category": "section",
    "text": "This section provides a description on all the individual components that comprise the Forcefield.HarmonicBond\nHarmonicAngle\nDihedralCos\nAtom\nTopology\nEnergy"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.evalbond!",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.evalbond!",
    "category": "function",
    "text": "evalbond!(bonds::Array{Forcefield.HarmonicBond}, state::Common.State, do_forces::Bool = false) -> Float64\n\nEvaluate an array of Bonds using the given state.\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#Evaluators-1",
    "page": "Forcefield",
    "title": "Evaluators",
    "category": "section",
    "text": "This section provides a description on how to use the previously documented Components to calculate their energy.evalbond!"
},

]}
