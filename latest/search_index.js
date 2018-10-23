var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "(Image: header)"
},

{
    "location": "index.html#PROTOSYN-1",
    "page": "Home",
    "title": "PROTOSYN",
    "category": "section",
    "text": "WORK IN PROGRESS!"
},

{
    "location": "common.html#",
    "page": "Common",
    "title": "Common",
    "category": "page",
    "text": ""
},

{
    "location": "common.html#Common-1",
    "page": "Common",
    "title": "Common",
    "category": "section",
    "text": "CurrentModule = Common"
},

{
    "location": "common.html#Components-1",
    "page": "Common",
    "title": "Components",
    "category": "section",
    "text": "This section provides a description of the Common components."
},

{
    "location": "common.html#ProtoSyn.Common.NullEnergy",
    "page": "Common",
    "title": "ProtoSyn.Common.NullEnergy",
    "category": "type",
    "text": "NullEnergy()\n\nEmpty placeholder energy container.\n\nExamples\n\njulia> Common.NullEnergy()\nNull\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.Energy",
    "page": "Common",
    "title": "ProtoSyn.Common.Energy",
    "category": "type",
    "text": "Energy(eTotal::Float64)\n\nSimple energy container.\n\nExamples\n\njulia> Common.Energy(1.15)\nEnergy(eTotal=1.15)\n\njulia> Common.Energy()\nEnergy(eTotal=0.0)\n\n\n\n\n\n"
},

{
    "location": "common.html#Energy-1",
    "page": "Common",
    "title": "Energy",
    "category": "section",
    "text": "Contains common and simple energy representations. More specific energy structures can be used from other modules such as AmberNullEnergy\nEnergy"
},

{
    "location": "common.html#ProtoSyn.Common.SSTYPE",
    "page": "Common",
    "title": "ProtoSyn.Common.SSTYPE",
    "category": "type",
    "text": "SSTYPE\n\nEnum: holds information regarding the secondary structure of each residue in the simulation.\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.Residue",
    "page": "Common",
    "title": "ProtoSyn.Common.Residue",
    "category": "type",
    "text": "Residue(atoms::Array{Int64, 1}, next::Union{Residue, Int64, Nothing}, name::String, ss::SecondaryStructureType)\n\nDefine a residue as part of the system.\n\nArguments\n\natoms::Vector{Int64}: list of atom global indices in this residue.\nnext::Union{Residue, Int64, Nothing}: Next residue in the system, attached to this. Is preferably a Residue instance, but can in certain cases be the index in a Residue list or empty (Nothing).\nname::String: Name of the residue. If the residue describes an aminoacid, the correspondent letter is suggested.\nss::SSTYPE: The intitial secondary structure type SSTYPE of this residue.\n\nExamples\n\njulia> Common.Residue([1, 2, 3, 4], Common.Residue([5, 6, 7, 8], nothing, \"V\", Common.coil), \"E\", Common.alpha)\nCommon.Residue(atoms=[1, 2, 3, 4], next=V, name=E, ss=alpha)\n\n\n\n\n\n"
},

{
    "location": "common.html#Residue-1",
    "page": "Common",
    "title": "Residue",
    "category": "section",
    "text": "In ProtoSyn, a Residue object is a collection of atoms (normally an aminoacid) that are identified by a name, a secondary structure (SSTYPE) and are part of a continuous tree of other Residues (have a next Residue).SSTYPE\nResidue"
},

{
    "location": "common.html#ProtoSyn.Common.DIHEDRALTYPE",
    "page": "Common",
    "title": "ProtoSyn.Common.DIHEDRALTYPE",
    "category": "type",
    "text": "DIHEDRALTYPE\n\nEnum: holds information regarding the dihedral type. ```\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.Dihedral",
    "page": "Common",
    "title": "ProtoSyn.Common.Dihedral",
    "category": "type",
    "text": "Dihedral(a1::Int64, a2::Int64, a3::Int64, a4::Int64, movable::Array{Int64, 1}, residue::Union{Common.Residue, Int64}, dtype::DIHEDRALTYPE)\n\nDefine a dihedral.\n\nArguments\n\na1::Int64, a2::Int64, a3::Int64, a4::Int64: global atom indices.\nmovable::Array{Int64, 1}: List of global atom indices that will be moved during the dihedral movement in this residue.\nresidue::{Residue}: Residue object that this dihedral belongs to.\ndtype::DIHEDRALTYPE: DIHEDRALTYPE.\n\nExamples\n\njulia> Dihedral(2, 3, 5, 7, [5, 6], Common.Residue([1, 2, 3, 4, 5, 6], (...), \"A\", coil), phi)\nDihedral(a1=2, a2=3, a3=5, a4=7, movable=[5, 6], residue=Common.Residue(atoms=[1, 2, 3, 4, 5, 6], next=V, name=A), type=phi)\n\nSee also: Mutators.DihedralMutator\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.rotate_dihedral!",
    "page": "Common",
    "title": "ProtoSyn.Common.rotate_dihedral!",
    "category": "function",
    "text": "rotate_dihedral!(xyz::Array{Float64, 2}, dihedral::Dihedral, angle::Float64)\n\nPerform a dihedral movement, adding the provided angle (in radians). If the dihedral.dtype is \"PHI\" or \"PSI\" the dihedral.residue.next is also rotated and this is propagated recursively until the end of the molecule. \n\nExamples\n\njulia> Mutators.Dihedral.rotate_dihedral!(state.xyz, dihedral, π/2)\n\nSee also: Dihedral\n\n\n\n\n\nrotate_dihedral!(xyz::Array{Float64, 2}, a2::Int64, a3::Int64, angle::Float64, dtype::DIHEDRALTYPE, movable::Vector{Int64}[, residue::Union{Residue, Nothing} = nothing])\n\nBase dihedral movement function. Especifies all arguments used in dihedral rotation movement. \n\nExamples\n\njulia> Mutators.Dihedral.rotate_dihedral!(xyz, dihedral.a2, dihedral.a3, π/2, dihedral.dtype, dihedral.movable, dihedral.residue)\n\nSee also: Aux.rotation_matrix_from_axis_angle Dihedral Crankshaft\n\n\n\n\n\n"
},

{
    "location": "common.html#Dihedral-1",
    "page": "Common",
    "title": "Dihedral",
    "category": "section",
    "text": "A Dihedral is a collection of 4 atoms that define a dihedral in the simulated molecule. Many Mutators operate over this dihedrals, changing them in order to explore the conformational space of the system. A Dihedral is part of a Residue and has a defined DIHEDRALTYPE.DIHEDRALTYPE\nDihedral\nrotate_dihedral!"
},

{
    "location": "common.html#ProtoSyn.Common.AtomMetadata",
    "page": "Common",
    "title": "ProtoSyn.Common.AtomMetadata",
    "category": "type",
    "text": "AtomMetadata(name::String[, elem::String = name, res_num::Int64 = 1, res_name::String = \"UNK\", chain_id::Union{String, Nothing} = nothing, connects::Union{Vector{Int64}, Nothing} = nothing])\n\nDefine an atom metadata, containing extra information pertaining the State.\n\nArguments\n\nname::String: Name of the atom.\nelem::String: (Optional) Element of the atom (Default: name).\nres_num::Int64: (Optional) Number of the residue this atom belongs to (Default: 1).\nres_name::Union{String, Nothing}: (Optional) Name of the residue this atom belongs to (Default: \"UNK\").\nchain_id::String: (Optional) Name of the chain that contains the residue this atom belongs to (Default: nothing).\nconnects::Union{Vector{Int64}, Nothing}: (Optional) List of global atom indices that this atom is connected to (Default: nothing). \n\nExamples\n\njulia> AtomMetadata(\"H1\", \"H\", 2, \"VAL\", \"A\", [4])\nAtomMetadata(name=H1, elem=H, res_num=2, res_name=VAL, chain_id=A, connects=[4])\n\njulia> AtomMetadata(\"H1\")\nAtomMetadata(name=H1, elem=H1, res_num=1, res_name=UNK, chain_id=nothing, connects=nothing)\n\nSee also: iter\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.iter",
    "page": "Common",
    "title": "ProtoSyn.Common.iter",
    "category": "function",
    "text": "iter(data::Vector{AtomMetadata}; property::Symbol = :res_num)\n\nIterate over an array of AtomMetadata objects, grouping them based on property (Default: :res_num)\n\nExamples\n\njulia> for residue in iter(state.metadata)\n    println(residue)\nend\n[AtomMetadata(name=H1, elem=H1, res_num=1, res_name=UNK, chain_id=nothing, connects=nothing), AtomMetadata(name=H2, elem=H2, res_num=1, res_name=UNK, chain_id=nothing, connects=nothing)]\n[AtomMetadata(name=H3, elem=H3, res_num=2, res_name=UNK, chain_id=nothing, connects=nothing), AtomMetadata(name=H4, elem=H4, res_num=2, res_name=UNK, chain_id=nothing, connects=nothing)]\n\nSee also: AtomMetadata\n\n\n\n\n\n"
},

{
    "location": "common.html#Metadata-1",
    "page": "Common",
    "title": "Metadata",
    "category": "section",
    "text": "ProtoSyn Metadata defines additional information of the system that is not necessarily necessary for the basic functions of the library, but allows for a better representation of the system. The AtomMetadata structure holds information related to each atom in the system, such as its element, connections, etcAtomMetadata\niter"
},

{
    "location": "common.html#ProtoSyn.Common.State",
    "page": "Common",
    "title": "ProtoSyn.Common.State",
    "category": "type",
    "text": "State(size::Int64, energy::AbstractEnergy, xyz::Array{Float64, 2}, forces::Array{Float64, 2}, metadata::Vector{AtomMetadata})\n\nDefine the current state of the system, containing the atoms positions, energy and forces applied. If only size::Int64 is provided, an empty State with the given size is created with zeros.\n\nArguments\n\nsize::Int64: Atom count in system.\nenergy::AbstractEnergy: Current energy of the system (kJ mol⁻¹).\nxyz::Array{Float64, 2}: Atom positions in 3 dimensions.\nforces::Array{Float64, 2}: Forces applied in each dimension to each atom (kJ mol⁻¹ nm⁻¹)\nmetadata::Vector{AtomMetadata}: List of atom names.\n\nExamples\n\njulia> Common.State(3)\nCommon.State(size=3, energy=Null, xyz=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], forces=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], metadata=(...))\n\njulia> Common.State(2, Common.NullEnergy(), [1.1 1.1 1.1; 2.2 2.2 2.2], zeros(2, 3), [AtomMetadata(...), AtomMetadata(...), ...])\nCommon.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], metadata=(...))\n\nSee also: AtomMetadata\n\n\n\n\n\n"
},

{
    "location": "common.html#State-1",
    "page": "Common",
    "title": "State",
    "category": "section",
    "text": "The system state holds information about the current coordinates, energy and forces, aswell as any additional metadata. If iterated over, it returns atom by atom position and metadata.State"
},

{
    "location": "common.html#ProtoSyn.Common.CallbackObject",
    "page": "Common",
    "title": "ProtoSyn.Common.CallbackObject",
    "category": "type",
    "text": "CallbackObject(freq::Int64, callback::Function)\n\nDefine the callback function parameters.\n\nArguments\n\nfreq: Frequency (in steps) that the callback function is called.\ncallback: Actual callback function. This function should have the following signature:\n\ncallback(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, args...)\n\nExamples\n\njulia> Common.CallbackObject(100, Print.as_xyz)\nCallbackObject(freq=100, callback=Print.as_xyz)\n\njulia> Common.CallbackObject(Print.as_xyz)\nCallbackObject(freq=1, callback=Print.as_xyz)\n\nSee also: Print.as_xyz @cbcall\n\n\n\n\n\n"
},

{
    "location": "common.html#Callback-1",
    "page": "Common",
    "title": "Callback",
    "category": "section",
    "text": "The CallbackObject allows for independent calls to various functions with individual frequency of output. CallbackObject"
},

{
    "location": "common.html#ProtoSyn.Common.load_from_gro",
    "page": "Common",
    "title": "ProtoSyn.Common.load_from_gro",
    "category": "function",
    "text": "load_from_gro(i_file::String)::Common.State\n\nReturn a new Common.State by loading the atom positions and names from the input .gro file. As a default, state.energy is NullEnergy and state.forces are set to zero.\n\nExamples\n\njulia> Common.load_from_gro(\"molecule.gro\")\nCommon.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], atnames=[\"C\", \"O\"])\n\nSee also: load_from_pdb\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.load_from_pdb",
    "page": "Common",
    "title": "ProtoSyn.Common.load_from_pdb",
    "category": "function",
    "text": "load_from_pdb(i_file::String)::Common.State\n\nReturn a new Common.State by loading the atom positions and names from the input .pdb file. As a default, state.energy is NullEnergy and state.forces are set to zero.\n\nExamples\n\njulia> Common.load_from_pdb(\"molecule.pdb\")\nCommon.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], metadata=(...))\n\nSee also: load_from_gro\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.load_topology",
    "page": "Common",
    "title": "ProtoSyn.Common.load_topology",
    "category": "function",
    "text": "load_topology(p::Dict{String, Any})\n\nParse a dictionary containing the dihedral and residue topology. Return a Dihedral array and a Residue array.\n\nExamples\n\njulia> Mutators.Diehdral.load_topology(p)\n(ProtoSyn.Mutators.Dihedral.NewDihedral[...], ProtoSyn.Common.Residue[...])\n\nSee also: Aux.read_JSON\n\n\n\n\n\n"
},

{
    "location": "common.html#Loaders-1",
    "page": "Common",
    "title": "Loaders",
    "category": "section",
    "text": "This section provides a description on how to load a new State, Residue and Dihedral arrays.load_from_gro\nload_from_pdb\nload_topology"
},

{
    "location": "common.html#ProtoSyn.Common.apply_initial_conf!",
    "page": "Common",
    "title": "ProtoSyn.Common.apply_initial_conf!",
    "category": "function",
    "text": "apply_initial_conf!(state::State, dihedrals::Vector{Dihedral})\n\nApply predefined angles to all dihedrals defined in dihedrals, based on the Dihedral.residue.ss, changing the State.xyz to apply the secondary structure. The applied angles (in degrees) are the following:\n\nBeta sheet: PHI = -139.0 | PSI = 135.0 Alpha helix: PHI = -57.0  | PSI = -47.0\n\nExamples\n\njulia> Common.apply_initial_conf(state, dihedrals)\n\n\n\n\n\n"
},

{
    "location": "common.html#Conformation-Generators-1",
    "page": "Common",
    "title": "Conformation Generators",
    "category": "section",
    "text": "Conformation generators are responsible to change the system State in a defined way.apply_initial_conf!"
},

{
    "location": "common.html#ProtoSyn.Common.@cbcall",
    "page": "Common",
    "title": "ProtoSyn.Common.@cbcall",
    "category": "macro",
    "text": "@Common.cbcall callbacks::Tuple{CallbackObject, N} step::Int64 Vararg::Any\n\n(Macro) Call the CallbackObject.function of each CallbackObject in the callbacks Tuple{CallbackObject, N} independently. Each CallbackObject.function is ran depending on the defined CallbackObject.freq and the given step. Vararg holds all the arguments necessary to run the callback function itself.\n\nExamples\n\njulia> @Common.cbcall (callback_object1, callback_object2) 1\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.@callback",
    "page": "Common",
    "title": "ProtoSyn.Common.@callback",
    "category": "macro",
    "text": "@Common.callback f::Function freq::Int64\n\n(Macro) Create a CallbackObject with the given function f and output frequency freq. \n\nExamples\n\njulia> @Common.callback 10 my_callback1\nCallbackObject(freq=10, callback=my_callback1)\n\njulia> @Common.callback my_callback1\nCallbackObject(freq=1, callback=my_callback1)\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.@faggregator",
    "page": "Common",
    "title": "ProtoSyn.Common.@faggregator",
    "category": "macro",
    "text": "@Common.faggregator name::String f::function Vararg::Any\n\n(Macro) Aggregate multiple functions f in a single variable name. Vararg contains the arguments used by function f\n\nExamples\n\njulia> @faggregator myeval f top1\n@faggregator myeval g top2\n@faggregator myeval h top3\n\nenergy = myevalf(state, false)\n\n\n\n\n\n"
},

{
    "location": "common.html#Macros-1",
    "page": "Common",
    "title": "Macros",
    "category": "section",
    "text": "Auxiliary functions that help speed up the system\'s performance.@cbcall\n@callback\n@faggregator"
},

{
    "location": "forcefield.html#",
    "page": "Forcefield",
    "title": "Forcefield",
    "category": "page",
    "text": ""
},

{
    "location": "forcefield.html#Forcefield-1",
    "page": "Forcefield",
    "title": "Forcefield",
    "category": "section",
    "text": "Currently, ProtoSyn only supports the Amber forcefield.CurrentModule = Forcefield"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Amber.HarmonicBond",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Amber.HarmonicBond",
    "category": "type",
    "text": "HarmonicBond(a1::Int64, a2::Int64, k::Float64, b0::Float64)\n\nHarmonic Bond of the form\n\nE(r_ab) = frac12k_ab(r_ab - b_0)^2\n\nwhere\n\nr_ab = vecr_ab = vecr_b - vecr_a\n\nArguments\n\na1::Int64, a2::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹ nm⁻²).\nb0::Float64: equilibrium bond length (nm).\n\nExamples\n\njulia> Forcefield.HarmonicBond(1, 2, 2500, 0.19)\nForcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19)\n\nSee algo: Amber.evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Amber.HarmonicAngle",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Amber.HarmonicAngle",
    "category": "type",
    "text": "HarmonicAngle(a1::Int64, a2::Int64, a3::Int64, k::Float64, θ::Float64)\n\nHarmonic Angle of the form\n\nE(θ_abc)=frac12k_abc(theta_abc-theta)^2\n\nArguments\n\na1::Int64, a2::Int64, a3::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹ rad⁻²).\nθ::Float64: equilibrium angle value (rad).\n\nExamples\n\njulia> Forcefield.HarmonicAngle(1, 2, 3, 670.0, 1.92)\nForcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92)\n\nSee algo: Amber.evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Amber.DihedralCos",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Amber.DihedralCos",
    "category": "type",
    "text": "DihedralCos(a1::Int64, a2::Int64, a3::Int64, a4::Int64, k::Float64, θ::Float64, mult::Float64)\n\nPeriodic Dihedral of the form\n\nE(phi_abcd)=K_phi(1+cos(nphi_abcd-phi))\n\nArguments\n\na1::Int64, a2::Int64, a3::Int64, a4::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹).\nθ::Float64: equilibrium angle value (rad).\nmult::Float64: multiplicity.\n\nExamples\n\njulia> Forcefield.DihedralCos(1, 2, 3, 4, 10.46, 180.0, 2.0)\nForcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0)\n\nSee algo: Amber.evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Amber.Atom",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Amber.Atom",
    "category": "type",
    "text": "Atom(name::String, σ::Float64, ϵ::Float64, q::Float64, excls::Array{Int64, 1}, pairs::Array{Int64, 1})\n\nDefine an atom.  σ, ϵ and q describe the non-bonded interactions between atoms:\n\nThe Lennard-Jones interaction is in the form:\n\nE(r_ab) = 4ϵ_ableft ( (fracσ_abr_ab)^12-(fracσ_abr_ab)^6right )\n\nwhere the Lorentz-Berthelot rule is applied. σ is the arithmetic average and ϵ is the geometric average:\n\nσ_ab=fracσ_a+σ_b2\n\nϵ_ab=sqrt(ϵ_aϵ_b)\n\nFor this reason, σ and ϵ are applied here in the reduced form: fracσ2 and sqrtϵ.\n\nThe Coulomb interation is in the form:\n\nE(r_ab)=k_ϵfracq_aq_br_ab^2\n\nwhere\n\nk_ϵ=frac14πϵ_0=138935485kJnmmol¹e¹\n\nFor this reason, q is applied here in the reduced form: qtimes sqrtk_ϵ\n\nExclusion list contains all atom indices who are excluded from non-bonded interactions (i.e. are at 3 or less connections from this atom - includes pairs). Pair list contains atoms that are at 3 connections from this atom, and are involved in 1-4 interactions (and have a different combination rule as a result).\n\nArguments\n\nname::String: Atom name (example: \"C\", \"H\", etc).\nσ::Float64: finite distance at which the inter-particle potential is zero (nm).\nϵ::Float64: depth of the potential well (kJ mol⁻¹).\nq::Float64: atom charge (eletron).\nexcls::Array{Int64, 1}: exclusion list (as global atom indices).\npairs::Array{Int64, 1}: pair list containing atoms that interfere in 1-4 interations (as global atom indices)\n\nExamples\n\njulia> Forcefield.Atom(\"N\", 0.325, 0.711, 0.0017, [0, 1, 2, 3, 4, 5], [4, 5])\nForcefield.Atom(name=\"N\", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5])\n\nSee algo: Amber.evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Amber.Topology",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Amber.Topology",
    "category": "type",
    "text": "Topology(atoms::Array{Atom}, bonds::Array{HarmonicBond}, angles::Array{HarmonicAngle}, dihedralsCos::Array{DihedralCos})\n\nGather all topology components.\n\nArguments\n\natoms::Array{Atoms}\nbonds::Array{HarmonicBond}\nangles::Array{HarmonicAngle}\ndihedralsCos::Array{DihedralCos}\n\nExamples\n\njulia> Forcefield.Forcefield(atoms, bonds, angles, dihedrals)\nForcefield.Topology(\n atoms=ProtoSyn.Forcefield.Atom[Forcefield.Atom(name=\"N\", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5]), ...],\n bonds=ProtoSyn.Forcefield.HarmonicBond[Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19), ...],\n angles=ProtoSyn.Forcefield.HarmonicAngle[Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92), ...],\n dihedralsCos=ProtoSyn.Forcefield.DihedralCos[Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0), ...])\n\nSee also: Amber.load_from_json\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Amber.Energy",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Amber.Energy",
    "category": "type",
    "text": "Energy(eBond::Float64, eAngle::Float64, eDihedral::Float64, eLJ::Float64, eLJ14::Float64, eCoulomb::Float64, eCoulomb14::Float64, eTotal::Float64)\n\nEnergy components.\n\nExamples\n\njulia> Forcefield.Energy()\nForcefield.Energy(eBond=0.0, eAngle=0.0, eDihedral=0.0, eLJ=0.0, eLJ14=0.0, eCoulomb=0.0, eCoulomb14=0.0, eTotal=0.0)\n\njulia> Forcefield.Energy(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 2.8)\nForcefield.Energy(eBond=0.1, eAngle=0.2, eDihedral=0.3, eLJ=0.4, eLJ14=0.5, eCoulomb=0.6, eCoulomb14=0.7, eTotal=2.8)\n\nSee also: Amber.evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#Components-1",
    "page": "Forcefield",
    "title": "Components",
    "category": "section",
    "text": "This section provides a description on all the individual components that comprise the Forcefield.Amber.HarmonicBond\nAmber.HarmonicAngle\nAmber.DihedralCos\nAmber.Atom\nAmber.Topology\nAmber.Energy"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Amber.evaluate!",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Amber.evaluate!",
    "category": "function",
    "text": "evaluate!(bonds::Array{Forcefield.HarmonicBond}, state::Common.State[, do_forces::Bool = false])::Float64\n\n\n\n\n\nevaluate!(angles::Array{Forcefield.HarmonicAngle}, state::Common.State, do_forces::Bool = false)::Float64\n\n\n\n\n\nevaluate!(dihedralsCos::Array{Forcefield.DihedralCos}, state::Common.State, do_forces::Bool = false)::Float64\n\n\n\n\n\nevaluate!(atoms::Array{Forcefield.Atom}, state::Common.State, do_forces::Bool = false)::Float64\n\nEvaluate an array of Forcefield.Components using the current Common.State, calculate and update state.energy according to the equations defined in each component. If do_forces flag is set to true, calculate and update state.forces. Non-bonded interactions are only assessed if the distance between atoms is below the defined cut_off value. Return the component energy value (kJ mol⁻¹).\n\nExamples\n\njulia> Forcefield.evaluate!(bonds, state)\n0.500\n\nSee also: evaluate! Amber.HarmonicBond Amber.HarmonicAngle Amber.DihedralCos Amber.Atom\n\n\n\n\n\nevaluate!(topology::Forcefield.Topology, state::Common.State[, cut_off::Float64 = 2.0, do_forces::Bool = false])::Float64\n\nEvaluate the current Common.State energy according to the defined Amber.Topology. If do_forces bool is set to true, calculate and update state.forces. Non-bonded interactions are only assessed if the distance between atoms is below the defined cut_off value. Return state.energy.eTotal value (kJ mol⁻¹).\n\nExamples\n\njulia> Forcefield.evalenergy!(topology, state, cut_off = Inf)\n0.500\n\nSee also: Amber.evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#Evaluators-1",
    "page": "Forcefield",
    "title": "Evaluators",
    "category": "section",
    "text": "This section provides a description on how to use Forcefield.Components to calculate their energy.Amber.evaluate!"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Amber.load_from_json",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Amber.load_from_json",
    "category": "function",
    "text": "load_from_json(i_file::String)::Forcefield.Topology\n\nGather all topology components and return a Amber.Topology object, parsing a JSON file.\n\nExamples\n\njulia> Forcefield.load_from_json(json_file)\nForcefield.Topology(\n atoms=ProtoSyn.Forcefield.Atom[Forcefield.Atom(name=\"N\", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5]), ...],\n bonds=ProtoSyn.Forcefield.HarmonicBond[Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19), ...],\n angles=ProtoSyn.Forcefield.HarmonicAngle[Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92), ...],\n dihedralsCos=ProtoSyn.Forcefield.DihedralCos[Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0), ...])\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#Loaders-1",
    "page": "Forcefield",
    "title": "Loaders",
    "category": "section",
    "text": "This section provides a description on how to load Forcefield.Components from external sources.Amber.load_from_json"
},

{
    "location": "mutators.html#",
    "page": "Mutators",
    "title": "Mutators",
    "category": "page",
    "text": ""
},

{
    "location": "mutators.html#Mutators-1",
    "page": "Mutators",
    "title": "Mutators",
    "category": "section",
    "text": "Mutators are modules that introduce a single conformational change on the system, according to a set of runtime parameters. These can be combined into multiple steps, defining a Driver."
},

{
    "location": "mutators.html#ProtoSyn.Mutators.Dihedral.DihedralMutator",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Dihedral.DihedralMutator",
    "category": "type",
    "text": "DihedralMutator(dihedrals::Vector{Common.Dihedral}, angle_sampler::Function, p_mut::Float64, step_size::Float64)\n\nHolds all necessary parameters for the correct simulation of dihedral movements.\n\nArguments\n\ndihedrals::Vector{Common.Dihedral}: List of dihedrals avaliable to be rotated in dihedral movements.\nangle_sampler::Function: Function responsible for defining the new angle for the dihedral. Should return a Float64.\np_mut::Float64: Probability of rotation of each dihedral.\nstep_size::Float64: Scalar that defines the amount of change resulting for a dihedral movement.\n\nExamples\n\njulia> Mutators.Diehdral.DihedralMutator(dihedrals, 0.05, randn, 0.25)\nDihedralMutator(dihedrals=68, p_pmut=0.05, angle_sampler=randn, step_size=0.25)\n\nSee also: run!\n\n\n\n\n\n"
},

{
    "location": "mutators.html#ProtoSyn.Mutators.Dihedral.run!",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Dihedral.run!",
    "category": "function",
    "text": "run!(state::Common.State, mutator::DihedralMutator)\n\nIterate over a list of Common.Dihedral (dihedrals) and perform dihedral movements on the current Common.State. The probability of each dihedral undergo movements is defined in the DihedralMutator.p_mut. The new angle is obtained from DihedralMutator.angle_sampler, who should return a Float64 in radians. After movement, the Common.State is updated with the new conformation.\n\nExamples\n\njulia> Mutators.Dihedral.run!(state, mutator)\n\nSee also: Common.rotate_dihedral!\n\n\n\n\n\n"
},

{
    "location": "mutators.html#Dihedral-1",
    "page": "Mutators",
    "title": "Dihedral",
    "category": "section",
    "text": "CurrentModule = Mutators.DihedralThis section provides a description on the Dihedral mutator, responsible for performing a single dihedral movement according to a DihedralMutator set of parameters.DihedralMutator\nrun!"
},

{
    "location": "mutators.html#ProtoSyn.Mutators.Crankshaft.CrankshaftMutator",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Crankshaft.CrankshaftMutator",
    "category": "type",
    "text": "CrankshaftMutator(dihedrals::Vector{Common.Dihedral}, angle_sampler::Function, p_mut::Float64, step_size::Float64)\n\nHolds all necessary parameters for the correct simulation of crankshaft movements.\n\nArguments\n\ndihedrals::Vector{Common.Dihedral}: List of dihedrals avaliable to be rotated in crankshaft movements.\nangle_sampler::Function: Function responsible for defining the rotation angle. Should return a Float64.\np_mut::Float64: Probability of rotation of each pair of alpha carbons.\nstep_size::Float64: Scalar that defines the amount of change resulting for a crankshaft movement.\n\nExamples\n\njulia> Mutators.Crankshaft.CrankshaftMutator(dihedrals, randn, 0.05, 0.25)\nCrankshaftMutator(dihedrals=68, angle_sampler=randn, p_pmut=0.05, step_size=0.25)\n\njulia> Mutators.Crankshaft.CrankshaftMutator(dihedrals, randn)\nCrankshaftMutator(dihedrals=68, angle_sampler=randn, p_pmut=0.0, step_size=0.0)\n\nSee also: run!\n\n\n\n\n\n"
},

{
    "location": "mutators.html#ProtoSyn.Mutators.Crankshaft.run!",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Crankshaft.run!",
    "category": "function",
    "text": "run!(state::Common.State, mutator::CrankshaftMutator)\n\nIterate over a list of Common.Dihedral (dihedrals) and perform crankshaft movements on the current Common.State. The probability of each pair of alpha carbons undergo movement is defined in the CrankshaftMutator.p_mut. The new angle is obtained from CrankshaftMutator.angle_sampler, who should return a Float64 in radians. After movement, the Common.State is updated with the new conformation.\n\nExamples\n\njulia> Mutators.Crankshaft.run!(state, mutator)\n\nSee also: rotate_crankshaft!\n\n\n\n\n\n"
},

{
    "location": "mutators.html#ProtoSyn.Mutators.Crankshaft.rotate_crankshaft!",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Crankshaft.rotate_crankshaft!",
    "category": "function",
    "text": "rotate_crankshaft!(xyz::Array{Float64, 2}, dihedral1::Common.Dihedral, dihedral2::Common.Dihedral, angle::Float64)\n\nPerform a crankshaft movement, adding the provided angle (in radians) to all atoms between the alpha carbon in dihedral1 and dihedral2. The last residue side_chain is also rotated and should be relaxed back to a equilibrium state (See: SteepestDescent).\n\nExamples\n\njulia> Mutators.Dihedral.rotate_crankshaft!(state.xyz, dihedral1, dihedral2, π/2)\n\n\n\n\n\n"
},

{
    "location": "mutators.html#Crankshaft-1",
    "page": "Mutators",
    "title": "Crankshaft",
    "category": "section",
    "text": "CurrentModule = Mutators.CrankshaftA Crankshaft movement is a rotation of all atoms between two randomly chosen alpha carbons according to a set of CrankshaftMutator parameters.CrankshaftMutator\nrun!\nrotate_crankshaft!"
},

{
    "location": "drivers.html#",
    "page": "Drivers",
    "title": "Drivers",
    "category": "page",
    "text": ""
},

{
    "location": "drivers.html#Drivers-1",
    "page": "Drivers",
    "title": "Drivers",
    "category": "section",
    "text": "Drivers are a set set of functions that drive the simulation to new states, often exploring the conformational space of the system in Monte Carlo algorithms (with the combination of Mutators with acceptance/rejection rules) or the application of minimizers (such as the Steppest Descent Algorithm)."
},

{
    "location": "drivers.html#ProtoSyn.Drivers.MonteCarlo.MonteCarloDriver",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.MonteCarlo.MonteCarloDriver",
    "category": "type",
    "text": "MonteCarloDriver(sampler!::Function, evaluator!::Function, [, temperature::Float64 = 1.0, n_steps::Int64 = 0])\n\nDefine the runtime parameters for the Monte Carlo simulation. No sampler! movement is performed by default, since n_steps = 0.\n\nArguments\n\nsampler!::Function: Responsible for generating a new structure to be evaluated. This function should have the following signature:\n\nsampler!(state::Common.State)\n\nevaluator!::Function: Responsible for evaluating the system energy. This function should have the following signature:\n\nevaluator!(state::Common.State, do_forces::Bool)\n\ntemperature::Float64: (Optional) Temperature of the system, determines acceptance in the Metropolis algorithm (Default: 1.0).\nn_steps: (Optional) Total amount of steps to be performed (Default: 0).\n\nExamples\n\njulia> Drivers.MonteCarlo.MonteCarloDriver(my_sampler!, my_evaluator!, 10.0, 1000)\nMonteCarloDriver(sampler=my_sampler!, evaluator=my_evaluator!, temperature=10.0, n_steps=1000)\n\njulia> Drivers.MonteCarlo.MonteCarloDriver(my_sampler!, my_evaluator!)\nMonteCarloDriver(sampler=my_sampler!, evaluator=my_evaluator!, temperature=1.0, n_steps=0)\n\ntip: Tip\nBoth my_sampler! and my_evaluator! functions often contain pre-defined function avaliable in Mutators and Forcefield modules, respectively.\n\nSee also: run!\n\n\n\n\n\n"
},

{
    "location": "drivers.html#ProtoSyn.Drivers.MonteCarlo.run!",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.MonteCarlo.run!",
    "category": "function",
    "text": "run!(state::Common.State, driver::MonteCarloDriver[, callbacks::Tuple{Common.CallbackObject}...])\n\nRun the main body of the driver. Creates a new conformation based on driver.sampler!, evaluates the new conformation energy using driver.evaluator!, accepting it or not depending on the driver.temperature in a Metropolis algorithm. This Monte Carlo process is repeated for driver.n_steps, saving the accepted structures to state and calling all the callbacks. \n\nExamples\n\njulia> Drivers.MonteCarlo.run!(state, driver, my_callback1, my_callback2, my_callback3)\n\n\n\n\n\n"
},

{
    "location": "drivers.html#Monte-Carlo-1",
    "page": "Drivers",
    "title": "Monte Carlo",
    "category": "section",
    "text": "CurrentModule = Drivers.MonteCarloThis section provides a description on the Monte Carlo Driver. This Driver iterates over a set amount of n_steps (defined in the MonteCarloDriver), sampling new conformations to the Common.State and accepting or rejecting them based on the Metropolis Algorithm.MonteCarloDriver\nrun!"
},

{
    "location": "drivers.html#ProtoSyn.Drivers.SteepestDescent.SteepestDescentDriver",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.SteepestDescent.SteepestDescentDriver",
    "category": "type",
    "text": "SteepestDescentDriver(evaluator!::Function[, n_steps::Int64 = 0, f_tol::Float64 = 1e-3, max_step:Float64 = 0.1])\n\nDefine the runtime parameters for the Steepest Descent simulation. If n_steps is zero, a single point energy calculation is performed.\n\nArguments\n\nevaluator!::Function: Responsible for evaluating the current state.energy and calculate the resulting forces. This function should have the following signature:\n\nevaluator!(state::Common.State, do_forces::Bool)\n\nn_steps: (Optional) Total amount of steps to be performed (if convergence is not achieved before) (Default: 0).\nf_tol: (Optional) Force tolerance. Defines a finalization criteria, as the steepest descent is considered converged if the maximum force calculated is below this value (Default = 1e-3).\nmax_step: (Optional) Defines the maximum value ɣ that the system can jump when applying the forces (Default: 0.1).\nostream: (Optional) Defines the output stream for logging\n\nExamples\n\njulia> Drivers.SteepestDescent.SteepestDescentDriver(my_evaluator!, 100, 1e-3, 0.1)\nSteepestDescentDriver(evaluator=my_evaluator!, n_steps=100, f_tol=1e-3, max_step=0.1)\n\njulia> Drivers.SteepestDescent.SteepestDescentDriver(my_evaluator!, f_tol = 1e-6)\nSteepestDescentDriver(evaluator=my_evaluator!, n_steps=0, f_tol=1e-6, max_step=0.1)\n\ntip: Tip\nThe my_evaluator! function often contains an aggregation of pre-defined functions avaliable in Forcefield. It is possible to combine such functions using the @faggregator macro.\n\nSee also: Amber.evaluate! run!\n\n\n\n\n\n"
},

{
    "location": "drivers.html#ProtoSyn.Drivers.SteepestDescent.run!",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.SteepestDescent.run!",
    "category": "function",
    "text": "run!(state::Common.State, driver::SteepestDescentDriver[, callback::Union{Common.CallbackObject, Nothing} = nothing])\n\nRun the main body of the Driver. If driver.n_steps is zero, a single point energy calculation is performed.\n\nArguments\n\nstate::Common.State: Current state of the system to be modified.\ndriver::SteepestDescentDriver: Defines the parameters for the SteepestDescent simulation. See SteepestDescentDriver.\ncallbacks::Vararg{Common.CallbackObject, N}: (Optional) Tuple of CallbackObjects (Default: empty).\n\ntip: Tip\nThe callback function often contains a Print function.\n\nExamples\n\njulia> Drivers.SteepestDescent.run(state, steepest_descent_driver, callback1, callback2, callback3)\n\n\n\n\n\n"
},

{
    "location": "drivers.html#Steepest-Descent-1",
    "page": "Drivers",
    "title": "Steepest Descent",
    "category": "section",
    "text": "CurrentModule = Drivers.SteepestDescentThis section provides a description on the Steepest Descent Driver. This Driver attempts to minimize the system energy based on the provided Topology, as it calculates the forces acting on each atom according to the defined SteepestDescentDriver. The finalization criteria is:Maximum number of steps was achieved (SteepestDescentDriver.n_steps).\nMaximum force calculated is below the force tolerance (SteepestDescentDriver.f_tol).\nThe gamma (γ) applied in the next step of the Steepest Descent is below machine precision.SteepestDescentDriver\nrun!"
},

{
    "location": "print.html#",
    "page": "Print",
    "title": "Print",
    "category": "page",
    "text": ""
},

{
    "location": "print.html#ProtoSyn.Print.as_xyz",
    "page": "Print",
    "title": "ProtoSyn.Print.as_xyz",
    "category": "function",
    "text": "as_xyz(io:IO, state::Common.State[, title::String = \"mol\"])\n\nPrint the current Common.State as a .xyz file to the output io.\n\nExamples\n\njulia> Drivers.MonteCarlo.load_parameters(file_xyz, state, title = \"molecule\")\n\n\n\n\n\nas_xyz(state::Common.State[, title::String = \"mol\"])::String\n\nPrint the current Common.State in .xyz format and returns a String.\n\nExamples\n\njulia> Drivers.MonteCarlo.load_parameters(state, title = \"molecule\")\n2\n molecule\n N      -0.0040    0.2990    0.0000\n H1      0.1200    1.3010    0.0000\n\n\n\n\n\n"
},

{
    "location": "print.html#Print-1",
    "page": "Print",
    "title": "Print",
    "category": "section",
    "text": "Print functions serve as output of the system, writting information in standardized format.CurrentModule = Printas_xyz"
},

{
    "location": "aux.html#",
    "page": "Aux",
    "title": "Aux",
    "category": "page",
    "text": ""
},

{
    "location": "aux.html#ProtoSyn.Aux.read_JSON",
    "page": "Aux",
    "title": "ProtoSyn.Aux.read_JSON",
    "category": "function",
    "text": "read_JSON(i_file::String)\n\nRead a JSON file and return a dictionary with the content.\n\nExamples\n\njulia> Aux.read_JSON(\"i_file.json\")\nDict{Any, Any}()\n\n\n\n\n\n"
},

{
    "location": "aux.html#ProtoSyn.Aux.rotation_matrix_from_axis_angle",
    "page": "Aux",
    "title": "ProtoSyn.Aux.rotation_matrix_from_axis_angle",
    "category": "function",
    "text": "rotation_matrix_from_axis_angle(axis::Vector{Float64}, angle::Float64)\n\nReturn a rotation matrix based on the provided axis and angle (in radians).\n\nExamples\n\njulia> Aux.rotation_matrix_from_axis_angle([1.1, 2.2, 3.3], π/2)\n3×3 Array{Float64,2}:\n  0.0714286  -0.658927  0.748808\n  0.944641    0.285714  0.16131 \n -0.320237    0.695833  0.642857\n\nSee also: rotate_dihedral!\n\n\n\n\n\n"
},

{
    "location": "aux.html#ProtoSyn.Aux.calc_dih_angle",
    "page": "Aux",
    "title": "ProtoSyn.Aux.calc_dih_angle",
    "category": "function",
    "text": "calc_dih_angle(a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, a4::Vector{Float64})\n\nCalculates the dihedral angle produced between a1, a2, a3 and a4, in radians.\n\nExamples\n\njulia> Aux.calc_dih_angle([1.0, 1.0, 1.0], [2.1, 2.1, 2.1], [3.0, 2.0, 5.0], [5.0, 5.0, 5.0])\n3.141592653589793\n\nSee also: apply_initial_conf!\n\n\n\n\n\n"
},

{
    "location": "aux.html#Aux-1",
    "page": "Aux",
    "title": "Aux",
    "category": "section",
    "text": "This section provides a description of miscellaneous auxiliary functions.CurrentModule = Auxread_JSON\nrotation_matrix_from_axis_angle\ncalc_dih_angle"
},

{
    "location": "json.html#",
    "page": "Input JSON Schemas",
    "title": "Input JSON Schemas",
    "category": "page",
    "text": ""
},

{
    "location": "json.html#Input-JSON-Schema-1",
    "page": "Input JSON Schemas",
    "title": "Input JSON Schema",
    "category": "section",
    "text": "This section describes in detail the general schematics of the input JSON so it can be easily read by ProtoSyn. "
},

{
    "location": "json.html#Topology-1",
    "page": "Input JSON Schemas",
    "title": "Topology",
    "category": "section",
    "text": "Protosyn receives various types of topology, depending on the Driver chosen by the user."
},

{
    "location": "json.html#Forcefield-Topology-1",
    "page": "Input JSON Schemas",
    "title": "Forcefield Topology",
    "category": "section",
    "text": "Forcedield topology holds the information regarding the system bonds, angles, dihedrals and non-bonded interactions. ProtoSyn offers a tool (ProtoSyn.jl/tools/tpr2json.py) that easily reads .tpr file formats from GROMACS (dumped to human readable format using gmx dump) and converts the obtained information to a correct JSON format. <details>\n<summary style=\"background-color:#1abc9c;border-radius: 13px;height: 30px;display: flex;align-items: center\">&ensp;Detailed forcefield topology JSON description▼</summary>\n<p>{\n    \"dihedraltypes\": {\n        ID::String: {\n            \"phi\"    : Float64,    // Equilibrium angle value (degrees)\n            \"cp\"     : Float64,    // Force constant (kJ mol⁻¹)\n            \"mult\"   : Float64     // Multiplicity\n        },\n        ...\n    },\n    \"angletypes\": {\n        ID::String: {\n            \"th\"     : Float64,    // Equilibrium angle value (degrees)\n            \"ct\"     : Float64     // Force constant (kJ mol⁻¹ rad⁻²)\n        },\n        ...\n    },\n    \"bondtypes\": {\n        ID::String: {\n            \"cb\"     : Float64,    // Force constant (kJ mol⁻¹ nm⁻²)\n            \"b0\"     : Float64     // Equilibrium bond length (nm)\n        },\n        ...\n    },\n    \"atomtypes\": {\n        ID::String: {\n            \"sigma\"  : Float64,    // Finite distance at which the inter-particle potential is zero (nm)\n            \"epsilon\": Float64,    // Depth of the potential well (kJ mol⁻¹)\n            \"name\"   : String      // Name of the atomtype\n        },\n        ...\n    },\n    \"dihedrals\": [\n        {\n            \"a1\"     : Int64,      // Global atom 1 index\n            \"a2\"     : Int64,      // Global atom 2 index\n            \"a3\"     : Int64,      // Global atom 3 index\n            \"a4\"     : Int64,      // Global atom 4 index\n            \"type\"   : Int64       // ID of dihedraltypes\n        },\n        ...\n    ],\n    \"angles\": [\n        {\n            \"a1\"     : Int64,      // Global atom 1 index\n            \"a2\"     : Int64,      // Global atom 2 index\n            \"a3\"     : Int64,      // Global atom 3 index\n            \"type\"   : Int64       // ID of angletypes\n        },\n        ...\n    ],\n    \"bonds\": [\n        {\n            \"a1\"     : Int64,      // Global atom 1 index\n            \"a2\"     : Int64,      // Global atom 2 index\n            \"type\"   : Int64       // ID of bondtypes\n        },\n        ...\n    ]\n    \"atoms\": [\n        {\n            \"q\"      : Float64,    // Atom charge (eletron)\n            \"type\"   : Int64,      // ID of atomtypes\n        },\n        ...\n    ],\n    \"excls\": {\n        a1::Int64: Array{Int64, 1}, // Global indices of the atoms excluded from non-bonded interactions with a1\n        ...\n    },\n    \"pairs\": [\n        {\n            \"a1\"     : Int64,       // Global index of a1\n            \"a2\"     : Int64,       // Global index of a2 (At 3 connections away from a1)\n        },\n        ...\n    ]\n}</p>\n</details>"
},

{
    "location": "json.html#Monte-Carlo-Topology-1",
    "page": "Input JSON Schemas",
    "title": "Monte Carlo Topology",
    "category": "section",
    "text": "Monte Carlo topology holds the information regarding the rotatable dihedrals and the residues present in the system. This allows for correct rotation in Monte Carlo moves.<details>\n<summary style=\"background-color:#1abc9c;border-radius: 13px;height: 30px;display: flex;align-items: center\">&ensp;Detailed Monte Carlo JSON description▼</summary>\n<p>{\n    \"dihedrals\": [\n        {\n            \"parent\" : Int64,           // Index of parent residue in list of residues\n            \"movable\": Array{Int64, 1}, // Global indices of atoms moved in the parent residue when rotating this dihedral\n            \"a1\"     : Int64,           // Global atom 1 index\n            \"a2\"     : Int64,           // Global atom 2 index\n            \"a3\"     : Int64,           // Global atom 3 index\n            \"a4\"     : Int64,           // Global atom 4 index\n            \"type\"   : String           // Dihedral type (i.e. \"PHI\", \"PSI\", ...)\n        },\n        ...\n    ],\n    \"residues\": [\n        {\n            \"next\"   : Int64,           // Index of the next residue in list of residues\n            \"atoms\"  : Array{Int64, 1}, // List of global indices of the atoms included in this residue\n            \"type\"   : String,          // Residue type (i.e. \"V\", \"K\", ...)\n            \"n\"      : Int64            // Residue identifier\n        },\n        ...\n    ]\n}</p>\n</details>"
},

]}
