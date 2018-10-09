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
    "text": "PAI NATAL MORRE ENGASGADO A COMER COELHO DA PÁSCOA À CAÇADOR!!(Image: Santa)function fancyAlert(arg) {\n  if(arg) {\n    $.facebox({div:\'#foo\'})\n  }\n}"
},

{
    "location": "guide.html#",
    "page": "Guide",
    "title": "Guide",
    "category": "page",
    "text": ""
},

{
    "location": "guide.html#Package-guide-1",
    "page": "Guide",
    "title": "Package guide",
    "category": "section",
    "text": ""
},

{
    "location": "guide.html#Installation-1",
    "page": "Guide",
    "title": "Installation",
    "category": "section",
    "text": "Protosyn is still not a registered package. Please download or clone the source code from GitHub."
},

{
    "location": "guide.html#Usage-1",
    "page": "Guide",
    "title": "Usage",
    "category": "section",
    "text": "What follows is a step-by-step guide on how to use the ProtoSyn library."
},

{
    "location": "guide.html#General-overview-and-workflow-1",
    "page": "Guide",
    "title": "General overview and workflow",
    "category": "section",
    "text": "ProtoSyn is a stuctural sampling library designed to explore the conformational space of molecules. Specifically, ProtoSyn was developed to study proteins and how they fold in 3D space. With ProtoSyn, the user is able to easily integrate different modules in order to perform distinct actions that change and evaluate the system state.The ProtoSyn Flow is as follows:Load the system initial Common.State from a structural file, such as a .pdb or a .gro;\nLoad all the necessary topologies, describing the system bonded and non-bonded interactions;\nChoose a Driver. This will change the system conformation step after step;\nFor certain Drivers, a sampler is required. This function is responsible for mutating the system in a certain way. Several Mutators are available.\nDefine the evaluating function. Regardless of the Driver chosen, this function will evaluate the system fitness and determine the next step.\nDefine the callback function. At this point the user\'s program requires a way to comunicate with the outside. Callback functions are employed to output the produced information by the Driver.\nRun the program and do science!"
},

{
    "location": "guide.html#Examples-1",
    "page": "Guide",
    "title": "Examples",
    "category": "section",
    "text": "As an example, a Steepest Descent Algorithm written in Julia will be explained in detail. As explained in the ProtoSyn Flow we have the following steps:Load the system initial Common.State.julia> state = Common.load_from_pdb(\"protein.pdb\")\nCommon.State(size=56, energy=Null, xyz=[1.1 2.1 1.2; 2.2 1.2 5.2, ...], forces=[0.0 0.0 0.0; 0.0 0.0 0.0, ...], atnames=[\"C\", \"O\", ...])\njulia> state.energy = Forcefield.Energy()Load all the necessary topologies.\nThe topology can be loaded by several ways. ProtoSyn by default expects a JSON file depicting the necessary information. For the Steepest Descent Algorithm, bonds, angles, dihedrals and non-bonded parameters should be depicted in the input JSON file.julia> topology = Forcefield.load_from_json(\"topology.json\")\nForcefield.Topology(\natoms=ProtoSyn.Forcefield.Atom[Forcefield.Atom(name=\"N\", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5]), ...],\nbonds=ProtoSyn.Forcefield.HarmonicBond[Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19), ...],\nangles=ProtoSyn.Forcefield.HarmonicAngle[Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92), ...],\ndihedralsCos=ProtoSyn.Forcefield.DihedralCos[Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0), ...])Choose a Driver.\nThis step includes loading the necessary runtime parameters. These are specific to the chosen Driver.julia> params = Drivers.SteepestDescent.ConfigParameters(10000, 100, 1e-3, 0.1)\nDrivers.SteepestDescent.ConfigParameters(n_steps=10000, log_freq=100, f_tol=1e-3, max_step=0.1)Define the sampler.\nAs stated in the documentation, the Steepest Descent Driver does not require a sampler, as the Driver itself is responsible for choosing the next step without the aid of a Mutator.\nDefine the evaluating function.\nWhen defining the necessary functions in ProtoSyn, careful care needs to be taken to match the expected function signature, described in the documentation.julia> function my_evaluator!(st::Common.State, do_forces::Bool)\n        return Forcefield.evalenergy!(topology, st, cut_off=1.2, do_forces=do_forces)\n    end\nmy_evaluator! (generic function with 1 method)tip: Tip\nEven though the my_evaluator! function does not directly receive the topology as an argument, it is still able to access it as it was defined in the main body of our program.Define the callback function.\nAltough optional, this function allows the user to easily retrieve the information being produced by the Driver. ProtoSyn includes some functions for this, such as to Print the current state to a structural file.julia> output = open(\"output.xyz\", \"w\")\nIOStream(<output.xyz>)\njulia> function my_callback(st::Common.State, step::Int)\n            Print.as_xyz(st, ostream = output, title = \"Step $step\")\n        end\nmy_callback (generic function with 1 method)Run the program and do science!\nAll the necessary variables and functions have been defined in the main body of our program and it is ready to be deployed.julia> Drivers.SteepestDescent.run!(state, my_evaluator!, params, callback = my_callback)For more detailed information, please reference to the Manual."
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
    "location": "common.html#ProtoSyn.Common.NullEnergy",
    "page": "Common",
    "title": "ProtoSyn.Common.NullEnergy",
    "category": "type",
    "text": "NullEnergy()\n\nEmpty placeholder energy container.\n\nExamples\n\njulia> Common.NullEnergy()\nNull\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.State",
    "page": "Common",
    "title": "ProtoSyn.Common.State",
    "category": "type",
    "text": "State(size::Int64, energy::AbstractEnergy, xyz::Array{Float64, 2}, forces::Array{Float64, 2}, atnames::Array{String, 1})\n\nDefine the current state of the system, containing the atoms positions, energy and forces applied. If only size::Int64 is provided, an empty State with the given size is created with zeros.\n\nArguments\n\nsize::Int64: Atom count in system.\nenergy::AbstractEnergy: Current energy of the system (kJ mol⁻¹).\nxyz::Array{Float64, 2}: Atom positions in 3 dimensions.\nforces::Array{Float64, 2}: Forces applied in each dimension to each atom (kJ mol⁻¹ nm⁻¹)\natnames::Array{String, 1}: List of atom names.\n\nExamples\n\njulia> Common.State(3)\nCommon.State(size=3, energy=Null, xyz=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], forces=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], atnames=String[])\n\njulia> Common.State(2, Common.NullEnergy(), [1.1 1.1 1.1; 2.2 2.2 2.2], zeros(2, 3), [\"C\", \"O\"])\nCommon.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], atnames=[\"C\", \"O\"])\n\n\n\n\n\n"
},

{
    "location": "common.html#ProtoSyn.Common.Residue",
    "page": "Common",
    "title": "ProtoSyn.Common.Residue",
    "category": "type",
    "text": "Residue(atoms::Array{Int64, 1}, next::Union{Residue, Int64, Nothing}, name::String)\n\nDefine a residue as part of the system. \n\nArguments\n\natoms::Array{Int64, 1}: list of atom global indices in this residue.\nnext::Union{Residue, Int64, Nothing}: Next residue in the system, attached to this. Is preferably a Residue instance, but can in certain cases be the index in a Residue list or empty (Nothing).\nname::String: Name of the residue. If the residue describes an aminoacid, the correspondent letter is suggested.\n\nExamples\n\njulia> Common.Residue([1, 2, 3, 4], Common.Residue([5, 6, 7, 8], nothing, \"V\"), \"E\")\nCommon.Residue(atoms=[1, 2, 3, 4], next=V, name=E)\n\n\n\n\n\n"
},

{
    "location": "common.html#Components-1",
    "page": "Common",
    "title": "Components",
    "category": "section",
    "text": "This section provides a description of the Common components, such as Common.State.NullEnergy\nState\nResidue"
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
    "text": "load_from_pdb(i_file::String)::Common.State\n\nReturn a new Common.State by loading the atom positions and names from the input .pdb file. As a default, state.energy is NullEnergy and state.forces are set to zero.\n\nExamples\n\njulia> Common.load_from_pdb(\"molecule.pdb\")\nCommon.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0], atnames=[\"C\", \"O\"])\n\nSee also: load_from_gro\n\n\n\n\n\n"
},

{
    "location": "common.html#Loader-1",
    "page": "Common",
    "title": "Loader",
    "category": "section",
    "text": "This section provides a description on how to load a new Common.State.load_from_gro\nload_from_pdb"
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
    "text": "CurrentModule = Forcefield"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.HarmonicBond",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.HarmonicBond",
    "category": "type",
    "text": "HarmonicBond(a1::Int64, a2::Int64, k::Float64, b0::Float64)\n\nHarmonic Bond of the form\n\nE(r_ab) = frac12k_ab(r_ab - b_0)^2\n\nwhere\n\nr_ab = vecr_ab = vecr_b - vecr_a\n\nArguments\n\na1::Int64, a2::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹ nm⁻²).\nb0::Float64: equilibrium bond length (nm).\n\nExamples\n\njulia> Forcefield.HarmonicBond(1, 2, 2500, 0.19)\nForcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19)\n\nSee algo: evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.HarmonicAngle",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.HarmonicAngle",
    "category": "type",
    "text": "HarmonicAngle(a1::Int64, a2::Int64, a3::Int64, k::Float64, θ::Float64)\n\nHarmonic Angle of the form\n\nE(θ_abc)=frac12k_abc(theta_abc-theta)^2\n\nArguments\n\na1::Int64, a2::Int64, a3::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹ rad⁻²).\nθ::Float64: equilibrium angle value (rad).\n\nExamples\n\njulia> Forcefield.HarmonicAngle(1, 2, 3, 670.0, 1.92)\nForcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92)\n\nSee algo: evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.DihedralCos",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.DihedralCos",
    "category": "type",
    "text": "DihedralCos(a1::Int64, a2::Int64, a3::Int64, a4::Int64, k::Float64, θ::Float64, mult::Float64)\n\nPeriodic Dihedral of the form\n\nE(phi_abcd)=K_phi(1+cos(nphi_abcd-phi))\n\nArguments\n\na1::Int64, a2::Int64, a3::Int64, a4::Int64: global atom indices.\nk::Float64: force constant (kJ mol⁻¹).\nθ::Float64: equilibrium angle value (rad).\nmult::Float64: multiplicity.\n\nExamples\n\njulia> Forcefield.DihedralCos(1, 2, 3, 4, 10.46, 180.0, 2.0)\nForcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0)\n\nSee algo: evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Atom",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Atom",
    "category": "type",
    "text": "Atom(name::String, σ::Float64, ϵ::Float64, q::Float64, excls::Array{Int64, 1}, pairs::Array{Int64, 1})\n\nDefine an atom.  σ, ϵ and q describe the non-bonded interactions between atoms:\n\nThe Lennard-Jones interaction is in the form:\n\nE(r_ab) = 4ϵ_ableft ( (fracσ_abr_ab)^12-(fracσ_abr_ab)^6right )\n\nwhere the Lorentz-Berthelot rule is applied. σ is the arithmetic average and ϵ is the geometric average:\n\nσ_ab=fracσ_a+σ_b2\n\nϵ_ab=sqrt(ϵ_aϵ_b)\n\nFor this reason, σ and ϵ are applied here in the reduced form: fracσ2 and sqrtϵ.\n\nThe Coulomb interation is in the form:\n\nE(r_ab)=k_ϵfracq_aq_br_ab^2\n\nwhere\n\nk_ϵ=frac14πϵ_0=138935485kJnmmol¹e¹\n\nFor this reason, q is applied here in the reduced form: qtimes sqrtk_ϵ\n\nExclusion list contains all atom indices who are excluded from non-bonded interactions (i.e. are at 3 or less connections from this atom - includes pairs). Pair list contains atoms that are at 3 connections from this atom, and are involved in 1-4 interactions (and have a different combination rule as a result).\n\nArguments\n\nname::String: Atom name (example: \"C\", \"H\", etc).\nσ::Float64: finite distance at which the inter-particle potential is zero (nm).\nϵ::Float64: depth of the potential well (kJ mol⁻¹).\nq::Float64: atom charge (eletron).\nexcls::Array{Int64, 1}: exclusion list (as global atom indices).\npairs::Array{Int64, 1}: pair list containing atoms that interfere in 1-4 interations (as global atom indices)\n\nExamples\n\njulia> Forcefield.Atom(\"N\", 0.325, 0.711, 0.0017, [0, 1, 2, 3, 4, 5], [4, 5])\nForcefield.Atom(name=\"N\", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5])\n\nSee algo: evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Topology",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Topology",
    "category": "type",
    "text": "Topology(atoms::Array{Atom}, bonds::Array{HarmonicBond}, angles::Array{HarmonicAngle}, dihedralsCos::Array{DihedralCos})\n\nGather all topology components.\n\nArguments\n\natoms::Array{Atoms}\nbonds::Array{HarmonicBond}\nangles::Array{HarmonicAngle}\ndihedralsCos::Array{DihedralCos}\n\nExamples\n\njulia> Forcefield.Forcefield(atoms, bonds, angles, dihedrals)\nForcefield.Topology(\n atoms=ProtoSyn.Forcefield.Atom[Forcefield.Atom(name=\"N\", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5]), ...],\n bonds=ProtoSyn.Forcefield.HarmonicBond[Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19), ...],\n angles=ProtoSyn.Forcefield.HarmonicAngle[Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92), ...],\n dihedralsCos=ProtoSyn.Forcefield.DihedralCos[Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0), ...])\n\nSee also: Forcefield.load_from_json\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.Energy",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.Energy",
    "category": "type",
    "text": "Energy(eBond::Float64, eAngle::Float64, eDihedral::Float64, eLJ::Float64, eLJ14::Float64, eCoulomb::Float64, eCoulomb14::Float64, eTotal::Float64)\n\nEnergy components.\n\nExamples\n\njulia> Forcefield.Energy()\nForcefield.Energy(eBond=0.0, eAngle=0.0, eDihedral=0.0, eLJ=0.0, eLJ14=0.0, eCoulomb=0.0, eCoulomb14=0.0, eTotal=0.0)\n\njulia> Forcefield.Energy(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 2.8)\nForcefield.Energy(eBond=0.1, eAngle=0.2, eDihedral=0.3, eLJ=0.4, eLJ14=0.5, eCoulomb=0.6, eCoulomb14=0.7, eTotal=2.8)\n\nSee also: Forcefield.evalenergy!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#Components-1",
    "page": "Forcefield",
    "title": "Components",
    "category": "section",
    "text": "This section provides a description on all the individual components that comprise the Forcefield.HarmonicBond\nHarmonicAngle\nDihedralCos\nAtom\nTopology\nEnergy"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.evaluate!",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.evaluate!",
    "category": "function",
    "text": "evaluate!(bonds::Array{Forcefield.HarmonicBond}, state::Common.State[, do_forces::Bool = false])::Float64\n\n\n\n\n\nevaluate!(angles::Array{Forcefield.HarmonicAngle}, state::Common.State, do_forces::Bool = false)::Float64\n\n\n\n\n\nevaluate!(dihedralsCos::Array{Forcefield.DihedralCos}, state::Common.State, do_forces::Bool = false)::Float64\n\n\n\n\n\nevaluate!(atoms::Array{Forcefield.Atom}, state::Common.State, do_forces::Bool = false)::Float64\n\nEvaluate an array of Forcefield.Components using the current Common.State, calculate and update state.energy according to the equations defined in each component. If do_forces flag is set to true, calculate and update state.forces. Non-bonded interactions are only assessed if the distance between atoms is below the defined cut_off value. Return the component energy value (kJ mol⁻¹).\n\nExamples\n\njulia> Forcefield.evaluate!(bonds, state)\n0.500\n\nSee also: evalenergy! Forcefield.HarmonicBond Forcefield.HarmonicAngle Forcefield.DihedralCos Forcefield.Atom\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.evalenergy!",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.evalenergy!",
    "category": "function",
    "text": "evalenergy!(topology::Forcefield.Topology, state::Common.State[, cut_off::Float64 = 2.0, do_forces::Bool = false])::Float64\n\nEvaluate the current Common.State energy according to the defined Forcefield.Topology. If do_forces bool is set to true, calculate and update state.forces. Non-bonded interactions are only assessed if the distance between atoms is below the defined cut_off value. Return state.energy.eTotal value (kJ mol⁻¹).\n\nExamples\n\njulia> Forcefield.evalenergy!(topology, state, cut_off = Inf)\n0.500\n\nSee also: evaluate!\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#Evaluators-1",
    "page": "Forcefield",
    "title": "Evaluators",
    "category": "section",
    "text": "This section provides a description on how to use Forcefield.Components to calculate their energy.evaluate!\nevalenergy!"
},

{
    "location": "forcefield.html#ProtoSyn.Forcefield.load_from_json",
    "page": "Forcefield",
    "title": "ProtoSyn.Forcefield.load_from_json",
    "category": "function",
    "text": "load_from_json(i_file::String)::Forcefield.Topology\n\nGather all topology components and return a Forcefield.Topology object, parsing a JSON file.\n\nExamples\n\njulia> Forcefield.load_from_json(json_file)\nForcefield.Topology(\n atoms=ProtoSyn.Forcefield.Atom[Forcefield.Atom(name=\"N\", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5]), ...],\n bonds=ProtoSyn.Forcefield.HarmonicBond[Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19), ...],\n angles=ProtoSyn.Forcefield.HarmonicAngle[Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92), ...],\n dihedralsCos=ProtoSyn.Forcefield.DihedralCos[Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0), ...])\n\n\n\n\n\n"
},

{
    "location": "forcefield.html#Loaders-1",
    "page": "Forcefield",
    "title": "Loaders",
    "category": "section",
    "text": "This section provides a description on how to load Forcefield.Components from external sources.load_from_json"
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
    "location": "mutators.html#ProtoSyn.Mutators.Dihedral.ConfigParameters",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Dihedral.ConfigParameters",
    "category": "type",
    "text": "ConfigParameters(p_mut::Float64 = 0.1)\n\nDefine the runtime parameters for Dihedral movements.\n\nArguments\n\np_mut::Float64: Probability of mutation of this dihedral (Default: 0.1).\n\nExamples\n\njulia> Mutators.Diehdral.ConfigParameters(0.2)\nConfigParameters(p_mut=0.2)\n\njulia> Mutators.Diehdral.ConfigParameters()\nConfigParameters(p_mut=0.1)\n\n\n\n\n\n"
},

{
    "location": "mutators.html#ProtoSyn.Mutators.Dihedral.NewDihedral",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Dihedral.NewDihedral",
    "category": "type",
    "text": "NewDihedral(a1::Int64, a2::Int64, a3::Int64, a4::Int64, movable::Array{Int64, 1}, residue::Union{Common.Residue, Int64}, dtype::String)\n\nDefine a dihedral.\n\nArguments\n\na1::Int64, a2::Int64, a3::Int64, a4::Int64: global atom indices.\nmovable::Array{Int64, 1}: List of global atom indices that will be moved during the dihedral movement in this residue.\nresidue::Union{Common.Residue, Int64}: Residue that this dihedral belongs to. Should be a Common.Residue object.\ndtype::String: Dihedral type (i.e. \"PHI\", \"PSI\", ...)\n\nExamples\n\njulia> Mutators.Diehdral.NewDihedral(2, 3, 5, 7, [5, 6], Common.Residue([1, 2, 3, 4, 5, 6], (...), \"A\"), \"PSI\")\nDihedral(a1=2, a2=3, a3=5, a4=7, movable=[5, 6], residue=Common.Residue(atoms=[1, 2, 3, 4, 5, 6], next=V, name=A), type=\"PSI\")\n\nSee also: load_topology\n\n\n\n\n\n"
},

{
    "location": "mutators.html#ProtoSyn.Mutators.Dihedral.load_topology",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Dihedral.load_topology",
    "category": "function",
    "text": "load_topology(p::Dict{String, Any})\n\nParse a dictionary containing the dihedral and residue topology. Return a NewDihedral array and a Common.Residue array.\n\nExamples\n\njulia> Mutators.Diehdral.load_topology(p)\n(ProtoSyn.Mutators.Dihedral.NewDihedral[...], ProtoSyn.Common.Residue[...])\n\nSee also: Aux.read_JSON\n\n\n\n\n\n"
},

{
    "location": "mutators.html#ProtoSyn.Mutators.Dihedral.run!",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Dihedral.run!",
    "category": "function",
    "text": "run!(state::Common.State, dihedrals::Array{NewDihedral, 1}, params::Config.Parameters, angle_sampler::Function[, ostream::IO = stdout])\n\nIterate over a list of NewDihedral (dihedrals) and perform dihedral movements on the current Common.State. The probability of each dihedral undergo movements is defined in the ConfigParameters (params.pmut). The new angle is obtained from angle_sampler, who should return a Float64 in radians. After movement, the Common.State is updated with the new conformation. Any logging is written to ostream (Default: stdout).\n\nExamples\n\njulia> Mutators.Diehdral.run!(state, dihedrals, params, () -> randn())\n\nSee also: rotate_dihedral!\n\n\n\n\n\n"
},

{
    "location": "mutators.html#ProtoSyn.Mutators.Dihedral.rotate_dihedral!",
    "page": "Mutators",
    "title": "ProtoSyn.Mutators.Dihedral.rotate_dihedral!",
    "category": "function",
    "text": "rotate_dihedral!(xyz::Array{Float64, 2}, dihedral::NewDihedral, angle::Float64)\n\nPerform a dihedral movement, adding the provided angle (in radians). If the dihedral.dtype is \"PHI\" or \"PSI\" the dihedral.residue.next is also rotated and this is propagated recursively until the end of the molecule. \n\nExamples\n\njulia> Mutators.Diehdral.rotate_dihedral(state.xyz, dihedral, π/2)\n\nSee also: run! Aux.rotation_matrix_from_axis_angle\n\n\n\n\n\n"
},

{
    "location": "mutators.html#Dihedral-1",
    "page": "Mutators",
    "title": "Dihedral",
    "category": "section",
    "text": "CurrentModule = Mutators.DihedralThis section provides a description on the Dihedral mutator, responsible for performing a single dihedral movement according to a set of ConfigParameters.ConfigParameters\nNewDihedral\nload_topology\nrun!\nrotate_dihedral!"
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
    "text": "Drivers are a set set of functions that drive the simulation to new states, often exploring the conformational space of the system in Monte Carlo algorithms (with the combination of Mutators with acceptance/rejection rules) or the application of minimizers (such as the Steppest Descent Algorithm). For each Driver a set of required callback functions is further explained, detailing the necessary signature."
},

{
    "location": "drivers.html#ProtoSyn.Drivers.MonteCarlo.ConfigParameters",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.MonteCarlo.ConfigParameters",
    "category": "type",
    "text": "ConfigParameters(n_steps::Int64 = 1000, temperature::Float64 = 0.5)\n\nDefine the runtime parameters for the Monte Carlo Driver.\n\nArguments\n\nn_steps: Total amount of steps to be performed (Default: 1000).\ntemperature: Temperature value to be used in the Metropolis Algorithm (Default: 0.5).\n\nExamples\n\njulia> Drivers.MonteCarlo.ConfigParameters(1000, 0.5)\nDrivers.MonteCarlo.ConfigParameters(n_steps=1000, temperature=0.5)\n\njulia> Drivers.MonteCarlo.ConfigParameters(temperature = 1.2)\nDrivers.MonteCarlo.ConfigParameters(n_steps=1000, temperature=1.2)\n\nSee also: load_parameters\n\n\n\n\n\n"
},

{
    "location": "drivers.html#ProtoSyn.Drivers.MonteCarlo.load_parameters",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.MonteCarlo.load_parameters",
    "category": "function",
    "text": "load_parameters(p::Dict{String, Any})::ConfigParameters\n\nLoad the ConfigParameters from a dictionary.\n\nExamples\n\njulia> Drivers.MonteCarlo.load_parameters(p)\nDrivers.MonteCarlo.ConfigParameters(n_steps=1000, temperature=0.5)\n\nSee also: Aux.read_JSON\n\n\n\n\n\n"
},

{
    "location": "drivers.html#ProtoSyn.Drivers.MonteCarlo.run!",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.MonteCarlo.run!",
    "category": "function",
    "text": "run!(state::Common.State, sampler!::Function, evaluator!::Function, params::ConfigParameters[, ostream::IO = stdout, callback::Union{Function, Nothing} = nothing])\n\nRun the main body of the Driver.\n\nArguments\n\nstate::Common.State: Current state of the system to be modified.\nsampler!::Function: Responsible for mutating the current state. This function should have the following signature:\n\nsampler!(state::Common.State)\n\nevaluator!::Function: Responsible for evaluating the current state.energy for the Metropolis Algorithm. This function should have the following signature:\n\nevaluator!(state::Common.State, do_forces::Bool)\n\nparams::ConfigParameters: Hold the runtime parameters of the Monte Carlo Driver. \nostream::IO: (Optional) Any logging will be written to the supplied ostream.\ncallback::Union{Function, Nothing}: (Optional) If present, this function will be called if the new conformation was accepted. This function should have the following signature:\n\ncallback(state::Common.State, step::Int64)\n\ntip: Tip\nThe callback function is often a Print function.\n\nExamples\n\njulia> Drivers.MonteCarlo.run(state, my_sampler!, my_evaluator!, params!, callback = my_callback)\n\nSee also: Mutators Forcefield.evalenergy!\n\n\n\n\n\n"
},

{
    "location": "drivers.html#Monte-Carlo-1",
    "page": "Drivers",
    "title": "Monte Carlo",
    "category": "section",
    "text": "CurrentModule = Drivers.MonteCarloThis section provides a description on the Monte Carlo Driver. This Driver iterates over a set amount of steps (defined in the ConfigParameters), sampling new conformations to the Common.State and accepting or rejecting them based on the Metropolis Algorithm.ConfigParameters\nload_parameters\nrun!"
},

{
    "location": "drivers.html#ProtoSyn.Drivers.SteepestDescent.ConfigParameters",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.SteepestDescent.ConfigParameters",
    "category": "type",
    "text": "ConfigParameters(n_steps::Int64 = 0, log_freq::Int64 = 1, f_tol::Float64 = 1e-3, max_step:Float64 = 0.1)\n\nDefine the runtime parameters for the Monte Carlo Driver. If n_steps is zero, a single point energy calculation is performed.\n\nArguments\n\nn_steps: Total amount of steps to be performed (if convergence is not achieved before) (Default: 0).\nlog_freq: Defines the frequency (in steps) of the logs (Default: 1).\nf_tol: Force tolerance. Defines a finalization criteria, as the steepest descent is considered converged if the maximum force calculated is below this value (Default = 1e-3).\nmax_step: Defines the maximum value ɣ that the system can jump when applying the forces (Default: 0.1).\n\nExamples\n\njulia> Drivers.SteepestDescent.ConfigParameters(100, 5, 1e-3, 0.1)\nDrivers.SteepestDescent.ConfigParameters(n_steps=100, log_freq=5, f_tol=1e-3, max_step=0.1)\n\njulia> Drivers.SteepestDescent.ConfigParameters(f_tol = 1e-6)\nDrivers.SteepestDescent.ConfigParameters(n_steps=0, log_freq=1, f_tol=1e-6, max_step=0.1)\n\nSee also: load_parameters\n\n\n\n\n\n"
},

{
    "location": "drivers.html#ProtoSyn.Drivers.SteepestDescent.load_parameters",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.SteepestDescent.load_parameters",
    "category": "function",
    "text": "load_parameters(p::Dict{String, Any})::ConfigParameters\n\nLoad the ConfigParameters from a dictionary.\n\nExamples\n\njulia> Drivers.SteepestDescent.load_parameters(p)\nDrivers.SteepestDescent.ConfigParameters(n_steps=100, log_freq=5, f_tol=1e-3, max_step=0.1)\n\nSee also: Aux.read_JSON\n\n\n\n\n\n"
},

{
    "location": "drivers.html#ProtoSyn.Drivers.SteepestDescent.run!",
    "page": "Drivers",
    "title": "ProtoSyn.Drivers.SteepestDescent.run!",
    "category": "function",
    "text": "run!(state::Common.State, evaluator!::Function, params::ConfigParameters[, ostream::IO = stdout, callback::Union{Function, Nothing} = nothing])\n\nRun the main body of the Driver. If params.n_steps is zero, a single point energy calculation is performed.\n\nArguments\n\nstate::Common.State: Current state of the system to be modified.\nevaluator!::Function: Responsible for evaluating the current state.energy and calculate the resulting forces. This function should have the following signature:\n\nevaluator!(state::Common.State, do_forces::Bool)\n\nparams::ConfigParameters: Hold the runtime parameters of the Steepest Descent Driver. \nostream::IO: (Optional) Any logging will be written to the supplied ostream.\ncallback::Union{Function, Nothing}: (Optional) If present, this function will be called every ConfigParameters.log_freq. This function should have the following signature:\n\ncallback(state::Common.State, step::Int64)\n\ntip: Tip\nThe callback function is often a Print function.\n\nExamples\n\njulia> Drivers.SteepestDescent.run(state, my_evaluator!, params!, callback = my_callback)\n\nSee also: Forcefield.evalenergy!\n\n\n\n\n\n"
},

{
    "location": "drivers.html#Steepest-Descent-1",
    "page": "Drivers",
    "title": "Steepest Descent",
    "category": "section",
    "text": "CurrentModule = Drivers.SteepestDescentThis section provides a description on the Steepest Descent Driver. This Driver attempts to minimize the system energy based on the provided Topology, as it calculates the forces acting on each atom according to the defined ConfigParameters. The finalization criteria is:Maximum number of steps was achieved (ConfigParameters.n_steps).\nMaximum force calculated is below the force tolerance (ConfigParameters.f_tol).ConfigParameters\nload_parameters\nrun!"
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
    "text": "as_xyz(state::Common.State[, ostream::IO = stdout, title::String = \"mol\"])\n\nPrint the current Common.State as a .xyz file.\n\nExamples\n\njulia> Drivers.MonteCarlo.load_parameters(state, title = \"molecule\")\n2\n molecule\n N      -0.0040    0.2990    0.0000\n H1      0.1200    1.3010    0.0000\n\n\n\n\n\n"
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
    "text": "rotation_matrix_from_axis_angle(axis::Vector{Float64}, angle::Float64)\n\nReturn a rotation matrix based on the provided axis and angle (in radians).\n\nExamples\n\njulia> Aux.rotation_matrix_from_axis_angle([1.1, 2.2, 3.3], π/2)\n3×3 Array{Float64,2}:\n  0.0714286  -0.658927  0.748808\n  0.944641    0.285714  0.16131 \n -0.320237    0.695833  0.642857\n\nSee also: Mutators.Dihedral.rotate_dihedral!\n\n\n\n\n\n"
},

{
    "location": "aux.html#Aux-1",
    "page": "Aux",
    "title": "Aux",
    "category": "section",
    "text": "This section provides a description of miscellaneous auxiliary functions.CurrentModule = Auxread_JSON\nrotation_matrix_from_axis_angle"
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
