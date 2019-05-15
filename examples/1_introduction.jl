using ProtoSyn
using Printf

#= -------------------------------------------------
1. Introduction Example

# Algorithm explanation:
ProtoSyn is a suite for molecular modeling, in specific, proteins.
Using a modular approach, ProtoSyn changes the conformation of a State
using Mutators, and new conformations are evaluated using Forcefield Components.
These structures can be combined in any custom variation, altough
ProtoSyn has default Drivers who carry out useful algorithms in an optimized
and reliable way. In the following examples, all these concepts will be
introduced and explained in more detail.

In this example, a simple Dihedral rotation will be performed on a structure
loaded from a PDB file.

# New ProtoSyn concepts:
- State & Metadata
- Mutators
- Printing
-------------------------------------------------=#

# Configuration
input_pdb                = "data/1i2t_no_sc.pdb"
ss                       = "CHHHHHHHHHHHHHHHCCCCHHHHHHHHHHCCHHHHHHHHHCHHHHHHHHHHHHHHHHHHC"

# State & Metadata
#=
* Note: On the basis of ProtoSyn lies a State. A State contains all information
relative to the current configuration of the molecule or system, such as the
atoms coordinates, forces being applied, energy, etc. A state is passed around
functions who carry actions over them: evaluators update the systems energy,
mutators change the coordinates in specific ways, etc. Usually a State is loaded
from a PDB file. Accompanying a State is usually a Metadata structure. This
structure contains pertinent static data that complements the State: 
a coordinate point can be an atom with a name, index, element, etc. This
information, however, is mostly not optional, as some functions require to know
additional information about the system to perform correctly: a Dihedral
Mutator requires a list of dihedrals (defined in the metadata) to perform
the movement, as an example.
=#
state, metadata = Common.load_from_pdb(input_pdb)
Common.apply_ss!(state, metadata, ss)

# Mutators
#=
* Note: A Mutator is a structure responsible from changing the State
coordinates in a smart way. Different Mutators have been made available
in ProtoSyn, and are explained in detail on the package documentation.
In this example, a simple Dihedral movement will be performed. Depending
on the nature of the Mutator, the required parameters may change. In this
case, only the dihedrals belonging to coil secondary structures of the
protein will be rotate.. Moreover, every Mutator has different configuration
parameters such as the probabily of mutation or step size.
=#
coil_dihedrals = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
dihedral_mutator   = Mutators.Dihedral.MutatorConfig(
    dihedrals = coil_dihedrals,
    angle_sampler = () -> (randn() * dihedral_mutator.step_size),
    p_mut = 0.025,
    step_size = π/8)


# Printing
#=
* Note: Several default functions have been made available in the Print module
of ProtoSyn, several of which allow to combine the State and Metadata information
and write structural files in useful formats (PDB, GRO, XYZ, etc)
=#
Print.as_pdb(output_file, state, metadata)

if ""!=PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    const output_file = open("1_introduction.pdb", "w")

    Mutators.apply!(state, dihedral_mutator)
    Print.as_pdb(output_file, state, metadata)
    
    close(output_file)
end