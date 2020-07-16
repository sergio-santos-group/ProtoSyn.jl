abstract type AbstractSelection end

abstract type AbstractStateMode end

struct Stateful <: AbstractStateMode end
struct Stateless <: AbstractStateMode end

export select

include("masks.jl")
include("promotion.jl")
include("binary.jl")
include("field.jl")
include("unary.jl")
include("cast.jl")
include("true.jl")
include("distance.jl")


# --- Resolve Function ---------------------------------------------------------

function (sele::AbstractSelection)(container::AbstractContainer)
    return select(sele, container)
end

function (sele::AbstractSelection)(container::AbstractContainer, state::State)
    if state_mode_type(sele) == Stateful
        return select(sele, container)(state)
    else
        return select(sele, container)
    end
end

(sele::AbstractSelection)(pose::Pose) = sele(pose.graph, pose.state)

# ---

# --- Polarity
# export PolarSelection
# mutable struct PolarSelection{M, T} <: AbstractSelection
#     is_exit_node::Bool

#     PolarSelection{}() = new{Stateless, Residue}(true)
# end

# # --- Select -------------------------------------------------------------------
# function select(::PolarSelection, container::AbstractContainer)

#     n_residues = count_residues(container)
#     mask = Mask{Residue}(n_residues)

#     for residue in eachresidue(container)
#         if residue.name in ["ARG", "ASN", "ASP", "GLU", "GLN", "HIS", "LYS", "SER", "THR", "TYR"]
#             mask[residue.index] = true
#         end
#     end
#     return mask
# end

# state_mode_type(::PolarSelection{M, T}) where {M, T} = M

export print_selection
function print_selection(io::IOStream, pose::Pose{Topology}, mask::Mask{T}) where {T <: AbstractContainer}

    if selection_type(mask) != Atom
        mask  = promote(mask, Atom, pose.graph)
    end
    
    Base.write(io, "MODEL\n")
    for (atom_index, atom) in enumerate(eachatom(pose.graph))
        sti = pose.state[atom.index] # returns and AtomState instance

        # In this file, selected atoms will be displayed in red while
        # non-selected atoms will be displayed in blue
        atom_symbol = mask[atom_index] ? "O" : "N"

        s = @sprintf("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%24s\n",
            atom.index, atom_symbol,
            atom.container.name, atom.container.container.code,
            atom.container.id,
            sti.t[1], sti.t[2], sti.t[3],
            atom_symbol)
            Base.write(io, s)
    end

    for atom in eachatom(pose.graph)
        Base.write(io, @sprintf("CONECT%5d", atom.index))
       foreach(n->Base.write(io, @sprintf("%5d",n.index)), atom.bonds)
       Base.write(io,"\n")
    end
    Base.write(io, "ENDMDL")
end