export AbstractSelection
abstract type AbstractSelection end

abstract type AbstractStateMode end

struct Stateful <: AbstractStateMode end  # * Needs state to get selection
struct Stateless <: AbstractStateMode end # * Doesn't need state to get selection

export select

include("masks.jl")
include("promotion.jl")
include("binary.jl")
include("field.jl")
include("serial.jl")
include("unary.jl")
include("true.jl")
include("distance.jl")
include("random.jl")
include("terminal.jl")
include("aromatic.jl")
include("bond.jl")
include("bond-count.jl")
include("charge.jl")


# --- Resolve Function ---------------------------------------------------------
function (sele::AbstractSelection)(container::AbstractContainer; gather::Bool = false)

    @assert selection_type(sele) <= typeof(container) "Can't select $(selection_type(sele))s from a container of type $(typeof(container))."

    mask = select(sele, container)
    if gather
        return ProtoSyn.gather(mask, container)
    end
    return mask
end

function (sele::AbstractSelection)(container::AbstractContainer, state::State; gather::Bool = false)

    @assert selection_type(sele) <= typeof(container) "Can't select $(selection_type(sele))s from a container of type $(typeof(container))."

    if state_mode_type(sele) == Stateful
        mask = select(sele, container)(state)
        if gather
            return ProtoSyn.gather(mask, container)
        end
        return mask
    else
        mask = select(sele, container)
        if gather
            return ProtoSyn.gather(mask, container)
        end
        return mask
    end
end

(sele::AbstractSelection)(pose::Pose; gather::Bool = false) = sele(pose.graph, pose.state; gather = gather)


# --- Gather Function ----------------------------------------------------------
"""
    ProtoSyn.gather(mask::Mask{T}, container::AbstractContainer) where {T <: AbstractContainer}

Gather all instances of type `T` from `container` whose relative position is
marked as `true` in the given `mask`.

# Examples
```jldoctest
julia> ProtoSyn.gather(rn"ALA"(pose), pose.graph)
4-element Vector{Residue}:
 Residue{/UNK:1/UNK:1/ALA:5}
 Residue{/UNK:1/UNK:1/ALA:12}
 Residue{/UNK:1/UNK:1/ALA:13}
 Residue{/UNK:1/UNK:1/ALA:20}
```
"""
function gather(mask::Mask{T}, container::AbstractContainer) where {T <: AbstractContainer}
    results = Vector{T}()
    for (index, item) in enumerate(iterator(T)(container))
        if mask[index]
            push!(results, item)
        end
    end
    
    return results
end


export print_selection
function print_selection(pose::Pose{Topology}, mask::Mask{T}, filename::String) where {T <: AbstractContainer}

    io = open(filename, "w")

    if selection_type(mask) != Atom
        mask  = promote(mask, Atom, pose.graph)
    end
    
    Base.write(io, "MODEL\n")
    for (atom_index, atom) in enumerate(eachatom(pose.graph))
        sti = pose.state[atom.index] # returns and AtomState instance

        # In this file, selected atoms will be displayed in chain A while
        # non-selected atoms will be displayed in chain B
        chain = mask[atom_index] ? "A" : "B"

        s = @sprintf("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%24s\n",
            atom.index, atom.name,
            atom.container.name, chain,
            atom.container.id,
            sti.t[1], sti.t[2], sti.t[3],
            atom.symbol)
            Base.write(io, s)
    end

    for atom in eachatom(pose.graph)
        Base.write(io, @sprintf("CONECT%5d", atom.index))
       foreach(n->Base.write(io, @sprintf("%5d",n.index)), atom.bonds)
       Base.write(io,"\n")
    end
    Base.write(io, "ENDMDL")
    close(io)
end