@doc raw"""
    NonBondedList(cut_off::Float64, buffer::Float64, list::Vector{Int64}, pointer::Vector{Int64})

Define the current list of interacting atoms through non bonded interactions.

# Arguments
- `cut_off::Float64`: Distance radius (in nm), from one given atom, where other atoms interactions are accounted for.
- `buffer::Float64`: Distance radius (in nm) that is ADDED to the `cut_off`, calculated based on the maximum movement of the atom in the given steps until update of the nblists.
- `list::Vector{Int64}`: List of all atoms that interact with every other atom.
- `pointer::Vector{Int64}`: List of all atoms considered for non-bonded interactions.

The `list` has all atoms in a single column, while `pointer` holds the index values at `list` where a new atoms neightbours are to be considered.
For example, is `list` has [1, 2, 3, 4, 5] and pointer has [1, 3], atom A's neighbours are [1, 2], atom B's neighbours are [3, 4, 5].

# Examples
```julia-repl
julia> Common.NonBondedList(3)
NonBondedList(list=[0, 0, 0, 0, 0, 0], pointer=[0, 0, 0])
   pointer size = 3 items
   list size = 6 items
   cut_off = -1.0 nm
   buffer = 0.0 nm
```
See also: [`State`](@ref)
"""
mutable struct NonBondedList
    cut_off::Float64
    buffer::Float64
    list::Vector{Int64}
    pointer::Vector{Int64}
end

function NonBondedList(n::Int)
    list_size = convert(Int, n*(n+1)/2)
    NonBondedList(-1.0, 0.0, zeros(Int, list_size), zeros(n))
end

function Base.show(io::IO, b::NonBondedList)
    print(io, "NonBondedList(list=$(b.list), pointer=$(b.pointer))")
    print(io, "\n   pointer size = $(length(b.pointer)) items")
    print(io, "\n   list size = $(length(b.list)) items")
    print(io, "\n   cut_off = $(b.cut_off) nm")
    print(io, "\n   buffer = $(b.buffer) nm")
end
 
function Base.copy!(dst::NonBondedList, src::NonBondedList)
    copy!(dst.pointer, src.pointer)
    copy!(dst.list, src.list)
    dst.cut_off = src.cut_off
    dst.buffer = src.buffer
    return dst
end


@doc raw"""
    State(size::Int64, energy::Energy, xyz::Array{Float64, 2}, forces::Array{Float64, 2}, metadata::Vector{Metadata})

Define the current state of the system, containing the atoms positions, energy and forces applied.
If only `size::Int64` is provided, an empty State with the given `size` is created with zeros.

# Arguments
- `size::Int64`: Atom count in system.
- `energy::Energy`: Current energy of the system (kJ mol⁻¹).
- `xyz::Array{Float64, 2}`: Atom positions in 3 dimensions.
- `forces::Array{Float64, 2}`: Forces applied in each dimension to each atom (kJ mol⁻¹ nm⁻¹)
- `nblist::Union{NonBondedList, Nothing}`: List of non bonded atoms interacting acvcording to current conditions.

# Examples
```julia-repl
julia> Common.State(3)
Common.State(size=3, energy=Null, xyz=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0], forces=[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])

julia> Common.State(2, Common.Energy(), [1.1 1.1 1.1; 2.2 2.2 2.2], zeros(2, 3))
Common.State(size=2, energy=Null, xyz=[1.1 1.1 1.1; 2.2 2.2 2.2], forces=[0.0 0.0 0.0; 0.0 0.0 0.0])
```
See also: [`Metadata`](@ref)
"""
mutable struct State
    size::Int64
    energy::Energy
    xyz::Array{Float64, 2}
    forces::Array{Float64, 2} # kJ mol⁻¹ nm⁻¹
    nblist::Union{NonBondedList, Nothing}
end
State(n::Int64, e::Energy, x::Array{Float64, 2}, f::Array{Float64, 2}) = State(n, e, x, f, NonBondedList(n))
State(n::Int64) = State(n, Energy(), zeros(n, 3), zeros(n, 3), NonBondedList(n))

function State(state::State)
    st = State(state.size)
    copy!(st, state)
end

Base.show(io::IO, b::State) = print(io, "State(size=$(b.size), energy=$(b.energy), xyz=$(b.xyz), forces=$(b.forces), nblist=$(b.nblist)")

function Base.iterate(st::State, idx = 1)

    if idx > size(st.xyz, 1)
        return nothing
    else
        return (st.xyz[idx, :], idx+1)
    end
end

function Base.copy!(dst::State, src::State)::State
    copy!(dst.xyz, src.xyz)
    copy!(dst.forces, src.forces)
    copy!(dst.energy, src.energy)
    copy!(dst.nblist, src.nblist)
    return dst
end


@doc raw"""
    update_nblist!(state::Common.State)::Common.State

Update the `state` non bonded list.

# Examples
```julia-repl
julia> Common.update_nblist!(state)
State(...)
```
See also: [`State`](@ref) [`NonBondedList`](@ref)
"""
function update_nblist!(state::Common.State)::Common.State
    # if no NonBondedList is attached to the input state,
    # the do nothing
    if state.nblist == nothing
        return state
    end
    
    nblist  = state.nblist
    n_atoms = state.size
    coords  = state.xyz

    if nblist.cut_off < 0.0
        # if all pairwise interactions are requested
        # (cut_off < 0), construction of the nblist does
        # not require any distance computations
        for i=1:n_atoms
            nblist.pointer[i] = i
            nblist.list[i] = i+1
        end
        nblist.list[n_atoms] = -1
    else
        # otherwise, the effective cut_off is the
        # sum of the cuttoff and the buffer.
        ptr = 1
        δx = 0.0
        Δx = 0.0
        cutSq = (nblist.cut_off + nblist.buffer)^2
        for i=1:n_atoms
            nblist.pointer[i] = ptr
            for j = (i+1):n_atoms
                Δx = 0.0
                for k=1:3
                    δx = coords[i,k] - coords[j,k]
                    Δx += δx^2
                end
                if Δx < cutSq
                    nblist.list[ptr] = j
                    ptr += 1
                end
            end
            nblist.list[ptr] = -1
            ptr += 1
        end
    end # endif

    return state
end # end fcn
