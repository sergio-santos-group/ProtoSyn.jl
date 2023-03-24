# This file should contain functions that work on the system Graph, such as 
# functions that deal with parent/child relations, indexes and bonds, among
# others.

export origin
"""
    origin(container::AbstractContainer)

Return the first `Atom` in `AbstractContainer` `container` that has no parent.
The iteration follows the [`Atom`](@ref) instance `:id` field, if correctly
indexed. If no [`Atom`](@ref) instance without parent is found (i.e.: circular
structures), return `nothing`. Note that the [`root`](@ref) atoms are not
considered. Note that if multiple origin [`Atom`](@ref) instances exists, this
method return only the first found, based on the current [`Atom`](@ref) order
in the given `AbstractContainer` `container`.

# See also
[`root`](@ref) [`reindex`](@ref) [`sort_atoms_by_graph!`](@ref)

"""
origin(container::AbstractContainer) = begin
    root = ProtoSyn.root(container)
    for atom in eachatom(container)
        !hasparent(atom) && return atom
        atom.parent === root && return atom
    end
    
    return nothing
end


export root
"""
    root(container::AbstractContainer)

Return the first [`Atom`](@ref) of the Root of the [Graph](@ref state-types) that given
`AbstractContainer` `container` belongs to. If the given `AbstractContainer`
`container` is not a [`Topology`](@ref) instance and has `:container` field set
to `nothing`, return `nothing`.

    root(topology::Topology)

Return the first [`Atom`](@ref) of the Root of the given [`Topology`](@ref)
`topology` instance.

# See also
[`origin`](@ref)
"""
@inline root(topology::Topology) = get(topology.root, "OO")
@inline root(container::AbstractContainer) = begin
    hascontainer(container) ? root(container.container) : nothing
end


export hasparent
"""
    hasparent(c::AbstractContainer) -> Bool

Test whether the given AbstractContainer `c` has a parent.

# See Also
[`isparent`](@ref)
"""
@inline hasparent(c::AbstractContainer) = c.parent !== nothing


export isparent
"""
    isparent(parent::AbstractContainer, child::AbstractContainer)

Test whether `parent` is the parent of `child`.

# See Also
[`hasparent`](@ref)
"""
@inline isparent(parent::AbstractContainer, child::AbstractContainer) = begin
    parent === child.parent
end

@inline isparent(::Nothing, child::AbstractContainer) = false


export setparent!
"""
    setparent!(child::T, parent::T) where {T <: AbstractContainer}

Set `parent` as the parent of `child`, while adding `child` to
`parent.children`.

# See also
[`popparent!`](@ref)
"""
function setparent!(child::T, parent::T) where {T <: AbstractContainer}
    @info "Setting $parent as the parent of $child ..."
    hasparent(child) && begin
        error("unable to setparent! of non-orphan item")
    end
    push!(parent.children, child)
    child.parent = parent
    child
end


export popparent!
"""
    popparent!(child::AbstractContainer}

Remove the parent from `child` (sets it to `nothing`) while removing `child`
from `parent.children` (only if `child` is a child of `parent`).

# See also
[`setparent!`](@ref)
"""
function popparent!(child::AbstractContainer)
    @info "Popping the parent of $child ..."
    if hasparent(child)
        parent = child.parent
        i = findfirst(x -> x === child, parent.children)
        if i !== nothing
            deleteat!(parent.children, i)
            child.parent = nothing
        end
    else
        @warn "Tried to pop parent of $child but no parent was found."
    end
    child
end


export hascontainer
"""
    hascontainer(c::AbstractContainer)

Return `true` if the given `AbstractContainer.container` is not `nothing`.
"""
@inline hascontainer(c::AbstractContainer) = c.container !== nothing


export genid
"""
    genid()

Return a random `UInt16` number.
"""
@inline genid() = Int(rand(UInt16))


export reindex
"""
    reindex(topology::Topology; set_ascendents = true)
    
Re-indexes the whole [`Topology`](@ref) `topology`, setting both the `:id` and
`:index` of instances inside the `topology` to the corresponding relative index
in the `container.items` which they belong to. If `set_ascendents` is set to
`true` (is, by default), each [`Atom`](@ref) instance `:ascendents` field will
be updated to reflect the new indices.

    reindex(segment::Segment)

Re-indexes a [`Segment`](@ref) `segment`, setting both the `:id` and `:index` of
instances inside the `topology` to the corresponding relative index in the
`container.items` which they belong to.

# See also
[`ascendents`](@ref) [`reindex(::State)`](@ref)

# Examples
```jldoctest
julia> reindex(pose.graph)
Topology{/UNK:1}

julia> reindex(pose.graph[1])
Segment{/UNK:1/UNK:1}
```
"""
function reindex(topology::Topology; set_ascendents::Bool = true)
    @info "Re-indexing topology. Set ascendents = $(set_ascendents)."
    aid = rid = sid = 0
    for segment in topology.items
        segment.id    = (sid += 1)
        segment.index = sid
        for residue in segment.items
            residue.id    = (rid += 1)
            residue.index = rid
            for atom in residue.items
                atom.id    = (aid += 1)
                atom.index = aid
            end
        end
    end

    # update ascendents (not possible before because
    # of possible problems with index assignment)
    if set_ascendents
        for atom in eachatom(topology)
            # println("\n Setting ascendents on atom: $atom")
            atom.ascendents = ascendents(atom, 4)
        end
    end

    return topology
end

function reindex(segment::Segment)
    aid = rid = 0
    for residue in segment.items
        residue.id    = (rid += 1)
        residue.index = rid
        for atom in residue.items
            atom.id    = (aid += 1)
            atom.index = aid
        end
    end
    return segment
end


export ascendents
"""
    ascedents(container::AbstractContainer, level::Int)
    
Return a `Tuple` containing the N (`level`) previous `:id` fields of the
`:parent` `AbstractContainer` instances of the given `container` (recursivelly).

# Examples
```jldoctest
julia> ascendents(pose.graph[1][1][4], 4)
(4, 3, 1, 0)
```
"""
function ascendents(container::AbstractContainer, level::Int)
    # println("$container (Level $level) (Next: $(container.parent))")
    if level > 1
        @assert container.parent !== nothing "$container has no parent."
        return (container.index, ascendents(container.parent, level - 1)...)
    else
        (container.index,)
    end
end


"""
    rename!(atom::Atom, name::String; force_rename::Bool = false)

Rename the selected [`Atom`](@ref) instance to the given `name`. Also updates
the [`Atom`](@ref) container `:itemsbyname` field. If `force_rename` is set to
`true` (`false`, by default), will change existing [`Atom`](@ref) instances with
the same `name` trying to be introduced to adopt a temporary name (current name
with \"_o\" appendix).

# Examples
```
julia> ProtoSyn.rename!(pose.graph[1][1]["N"], "N1")
Atom{/2a3d:51894/A:1/MET:1/N1:1}
```
"""
function rename!(atom::Atom, name::String; force_rename::Bool = false)
    similar = findfirst(atom -> atom.name === name, atom.container.items)

    # Change to same name
    similar !== nothing && atom.container.items[similar] === atom && return atom 

    if !(similar === nothing)
        if force_rename
            @info "Tried to rename atom $atom to $name, but another atom in the same Residue was found with that name ($(atom.container.items[similar])). Changing existing atom to temporary placeholder."
            similar_atom = atom.container.items[similar]
            ProtoSyn.rename!(similar_atom, "$(similar_atom.name)_o") # Temporary
        else
            @assert similar === nothing "Tried to rename atom $atom to $name, but another atom in the same Residue was found with that name ($(atom.container.items[similar]))"
        end
    end

    if !force_rename
        pop!(atom.container.itemsbyname, atom.name)
    end
    atom.name = name
    atom.container.itemsbyname[name] = atom

    return atom
end


export unbond!
"""
    unbond!(pose::Pose, at1::Atom, at2::Atom; [keep_downstream_position::Bool = true])::Pose
    
Return a [Pose](@ref pose-types) instance with both given [`Atom`](@ref) instances unbonded
(removed from eachother `bonds` list, pops parenthood and sets the downstream
[`Residue`](@ref)`.parent` field to be the Root of the upstream
[`Topology`](@ref)). If `keep_downstream_position` is set to `true` (is, by
default), the downstream [`Residue`](@ref) position is maintained (by calling
[`request_c2i!`](@ref) and [`sync!`](@ref) methods). 

!!! ukw "Note:"
    Unbonding two atoms also removes any parenthood relationship, therefore
    making the returned [`Pose`](@ref) from this function un-usable without
    further changes (the internal coordinates graph is severed on the unbonding
    site).

# Examples
```jldoctest
julia> unbond!(pose, pose.graph[1][2]["C"], pose.graph[1][3]["N"])
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function unbond!(pose::Pose, at1::Atom, at2::Atom; keep_downstream_position::Bool = true)::Pose
    @assert (at2 in at1.bonds) & (at1 in at2.bonds) "Atoms $at1 and $at2 are not bonded and therefore cannot be unbonded."
    isparent(at1, at2) && return _unbond!(pose, at1, at2; keep_downstream_position = keep_downstream_position)
    isparent(at2, at1) && return _unbond!(pose, at2, at1; keep_downstream_position = keep_downstream_position)
    
    # The two atoms might be bonded but not have a parenthood relationship, in
    # which case we just remove eachother from the bond list and return the pose
    @info " Unbonding $at1 <-> $at2 (With NO parenthood relationships)"
    i = findfirst(at1, at2.bonds)
    i !== nothing && deleteat!(at2.bonds, i)
    j = findfirst(at2, at1.bonds)
    j !== nothing && deleteat!(at1.bonds, j)
    return pose
end

function _unbond!(pose::Pose, at1::Atom, at2::Atom; keep_downstream_position::Bool = true)::Pose
    
    @info " Unbonding $at1 <-> $at2"

    # Following this two step example:
    # 1) ProtoSyn.setdihedral!(pose.state, pose.graph[1][10]["CG"], 90°)
    # 2) ProtoSyn.unbond!(pose, pose.graph[1][10]["CG"],pose.graph[1][10]["CB"])
    # If the pending internal to cartesian changes were not synched at the
    # beginning of `unbond!`, the saved position (if `keep_downstream_position`
    # is set to true) would be the pre-rotation position.

    # Sync any pending internal to cartesian changes (except if it's a fragment)
    if keep_downstream_position && !isfragment(pose)
        sync!(pose)
    end
    
    i = findfirst(at1, at2.bonds)
    i !== nothing && deleteat!(at2.bonds, i)
    
    j = findfirst(at2, at1.bonds)
    j !== nothing && deleteat!(at1.bonds, j)
    
    # Detach from atom graph
    # Remove at2 from at1.children and set at2.parent to be the root of at1
    # (only if this this not a fragment)
    # This assumes at1 is parent of at2
    popparent!(at2)
    !isfragment(pose) && begin
        _root = ProtoSyn.root(at1.container.container.container)
        at1.container.container.container !== nothing && begin
            setparent!(at2, _root)
        end
        @info " $at2 parent is now Root."
        reindex(pose.graph) # To set new ascendents
    end

    # ? Should `unbond!` remove Residue level parenthood relationships in
    # ? inter-residue bonds, and connect the downstream residue to Root?
    if at1.container !== at2.container
        ProtoSyn.popparent!(at2.container)
        ProtoSyn.setparent!(at2.container, _root.container)
    end
    
    if keep_downstream_position

        ProtoSyn.request_c2i!(pose.state, all = true)
        !isfragment(pose) && begin
            @info " Fixating cartesian positions."
            sync!(pose)
        end
    end

    ProtoSyn.request_i2c!(pose.state)    
    return pose
end


"""
    bond(at1::Atom, at2::Atom)
    
Bond both given [`Atom`](@ref) instances (adds `at2` to `at1.bonds` and
vice-versa). Both [`Atom`](@ref) instances need to be in the same
[`Segment`](@ref).

# See also
[`join`](@ref) [`unbond!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.bond(pose.graph[1][1]["C"], pose.graph[1][2]["CA"])
```
"""
@inline function bond(at1::Atom, at2::Atom)
    @assert at1.container.container === at2.container.container "Can only bond atoms within the same segment.\n Tried to bond $at1 <-> $at2"
    !in(at2, at1.bonds) && push!(at1.bonds, at2)
    !in(at1, at2.bonds) && push!(at2.bonds, at1)
    return nothing
end


export join
"""
    join(at1::Atom, at2::Atom)

Join [`Atom`](@ref) `at1` with [`Atom`](@ref) `at2`.

    join(r1::Residue, s1::String, r2::Residue, s2::String)
    
Join [`Atom`](@ref) named `s1` from [`Residue`](@ref) `r1` with [`Atom`](@ref)
named `s2` from [`Residue`](@ref) `r2`.

Bond (add eachother to `other.bonds` field) and set parent/children relationship
of both the [`Atom`](@ref) instances and respective `atom.container`
([`Residue`](@ref)). Note that `at2` [`Atom`](@ref) will become parent at `at1`
(and `at2.container` [`Residue`](@ref) will become parent of `at1.container`).

# See also
[`bond`](@ref) [`unbond!`](@ref)

# Examples
```jldoctest
julia> Residue!(pose.graph[1], ProtoSyn.ResidueName("ALA"), 1);

julia> Atom!(pose.graph[1][end], "N", 1, 1, "N");

julia> ProtoSyn.join(pose.graph[1][1], "C", pose.graph[1][end], "N")
Residue{/UNK:1/UNK:1/ALA:1}
```
"""
function join(r1::Residue, s1::String, r2::Residue, s2::String) # ! IMPORTANT
    hasparent(r2) && error("r2 is already connected.")
    at1 = r1[s1]
    at2 = r2[s2]
    bond(at1, at2)          # at1 <-> at2
    setparent!(at2, at1)    # at1 -> at2
    setparent!(r2, r1)      # r1 -> r2
end

function join(at1::Atom, at2::Atom)
    hasparent(at2) && error("at2 is already connected.")
    bond(at1, at2)
    setparent!(at2, at1)
    setparent!(at2.container, at1.container)
end

# COUNTERS
count_segments(t::Topology) = length(t.items)

count_residues(c::AbstractContainer) = mapreduce(x -> count_residues(x), +, c.items; init=0)
count_residues(s::Segment) = length(s.items)
count_residues(r::Residue) = 1

count_atoms(c::AbstractContainer) = mapreduce(x -> count_atoms(x), +, c.items, init=0)
count_atoms(r::Residue) = r.size
count_atoms(a::Atom) = 1


include("travel-graph.jl")

export ids
"""
    ids(atoms::Vector{Atom})

Return a vector with the `:id` `Int` field for every [`Atom`](@ref) in the given
`atoms` vector.

# See also
[`travel_graph`](@ref)

# Examples
```
julia> ProtoSyn.ids(an"CA"(pose, gather = true))
21-element Vector{Int64}:
   3
  14
  29
  40
  55
  65
   ⋮
 243
 257
 281
 300
 317
 327
```
"""
function ids(atoms::Vector{Atom})::Vector{Int}
    idxs = Vector{Int}()
    for atom in atoms
        push!(idxs, atom.id)
    end
    return idxs
end


"""
    is_contiguous(pose::Pose, selection::AbstractSelection)

Returns `true` if all the [`Residue`](@ref) instances gathered from the
`selection` applied to the given `pose` are contiguous (have a parenthood
relationship connecting them all). Note that the given `selection` is always
promoted to [`Residue`](@ref) level.

# See also
[`ProtoSyn.promote`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.is_contiguous(pose, rid"1" | rid"3")
false

julia> ProtoSyn.is_contiguous(pose, rid"1:10")
true
```
"""
function is_contiguous(pose::Pose, selection::ProtoSyn.AbstractSelection)
    sele              = promote(selection, Residue)
    selected_residues = sele(pose, gather = true)
    root              = ProtoSyn.root(pose.graph).container # * Residue
    
    # Initialize
    stack = Vector{Residue}()
    for child in root.children
        push!(stack, child)
    end
    
    # Travel graph, search for first selected residue
    selected_stack = Vector{Residue}()
    while length(stack) > 0
        residue = pop!(stack)
        if residue in selected_residues
            push!(selected_stack, residue)
            break
        end
        for child in residue.children
            push!(stack, child)
        end
    end

    length(selected_stack) == 0 && begin
        @warn "Finished travelling the Pose.graph but the selected residues were not found."
        return nothing
    end

    # Travel graph, mark selected atoms found
    marks = 0
    while length(selected_stack) > 0
        marks += 1
        residue = pop!(selected_stack)
        for child in residue.children
            if child in selected_residues
                push!(selected_stack, child)
            end
        end
    end

    # Compare number of found marks with the number of selected residues
    return marks == length(selected_residues)
end


"""
    infer_parenthood!(container::ProtoSyn.AbstractContainer; overwrite::Bool = false, start::Opt{Atom} = nothing)

Infers parenthood of [`Atom`](@ref) instances on the given `AbstractContainer`
`container`, from bond information, using a custom algorithm similar to breath
first algorithm (atoms are sorted based on the size of the downstream graph and
aromaticity). By default, the [Graph](@ref graph-types) origin is set to the first
[`Atom`](@ref) instance in the `container`. This behaviour can be controlled by
setting a `start` [`Atom`](@ref) as the origin of the new infered parenthood
[Graph](@ref graph-types). If `overwrite` is set to `true` (`false`, by default), will
overwrite existing pranthood information. After infering parenthood, if changes
to the [Graph](@ref graph-types) occurred, the existing internal coordinates match different
cartesian coordinates. It's suggested to update internal coordinates
([`request_i2c!`](@ref ProtoSyn.request_i2c!) &
[`sync!`](@ref ProtoSyn.sync!)). For more details, see the
[Travelling the Graph](@ref) section. If the `linear_aromatics` flag is set to
`true` (is, by default), aromatic rings are treated as isolated structures is an
otherwise linear [Graph](@ref graph-types) (for example, in some protein aminoacids). More
complex structures (such as carbon sheets) have interlaced aromatic rings, and
the `linear_aromatics` should be set to `false` to ensure all [`Atom`](@ref)
instances are visited.

# See also
[`travel_bonds`](@ref) [`infer_bonds!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.infer_parenthood!(pose.graph[1], overwrite = true)
Segment{/2a3d:8625/A:1}
```
"""
function infer_parenthood!(container::ProtoSyn.AbstractContainer; overwrite::Bool = false, start::Opt{Atom} = nothing, linear_aromatics::Bool = true)

    w = Dict{String, Int}("H" => 1, "N" => 5,  "O" => 10, "C" => 15, "S" => 20, "FE" => 25)

    @assert typeof(container) > Atom "Can't infer parenthood of a single Atom!"

    atoms   = collect(eachatom(container))

    aromatics = AromaticSelection()(container)

    if start === nothing
        start = collect(eachatom(container))[1]
    end

    n_atoms = length(atoms)
    visited = ProtoSyn.Mask{Atom}(n_atoms)
    changed = ProtoSyn.Mask{Atom}(n_atoms)

    # First atom is always connected to root
    _root = ProtoSyn.root(container)
    if overwrite && hasparent(start)
        ProtoSyn.popparent!(start)
    end
    ProtoSyn.setparent!(start, _root)

    # Note the usage of breath-first search (BFS) algorithm. This is the most
    # chemically correct was of attributing parenthood relationships, in
    # accordance with generally accepted conventions.
    stack = Vector{Atom}([start])
    
    while length(stack) > 0

        atom_i = pop!(stack)
        i      = findfirst((a) -> a === atom_i, atoms)
        @info "Focus on atom $atom_i - $i"
        visited[i] = true

        atom_i_bonds = Vector{Atom}()
        for atom_j in atom_i.bonds
            atom_j_index = findfirst((a) -> a === atom_j, atoms)
            if atom_j_index !== nothing && !visited[atom_j_index]
                push!(atom_i_bonds, atom_j)
            end
        end

        if length(atom_i_bonds) === 0
            @info " $atom_i has no non-visited bonds."
            continue
        end

        atom_i_bond_weight = [w[a.symbol] for a in atom_i_bonds]
        atom_i_arom_weight = Vector{Int}()
        for atom_j in atom_i_bonds
            atom_j_index = findfirst((a) -> a === atom_j, atoms)
            if atom_j_index === nothing
                aw = 0
            else
                aw = aromatics[atom_j_index] ? 25 : 0
            end
            push!(atom_i_arom_weight, aw)
        end
        atom_i_weight      = atom_i_bond_weight .+ atom_i_arom_weight
        @info "Atom: $atom_i"
        @info " Non-visited Bonds: $([x.name for x in atom_i_bonds])"
        @info "     Symbol weight: $atom_i_bond_weight"
        @info "   Aromatic weight: $atom_i_arom_weight"

        # Case exists same symbol & same aromaticity
        atom_i_weight_set = collect(Set(atom_i_weight))
        if length(atom_i_weight_set) < length(atom_i_weight)
            for weight in atom_i_weight_set
                eq_weight = findall((w) -> w === weight, atom_i_weight)
                length(eq_weight) === 1 && continue
                eq_weight_atoms = atom_i_bonds[eq_weight]

                # Solve equalization with size of downstream bond graph
                equal_downstream = true
                range            = 1
                N                = length(eq_weight_atoms)
                lengths          = Vector{Int}()
                max_range        = 5
                while equal_downstream
                    lengths = Vector{Int}() 
                    for eq_weight_atom in eq_weight_atoms
                        atoms_in_range = ProtoSyn.travel_bonds(eq_weight_atom, range, atom_i)
                        @info "For range $range, atoms in range of $eq_weight_atom : ($(length(atoms_in_range))) $atoms_in_range"
                        push!(lengths, length(atoms_in_range))
                    end
                    
                    equal_downstream = false
                    for i in 1:(N-1)
                        for j in i+1:N
                            if lengths[i] === lengths[j]
                                equal_downstream = true
                                break
                            end
                        end
                    end

                    range += 1
                    range > max_range && break
                end

                # If first one is, all equal weighted atoms are
                # Aromatic     : longer first
                # Non-aromatic : shorter first
                id = findfirst((a) -> a.id === eq_weight_atoms[1].id, atoms)
                @info "          Aromatic: $(aromatics[id]) (ID: $id)"
                @info "         EQ weight: $eq_weight"
                @info "     Final lengths: $lengths"
                if all(y->y==lengths[1], lengths)
                    permvec = sortperm([x.name for x in atom_i_bonds[eq_weight]], rev = false)
                    @info "     ... Sorting by name ..."
                else
                    permvec = sortperm(lengths, rev = aromatics[id])
                end
                @info "       Size weight: $(["$sw ($(x.name))" for (sw, x) in zip(permvec, atom_i_bonds[eq_weight])])"
                atom_i_weight[eq_weight] .+= permvec
            end
        end

        permvec = sortperm(atom_i_weight, rev = false)
        atom_i_bonds = atom_i_bonds[permvec]
        @info repeat("-", 40)
        @info "      Final weight: $atom_i_weight"
        @info "       Final order: $([a.name for a in atom_i_bonds]) (Note: Atoms are consumed in reverse.)\n"

        for atom_j in atom_i_bonds
            j = findfirst(atom -> atom === atom_j, atoms)
            j !== nothing && @info " Bonded to atom $atom_j - $j (Visited? => $(visited[j]))"
            if j !== nothing && !(visited[j]) && !(changed[j])

                atom_j_index = findfirst((a) -> a === atom_j, atoms)
                if linear_aromatics && aromatics[atom_j_index]
                    push!(stack, atom_j)
                    @info "Linear aromatics! Pushing $atom_j to the top of the stack."

                    if overwrite && hasparent(atom_j)
                        ProtoSyn.popparent!(atom_j)
                    end
                    if !hasparent(atom_j)
                        ProtoSyn.setparent!(atom_j, atom_i)
                        changed[j] = true
                    else
                        @warn "  Tried to set $atom_i as parent of $atom_j, but $atom_j already had parent"
                    end

                    break
                end
                insert!(stack, 1, atom_j)

                if overwrite && hasparent(atom_j)
                    ProtoSyn.popparent!(atom_j)
                end
                if !hasparent(atom_j)
                    ProtoSyn.setparent!(atom_j, atom_i)
                    changed[j] = true
                else
                    @warn "  Tried to set $atom_i as parent of $atom_j, but $atom_j already had parent"
                end
            end
        end
    end

    # Any non-connected atoms (by bonds) are children of root
    for atom in atoms
        !hasparent(atom) && ProtoSyn.setparent!(atom, _root)
    end

    if typeof(container) > Residue
        reindex(container) # Doesn't set ascedents
    end

    @info "After infering parenthood, if changes to the Graph occured, the existing internal coordinates match different cartesian coordinates. It's suggested to update internal coordinates (request_c2i & sync!)."

    return container
end


"""
    infer_bonds!(pose::Pose; [threshold::T = 0.1]) where {T <: AbstractFloat}

Infers bonds for all [`Atom`](@ref) instances of the given [`Pose`](@ref). A new
bond is assigned when a pair of [`Atom`](@ref) instances are within a given
distance, as defined in `ProtoSyn.Units.bond_lengths`. The `threshold` value is
multiplied by the standard bond distance and added to the comparison value to
allow some leeway in the bond distance (0.1, by default).

# See also
[`travel_bonds`](@ref) [`infer_parenthood!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.infer_bonds!(pose)
Pose{Topology}(Topology{/CRV:54976}, State{Float64}:
 Size: 201
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function infer_bonds!(pose::Pose; threshold::T = 0.1) where {T <: AbstractFloat}

    @assert threshold > 0.0 "The given threshold value ($threshold) for `infer_bonds!` should be a value above 0.0!"

    dm      = collect(ProtoSyn.Calculators.full_distance_matrix(pose))
    atoms   = collect(eachatom(pose.graph))
    for (i, atom_i) in enumerate(atoms)
        for (j, atom_j) in enumerate(atoms)
            i == j && continue
            atom_j = atoms[j]
            atom_j in atom_i.bonds && continue
            putative_bond = "$(atom_i.symbol)$(atom_j.symbol)"

            if !(putative_bond in keys(ProtoSyn.Units.bond_lengths))
                continue
            end

            d = ProtoSyn.Units.bond_lengths[putative_bond]
            d += d * threshold
            if dm[i, j] < d && atom_i.container.container === atom_j.container.container
                ProtoSyn.bond(atom_i, atom_j)
            end
        end
    end

    return pose
end


"""
    identify_atom_by_bonding_pattern(container::AbstractContainer, pattern::Vector{String})

Returns one or more candidate [`Atom`](@ref) instances from the given
`AbstractContainer` `container` that match the provided `pattern` (a Vector of
[`Atom`](@ref) elements). This method follows the following hierarchical criteria:

# Examples
```
julia> ProtoSyn.identify_atom_by_bonding_pattern(pose.graph[1][1], ["H", "N", "C", "C"])
3-element Vector{Atom}:
 Atom{/2a3d:3900/A:1/MET:1/H1:2}
 Atom{/2a3d:3900/A:1/MET:1/H2:3}
 Atom{/2a3d:3900/A:1/MET:1/H3:4}

julia> ProtoSyn.identify_atom_by_bonding_pattern(pose.graph[1][1], ["C", "C", "C", "C", "H"])
Atom{/2a3d:3900/A:1/MET:1/C:18}
```

"""
function identify_atom_by_bonding_pattern(container::AbstractContainer, pattern::Vector{String})

    function seek_pattern(atom::Atom, inner_pattern::Vector{String}, previous::Opt{Atom})
        inner_candidates = [a for a in atom.bonds if (a.symbol === inner_pattern[1]) && (a !== previous)]
        length(inner_candidates) === 0 && return false
        deleteat!(inner_pattern, 1)
        length(inner_pattern) === 0 && return true
        return [seek_pattern(c, copy(inner_pattern), atom) for c in inner_candidates]
    end

    candidates = [a for a in eachatom(container) if a.symbol === pattern[1]]
    length(candidates) === 1 && return candidates[1]
    deleteat!(pattern, 1)
    length(pattern) === 0 && return candidates

    c = [seek_pattern(c, copy(pattern), nothing) for c in candidates]

    # Unravel results
    _c = Vector{Bool}([])
    for group in c
        for _ in 1:length(pattern)
            group = reduce(vcat, group)
        end
        push!(_c, any(group))
    end

    hits = candidates[_c]
    if length(hits) === 1
        return hits[1]
    else
        return hits
    end
end


"""
    sequence(container::ProtoSyn.AbstractContainer)::String
    sequence(pose::Pose)::String

Return the sequence of residues (in 1 letter code) of the given container/
[`Pose`](@ref)
as a String. Checks the `ProtoSyn.three_2_one` dictionary for name to 1 letter
code translation, uses '?' if no entry was found.

# Examples
```
julia> ProtoSyn.Peptides.sequence(pose)
"SESEAEFKQRLAAIKTRLQAL"
```
"""
function sequence(container::ProtoSyn.AbstractContainer)::String

    sequence = ""
    for residue in eachresidue(container)
        try
            sequence *= ProtoSyn.three_2_one[residue.name.content]
        catch KeyError
            sequence *= '?'
        end
    end

    return sequence
end

sequence(pose::Pose) = sequence(pose.graph)


# --- Visualize ----------------------------------------------------------------
# Note: This functions are incomplete, as they cannot work with ramifications.
# This may change in future iterations of ProtoSyn.

function visualize(io::IO, topology::ProtoSyn.Topology)
    _root = ProtoSyn.root(topology)
    for child in _root.container.children[1:(end - 1)]
        println(io, " ├── ", ProtoSyn.visualize(child))
    end
    child = _root.container.children[end]
    print(io, " └── ", ProtoSyn.visualize(child))
end

visualize(topology::Topology) = visualize(stdout, topology)
visualize(pose::Pose{Topology}) = visualize(stdout, pose.graph)

function visualize(residue::Residue)
    s = "|$(residue.id)"
    length(residue.children) == 0 && return s * "|"
    for child in [residue.children[1]]
        s *= visualize(child)
    end
    return s
end

visualize(segment::Segment) = visualize(ProtoSyn.origin(segment).container)

visualize(io::IO, frag::Fragment) = begin
    println(io, " ○ No Root")
    print(io, " └── ", ProtoSyn.visualize(frag.graph))
end

visualize(frag::Fragment) = visualize(stdout, frag)