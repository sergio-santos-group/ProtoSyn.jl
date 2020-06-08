# using ...ProtoSyn
# include("types.jl")
# include("potentials.jl")

export genbonded, genff


function genbonded_dfs(f::Function, pivot::Atom, path::Vector{Atom}, maxdepth::Int, depth::Int=0)
    push!(path, pivot)
    pivot.visited = true
    f(depth, path)
    if depth < maxdepth
        for atom in pivot.bonds
            !atom.visited && genbonded_dfs(f, atom, path, maxdepth, depth+1)
        end
    end
    pivot.visited = false
    pop!(path)
end


function genbonded(graph::ProtoSyn.AbstractContainer, depth::Int)
    
    blist = BondedList()
    byatom = eachatom(graph)
    path = Vector{ProtoSyn.Atom}()
    
    for atom in byatom
        atom.visited = false
    end

    for i = 1:depth
        blist[i] = Vector{NTuple{i,Atom}}()
    end
    for atom in byatom
        genbonded_dfs(atom, path, depth) do d,pth
            if pth[1].index < pth[end].index
               push!(blist[d], Tuple(pth))
            end
        end
    end

    haskey(blist, 1) && (blist[:bonds ] = blist[1])
    haskey(blist, 2) && (blist[:angles] = blist[2])
    haskey(blist, 3) && (blist[:proper] = blist[3])
    haskey(blist, depth) && (blist[:pairs] = blist[depth])

    blist
end


function genff(graph::ProtoSyn.AbstractContainer, ffparams::ForcefieldParameters, resmaps::Dict)

    T = eltype(ffparams)

    blist = genbonded(graph, ffparams.exclusion_depth)
    
    byres = eachresidue(graph)
    
    # verify if the residue map contains the required entries
    resnames = Set(map(r->r.name, byres))
    for name in resnames
        !haskey(resmaps, name) && error("unable to find map for residue $name")
    end

    # FILL ADDITIONAL FORCEFIELD COMPONENTS
    for residue in eachresidue(graph)
        resmap = resmaps[residue.name]
        for (compname, addcomps) in resmap
            compname=="atoms" && continue

            blistkey = Symbol(compname)
            blistitems = get!(blist, blistkey, [])
            for addcomp in addcomps
                item = cprod(residue, addcomp...)
                item !== nothing && push!(blistitems, Tuple(item))
            end
        end
    end
    
    ff = Forcefield()

    natoms = ProtoSyn.count_atoms(graph)
    atomtypes = Vector{String}(undef, natoms)
    nbatoms = []

    # Use the residue maps to assign atom types
    for residue in eachresidue(graph)

        # info map for the atoms in this residue
        atmap = resmaps[residue.name]["atoms"]
        for atom in eachatom(residue)
            
            if !haskey(atmap, atom.name)
                error("unable to find atomic details for atom $(atom.name)")
            end

            # get information for this atom. It should
            # contain the charge and atom type; save the atomtype
            # for latter use
            atinfo = atmap[atom.name]
            atomtype = atinfo["type"]
            atomchrg = atinfo["charge"]
            atomtypes[atom.index] = atomtype

            # look for the adequate potential for this atomtype
            typedef = ffparams[:atoms, atomtype]
            if typedef === nothing
                @warn "unable to find ff entry 'atoms' for type $atomtype with index $(atom.index)"
                continue
            end
            
            # instantiate a concrete object for this potential
            # and add it to the forcefield
            obj = concrete(typedef, (atom.index,), atomchrg)
            push!(ff, obj)
            push!(nbatoms, obj)
        end
    end

    
    # go over all force-field components in the forcefield definition
    for (compkey, comptypes) in ffparams.components
        
        # skip this component if it relates to :atoms or
        # the blist does not contain an entry for this component
        if (compkey===:atoms) || !haskey(blist, compkey)
            continue
        end

        for atmtuples in blist[compkey]
            # retrieve atom indices and types
            indxs = map(at->at.index, atmtuples)
            types = map(i->atomtypes[i], indxs)

            # look for the forcefield potential type definition
            typedef = ffparams[compkey, types...]
            if typedef === nothing
                @warn "unable to find ff entry '$compkey' for types $types with indices $indxs"
                continue
            end
            # instantiate a concrete object for this potential
            # and add it to the forcefield
            obj = concrete(typedef, indxs)
            push!(ff, obj)
        end
    end
    
    genexcluded!(ff, graph, ffparams.exclusion_depth)

    if ffparams.genpairs
        for item in blist[:pairs]
            i = item[1].index
            j = item[end].index
            obj = combine(nbatoms[i], nbatoms[j])
            obj = scale(obj; q=ffparams.fudgeQQ, lj=ffparams.fudgeLJ)
            push!(ff, obj)
        end
    end

    @warn "add pairs"

    ff
end


function genexcluded!(ff::Forcefield, graph::ProtoSyn.AbstractContainer, maxdepth::Int)
    lst = Int[]
    exclusions = ff.exclusions
    for atom in eachatom(graph)
        empty!(lst)
        index = atom.index
        exclusions_dfs(atom, maxdepth) do at
            at.index > index && push!(lst, at.index)
        end
        exclusions[index] = sort(lst)
    end
end

function exclusions_dfs(f::Function, pivot::Atom, maxdepth::Int, depth::Int=0)
    f(pivot)
    pivot.visited = true
    if depth < maxdepth
        for atom in pivot.bonds
            !atom.visited && exclusions_dfs(f, atom, maxdepth, depth+1)
        end
    end
    pivot.visited = false
end


"""
!!!
    NOTE: requires attention
"""
function cprod(residue::Residue, names::String...)
    stack = []
    for name in names
        subname = SubString(name,2)
        if startswith(name, '-') && hasparent(residue) && haskey(residue.parent, subname)
            push!(stack, residue.parent[subname])
        elseif startswith(name, '+') && haschildren(residue) && haskey(residue.children[1], subname)
            push!(stack, residue.children[1][subname])
        elseif haskey(residue, name)
            push!(stack, residue[name])
        else
            return
        end
    end
    stack
end
