
const BondedList = Dict{Union{Int,Symbol}, Vector{Tuple{Vararg{Int}}}}
const ExclusionList = Dict{Int, Vector{Int}}


function genbonded(mol::Molecule, maxdepth::Int=3)
    if mol.bonds !== nothing
        return genbonded(mol.bonds, maxdepth)
    end
    nothing
end


function genbonded(graph::ConnectGraph, maxdepth::Int=3)
    
    bonded = BondedList()

    visited = falses(length(graph))
    
    function traverse(pivot::Int, path::Vector{Int}, depth::Int=0)
        # do nothing if one is beyond the maximum
        # allowed depth
        if depth > maxdepth
            return
        end
        
        # add this item to the path
        push!(path, pivot)

        # save the path to the adequate container:
        #  by now, the depth correponds to the number
        #  of bonds away from the initial atom 
        if (depth > 0) && (path[1] < path[end])
            push!(get!(bonded, depth, []), Tuple(path))
        end
        
        # if we have reached the maximum depth, then return
        #  and prevent another cycle of traversal that will not
        #  be useful.
        if depth == maxdepth
            return
        end

        # traverse the graph 
        visited[pivot] = true
        for i in graph[pivot]
            if !visited[i]
                traverse(i, copy(path), depth+1)
            end
        end
        visited[pivot] = false
        
    end

    for i = 1:length(graph)
        traverse(i, Int[])
    end

    # if the 1 to 3 keys exist, create named alias
    haskey(bonded, 1) && (bonded[:bonds ] = bonded[1])
    haskey(bonded, 2) && (bonded[:angles] = bonded[2])
    haskey(bonded, 3) && (bonded[:proper] = bonded[3])
    haskey(bonded, maxdepth) && (bonded[:pairs] = bonded[maxdepth])

    return bonded
end





function genexcluded(mol::Molecule, maxdepth::Int=3)
    if mol.bonds !== nothing
        return genexcluded(mol.bonds, maxdepth)
    end
    nothing
end


function genexcluded(graph::ConnectGraph, maxdepth::Int=3)
    
    excluded = ExclusionList()
    visited = falses(length(graph))
    
    function traverse(pivot::Int, cur::Int, path::Vector{Int}, depth::Int=0)
        # do nothing if already beyond the max depth
        (depth > maxdepth) && return
        
        # add this index tp the path only if depth > 0.
        # Otherwise, the first entry of the exclusion list
        # for each atom is the atom itself, and we want the
        # exclusion list to contain only indices greater
        # than the key itsefl
        (depth > 0) && (cur > pivot) && push!(path, cur)
        
        # prevent the next function call overhead
        (depth == maxdepth) && return

        # mark this atom as visited , traverse the graph
        # and unmark the atom back
        visited[cur] = true
        for i in graph[cur]
            !visited[i] && traverse(pivot, i, path, depth+1)
        end
        visited[cur] = false
    end
    

    for k = 1:length(graph)
        # get the list for atom k
        excl = get!(excluded, k, [])
        
        # build the list and sort it
        traverse(k, k, excl)
        sort!(excl)

        # add a terminal index to flag end of list (this
        # index does not exist and ensures consistency of the list)
        push!(excl, -1)
    end

    excluded
end







@inline function genkey(types::T...) where {T<:AbstractString}
    join(types, ":")
end



# function load(ffield::String)
#     filename = joinpath(resource_dir, ffield, "forcefield.yml")
#     println("loading forcefield $filename")
#     open(filename) do io
#         load(io)
#     end
# end



# function load(io::IO)
#     raw = YAML.load(io)

#     ffield = Forcefield()

#     # iterate over forcefield components
#     #  (bonds, angles, proper dihedrals, improper ...)
#     for component in keys(raw)
        
#         # Instantiate an appropiate container for this type
#         #  this is a dictionary with keys given by
#         #  the atomtypes comprising the component
#         container = Dict{String, IForcefieldType}()
#         ffield[Symbol(component)] = container

        
#         for (typename, entries) in raw[component]
            
#             # get the adequate type for this component.
#             # These will be subtypes of IForcefieldType
#             # such as HarmonicBondType, HarmonicAngleType, ...
#             type = getfield(Forcefields, Symbol(typename))
            
#             # determine how many entry fields define the key
#             #    - bonds require 2 fields;
#             #    - angles require 3 fields,
#             #    - (im)proper require 4 fields
#             n = length(entries[1]) - fieldcount(type)
#             if type <: IAppendableForcefieldType
#                 n += 1
#             end
            
#             # iterate over all entries for the current
#             # component type. These will correspond to
#             # different combinations of atom types.
#             for entry in entries
                
#                 # convert all arguments to the required
#                 # types of the target type
#                 args = [
#                     convert(eltype(t),v)
#                     for (t,v) in zip(fieldtypes(type), entry[n+1:end])
#                 ]
                
#                 # generate a key and get the item from the container
#                 # (if it exists). If it is a subtype of 
#                 # IAppendableForcefieldType, then multiple definitions
#                 # in the raw data are appended to that instance. Otherwise,
#                 # it is replaced!
#                 key = genkey(entry[1:n]...)
#                 item = get(container, key, nothing)
#                 if isa(item, IAppendableForcefieldType)
#                     push!(item, args...)
#                 else
#                     container[key] = type(args...)
#                 end
                
#             end
#         end
#     end

#     return ffield
# end




function lookup(component::Symbol, comptypes::Dict{String, IForcefieldType}, types::String...)
    
    # direct order
    key = genkey(types...)
    haskey(comptypes, key) && return comptypes[key]
    
    if component === :improper
        key = genkey("X", types[2:end]...)
        haskey(comptypes, key) && return comptypes[key]
        
        key = genkey("X", "X", types[3:end]...)
        haskey(comptypes, key) && return comptypes[key]
    elseif component === :proper
        key = genkey(reverse(types)...)
        haskey(comptypes, key) && return comptypes[key]
    
        key = genkey("X", types[2], types[3], "X")
        haskey(comptypes, key) && return comptypes[key]

        key = genkey("X", types[3], types[2], "X")
        haskey(comptypes, key) && return comptypes[key]
    else
        key = genkey(reverse(types)...)
        haskey(comptypes, key) && return comptypes[key]
    end

    nothing
end



function loadmap(mapname)
    filename = joinpath(resource_dir, mapname)
    open(filename) do io
        loadmap(io)
    end
end

loadmap(io::IO) = YAML.load(io)

# function loadmap(io::IO)
#     raw = YAML.load(io)
    
#     d = Dict()

#     # reverse first and second levels
#     for (resname,components) in raw
#         for compname in keys(components)
#             key = Symbol(compname)
#             c = get!(d, key, Dict())
#             c[resname] = components[compname]
#         end
#     end

#     d
# end



function assign(mol::Molecule, ff::Forcefield, resinfo::Dict)
    
    if !mol.coherent
        error("molecule is not coherent")
    end
    
    # the topology itsef. It is currently empty (all container
    # attributes are empty)
    top = Topology()
    top2 = Topology2()

    # temporary array for storing additional interaction indices
    # for use with the cporduct function (prevents multiple memory
    # allocations ahead)
    add_indices = zeros(Int,10)


    # copy the connectivity graph to be used in the latter
    # appending of additional ff components 
    graph = ConnectGraph(
        k=>v[:]
        for (k,v) in mol.bonds
    )
    

    # Add additional bonds. These additional bonds are
    # used for exclusion generation
    for residue in mol.residues
        
        resmap = resinfo[residue.source.name]
        addbonds = get(resmap, "bonds", nothing)
        (addbonds===nothing) && continue

        for bond in addbonds
            ProtoSyn.cproduct(residue, bond, 0, add_indices) do ij
                !in(ij[2], graph[ij[1]]) && push!(graph[ij[1]], ij[2])
                !in(ij[1], graph[ij[2]]) && push!(graph[ij[2]], ij[1])
            end
        end
    end


    # generate bonded lists
    blist = genbonded(graph, ff.exclusion_depth)
    excluded = genexcluded(graph, ff.exclusion_depth)
    

    # assign atom types using the resinfo map
    atoms = map(ProtoSyn.iterbyatom(mol)) do at
        qt = resinfo[at.parent.name]["atoms"][at.name]
        q = qt["charge"]
        t = qt["type"]
        Atom(q, ff.components[:atoms][t])
    end
    top.atoms = AtomContainer(atoms)
    push!(top2, atoms...)
    
    # add exclusions
    for (k,excl) in excluded
        atoms[k].exclusions = excl
    end

    # deal with pairs (ex. 1-4 pairs for scaled nb interactions)
    # if ff.genpairs
    #     top.pairs = [
    #         (l[1], l[end])
    #         for l in blist[:pairs]
    #     ]
    # end
    println("do not forget about pairs")

    # add additional forcefield components (as defined by
    # the residue map info)
    for residue in mol.residues
        resmap = resinfo[residue.source.name]

        # iterate over additional ff component names/values
        for (compname, addcomps) in resmap
            # skip "atoms" (dealt with separately) and
            # "bonds" (already taken into account)
            in(compname, ("atoms", "bonds")) && continue

            blist_entry = get!(blist, Symbol(compname), [])

            # number of indices making up this component
            ilen = length(addcomps[1])

            for addcomp in addcomps
                ProtoSyn.cproduct(residue, addcomp, 0, add_indices) do indices
                    indices = Tuple(view(indices,1:ilen))
                    if indices[1] > indices[end]
                        indices = reverse(indices)
                    end
                    !in(indices, blist_entry) && push!(blist_entry, indices)
                end
            end
        end
    end


    # finally generate full Topology
    cache = Dict()
    
    for (compsymb, comptypes) in ff.components
        
        if (compsymb===:atoms) || !haskey(blist, compsymb)
            continue
        end

        for indices in blist[compsymb]
            
            # get atom types for these indices and look fot the
            # adequate component type in the forcefield. If not
            # found (nothing returned), skip to the next tuple of indices. 
            atomtypes = map(i -> atoms[i].type.name, indices)
            comptype = lookup(compsymb, comptypes, atomtypes...)
            if comptype === nothing
                println("unable to find ff entry \"$(compsymb)\" for types $(atomtypes) $(indices)")
                continue
            end

            typeofcomptype = typeof(comptype)
            if !haskey(cache, typeofcomptype)
                
                # get name of the required object type:
                #  ex: HarmonicBondType => HarmonicBond
                typename = string(nameof(typeofcomptype))
                typename = replace(typename, "Type"=>"")
                
                # # get the required object type:
                # #  ex: HarmonicBond (string) => HarmonicBond (DataType)
                # type = getfield(Forcefields, Symbol(typename))
                    
                # # get the name of the adequate container in the topology
                # #  ex: HarmonicBond (DataType) => harmonicBonds
                # propname = lowercasefirst(typename)*"s"

                # # initialize container (by default is nothing)
                # container_type = getfield(Forcefields, Symbol(typename*"Container"))
                # setproperty!(top, Symbol(propname), container_type())
                
                # # save info to cache
                # cache[typeof(comptype)] = (type, getproperty(top, Symbol(propname)))

                # get the required object type and cache it
                #  ex: HarmonicBond (string) => HarmonicBond (DataType)
                cache[typeofcomptype] = getfield(Forcefields, Symbol(typename))
            end

            # get target type and topology container from cache
            # type, container = cache[typeof(comptype)]
            type = cache[typeofcomptype]
            
            #instantiate a new object of type <type> and add
            # it to the container
            obj = type(indices..., comptype)
            # push!(container, obj)
            push!(top2, obj)
        end

    end

    top,top2
end


function aggregate(comps...)
    
    function(state::State, do_forces::Bool)
        if do_forces
            fill!(state.forces, 0.0)
        end
        e::Float64 = 0.0
        for comp in comps
            e += ProtoSyn.eval!(state, comp, do_forces)
        end
        state.energy.total = e
    end
end


