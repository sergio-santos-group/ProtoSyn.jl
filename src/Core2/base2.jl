
@doc """
finds atom indices corresponding to the requested atom names, based on
the given residue. Names starting with '-'/'+' are assumed to belong to the
previous/next residue. If multiple residues are connected to the central (given)
residue, all tuple are generated.
"""
function cproduct(f::Function, residue::Residue, names::Tuple{Vararg{String}})
    st = Vector{Atom}(undef, length(names))
    _cproduct(f, residue, names, 0, st, residue)
end

function cproduct(f::Function,
    residue::Residue, names::Tuple{Vararg{String}}, n::Int, state::Vector{Atom})
    _cproduct(f, residue, names, n, state, residue)
end

function _cproduct(f::Function,
    residue::Residue, names::Tuple{Vararg{String}}, n::Int, state::Vector{Atom}, pres::Residue)
    # pres -> previous residue

    # apply function if the current state is complete
    if n==length(names)
        f(state)
    else
        name = names[n+1]
        if startswith(name, '+')
            for link in residue.links
                # is link.residue1 the source (central) residue? if so, the
                # next ('+') residue is residue2.
                # It is assumed that every link is directed from
                # residue1 (source) -> residue2 (target)
                if link.r1 === residue
                    (n>0) && startswith(names[n],'+') && (link.r2 !== pres) && continue
                    atom = get(link.r2, SubString(name,2), nothing)
                    (atom===nothing) && continue
                    state[n+1] = atom
                    _cproduct(f, residue, names, n+1, state, link.r2)
                end
            end
        elseif startswith(name, '-')
            for link in residue.links
                # is link.residue2 the target residue? if so, the
                # previous ('-') residue is residue1
                if link.r2 === residue
                    (n>0) && startswith(names[n],'-') && (link.r1 !== pres) && continue
                    atom = get(link.r1, SubString(name,2), nothing)
                    (atom===nothing) && continue
                    state[n+1] = atom
                    _cproduct(f, residue, names, n+1, state, link.r1)
                end
            end
        else
            atom = get(residue, name, nothing)
            if atom!==nothing
                state[n+1] = atom
                _cproduct(f, residue, names, n+1, state, residue)
            end
        end
    end

end

# function splitchains(mol::Molecule)
#     if !isvalid(mol)
#         error("splitchains requires an updated molecule")
#     end
#     ts = mol.traverse_state
#     reset(ts)
    
#     for atom in iterbyatom(mol)
#         if
#     end

# end

function Base.replace!(f::F, mol::Molecule, state::State, old_new::Pair{Residue,Residue}, linker::Function) where {F<:Function}
    
    idx = findfirst(r->r===old_new.first, mol.residues)

    # do nothing if the old residue does not exist
    if idx === nothing
        return mol
    end
    

    old = old_new.first
    new = deepcopy(old_new.second)
    
    # update internals
    new.id = old.id
    new.molecule = mol
    old.molecule = nothing
    # mol.residues[idx] = new
    
    # reuse links
    while !isempty(old.links)
        l = pop!(old.links)
        #push!(new.links, link)
        #setproperty!(link, link.r1 === old ? :r1 : :r2, new)
        delete!(mol, l)
        if l.r1 === old
            delete!(l.r2, l)
            link(linker, new, l.r2)
        else
            delete!(l.r1, l)
            link(linker, l.r1, new)
        end
    end
    
    newlen = length(new)
    oldlen = length(old)

    f(
        old, view(state.coords,:,old.offset+1:old.offset+oldlen),
        new, new.coords,
        state
    )

    #copyto!(old.coords, 1, state.coords, 3*old.offset+1, 3*oldlen)

    # shifts coordinates left/right and allocate space if necessary
    if newlen != oldlen
        natoms = size(mol)
        newsize = natoms - oldlen + newlen
        if newsize > state.capacity
            resize!(state, newsize)
        end
        # ony shift coordinates if this is not the last residue
        if idx < length(mol)
            src = copy(state.coords)
            copyto!(state.coords,                 # destiny
                    3*(old.offset+newlen)+1,      # destiny offset
                    #state.coords,                 # source
                    src,
                    3*(old.offset+oldlen)+1,      # source offset
                    3*(natoms-old.offset-oldlen)  # number of elements to copy
                )
        end
        state.size = newsize
    end
    mol.residues[idx] = new
    # copy new coordinates
    copyto!(state.coords, 3*old.offset+1, new.coords, 1, 3*newlen)
    
    update!(mol)
    println("BEFORE EVERYTHING ", state.size)
    foreach(println, eachcol(state.coords))
    mol
end




export set!

function set!(state::State, mol::Molecule, rng::UnitRange{Int}, atnames::NTuple{4, String}, θ::Float64)
    
    # Allocate-once variables
    container = Vector{Atom}(undef, 4)

    for i in rng
        cproduct(mol.residues[i], atnames, 0, container) do atoms
            dihd = Dihedral(atoms[1], atoms[2], atoms[3], atoms[4], nothing)
            #println(dihd, " ", measure(dihd, state))
            rotate!(state, mol, dihd, θ-measure(dihd, state.coords))
        end
    end
    state
end

function set!(state::State, mol::Molecule, dihedral::Dihedral, θ::Float64)
    rotate!(state, mol, dihedral, θ-measure(dihedral, state.coords))
end

# function set!(mol::Molecule, state::State) = 
#     mol.coords = view
#     edn