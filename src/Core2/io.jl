const PDB = Val{1};
const XYZ = Val{2};


function save(fname::String, mol::Molecule, state::State, mode="w")
    open(fname, mode) do fout
        write(fout, mol, state)
    end
end


function Base.write(io::IO, mol::Molecule, st::State, ::Type{PDB})
    println(io, "MODEL")
    for (rid,lr) in enumerate(mol.residues)
        offset = mol.offset + lr.offset
        for at in lr.source.atoms
            i = at.id + offset
            s = @sprintf("ATOM %6d %4s %-4sA %3d    %8.3f%8.3f%8.3f%22s%2s",
                i, at.name,
                at.parent.name, rid,
                st.coords[i,1], st.coords[i,2], st.coords[i,3],
                # 10.0*st.coords[i,1], 10.0*st.coords[i,2], 10.0*st.coords[i,3],
                lr.source.name,
                at.symbol)
                println(io, s)
        end
    end

    # for key in sort(collect(keys(mol.bonds)))
    #     s = map(x->@sprintf("%5d", x+mol.offset), mol.bonds[key])
    #     println(io, @sprintf("CONECT%5d", key+mol.offset), s...)
    # end
    println(io, "ENDMDL")
end

Base.write(io::IO, mol::Molecule, state::State) = Base.write(io, mol, state, PDB)

function Base.write(io::IO, st::State)
    println(io, "MODEL")
    x = st.coords
    for i = 1:st.size
        println(io, @sprintf("ATOM %6d %4s UNK A %3d    %8.3f%8.3f%8.3f",
            # i, "X", i, 10.0*x[i,1], 10.0*x[i,2], 10.0*x[i,3])
            i, "X", i, x[i,1], x[i,2], x[i,3])
        )
    end
    println(io, "ENDMDL")
end



function Base.write(io::IO, mol::Molecule, st::State, ::Type{XYZ})
    println(io, "$(mol.size)\ntitle")
    for lr in iterbyresidue(mol)
        offset = mol.offset + lr.offset
        for at in iterbyatom(lr)
            i = at.id + offset
            s = @sprintf("%-4s    %8.3f%8.3f%8.3f",
                at.symbol,
                st.coords[i,1], st.coords[i,2], st.coords[i,3])
                # 10.0*st.coords[i,1], 10.0*st.coords[i,2], 10.0*st.coords[i,3])
                println(io, s)
        end
    end
end


# function loadresidues(filename::String)
    
#     residues = Dict{String, Residue}()
    
#     frag_regex = r"^REMARK\s+400\s+NAME:\s+(?<name>\w+)"

#     open(filename, "r") do fin
#         let residue::Residue
#         for line in eachline(fin)
            
#             match_res = match(frag_regex, line)
#             if match_res != nothing
#                 residue = Residue(name=match_res[:name])
#                 residues[residue.name] = residue
            
#             elseif startswith(line, "ATOM")

#                 # atom = Atom{Residue}(
#                 atom = Atom(
#                     id = parse(Int, strip(line[7:11])),
#                     name = string(strip(line[13:16])),
#                     symbol = length(line)>77 ? string(strip(line[77:78])) : "?",
#                     x = parse(Float64, line[31:38]),
#                     y = parse(Float64, line[39:46]),
#                     z = parse(Float64, line[47:54])
#                 )
#                 push!(residue, atom)
                
#             elseif startswith(line, "CONECT")
#                 (residue.bonds == nothing) && (residue.bonds = ConnectGraph())
                
#                 indices = map(s->parse(Int, s), [line[n:n+4] for n=7:5:length(line)])
#                 for id2 in indices[2:end]
#                     id1 = indices[1]
#                     if id1 > id2
#                         id1, id2 = id2, id1
#                     end
#                     bonds = get!(residue.bonds, id1, [])
#                     if id2 ∉ bonds
#                         push!(bonds, id2)
#                         push!(get!(residue.bonds, id2, []), id1)
#                     end
#                 end

#             end

#         end # end for line in
#         end # end let
#     end # end open

#     residues

# end






# function readpdb(filename::String, fragments::Dict{String, Fragment})
#     open(filename) do fin
#         readpdb(fin, fragments)
#     end
# end

# function readpdb(io::IO, fragments::Dict{String, Fragment})::Tuple{Vector{Molecule}, State}
    
#     # === HOUSEKEEPING VARIABLES ===
#     offset = 0
#     atoms_left = 0
#     atom_counter = 0
#     create_mol = true
#     create_frag = true
    
#     lfrag2mol  = IdDict{LinkedFragment, Molecule}()
#     atmid2frag = Dict{Int,LinkedFragment}()
#     pair2link  = IdDict{UInt, Link}()
    
#     # === OUTPUT VARIABLES ===
#     # determine the number of atoms in the PDB file and
#     # instantiate a new state to accomodate the currrent coordinates
#     state = State(
#         count(l->startswith(l, "ATOM"), eachline(io))
#     )
#     molecules = Molecule[]
    

#     # do some work ...
#     let lfrag::LinkedFragment,
#         molecule::Molecule
        
#         seekstart(io)
#         for line in eachline(io)
            
#             if startswith(line, "ATOM")

#                 # should a new molecule be instantiated?
#                 if create_mol
#                     molecule = Molecule(offset=offset)
#                     push!(molecules, molecule)
#                     create_mol = false
#                     # this molecule is empty, hence a new
#                     # l-frag should be instantiated.
#                     create_frag = true
#                 end

#                 # should a new l-fragment be instantiated?
#                 if create_frag
                    
#                     # Check if this a non-standard fragment:
#                     #  it is considered a non-standard fragment
#                     #  if columns 73:76 exist and are non empty.
#                     #  otherwise, the fragment name is given by the residue
#                     #  name in columns 18:20.
#                     if (length(line) >= 76) && !isempty(strip(line[73:76]))
#                         fragname = string(strip(line[73:76]))
#                     else
#                         fragname = string(strip(line[18:20]))
#                     end

#                     # instantiate a linked-fragment pointing
#                     # to the fragment template identified above. If
#                     # it does not exist, then raise an error!
#                     if !haskey(fragments, fragname)
#                         error("Unknown fragment name $(fragname)")
#                     end

#                     lfrag = LinkedFragment(
#                         source=fragments[fragname]
#                     )

#                     # push the linked-fragment to the molecule.
#                     # this process already assigns id and offset
#                     # according to the position of the l-fragment in
#                     # the molecule. However, and because of the latter
#                     # CONECT records, one should use global offsets so
#                     # that connectivities can be latter identified correctly.
#                     push!(molecule, lfrag)
#                     lfrag.offset = offset

#                     # add this l-fragment to the frag-to-molecule map
#                     lfrag2mol[lfrag] = molecule

#                     # update atom counter (used to identify when a single
#                     # fragment has been fully read), offsets and flags
#                     atoms_left = length(lfrag)
#                     offset += atoms_left
#                     create_frag = false
#                 end

                
#                 # add this atom id to the adequate map
#                 # so that one can latter retrieved the l-fragment
#                 # to which this atom belongs. Meanwhile, check if the atom names
#                 # in the current line and correponding fragment's atomname match 
#                 atmid = parse(Int, strip(line[7:11]))
#                 expname = lfrag.source.atoms[end-(atoms_left-1)].name
#                 curname = strip(line[13:16])
#                 if expname != curname
#                     error("name mismatch between atom $(atmid) and fragment ($curname vs $expname)")
#                 end
#                 atmid2frag[atmid] = lfrag
                
#                 # parse coordinates into state:
#                 #  one should convert PDB coordinates (Angstrom)
#                 #  to nanometers (hence the 0.1 factor). 
#                 atom_counter += 1
#                 state.coords[atom_counter, 1] = 0.1*parse(Float64, line[31:38])
#                 state.coords[atom_counter, 2] = 0.1*parse(Float64, line[39:46])
#                 state.coords[atom_counter, 3] = 0.1*parse(Float64, line[47:54])

#                 # decrement atom counter: when it reaches zero, then all
#                 # atoms for the current l-fragment have been read and one
#                 # should indicate that a new l-fragment should be 
#                 # instantiated, even if the next one (if existent) is of
#                 # similar nature.
#                 atoms_left -= 1
#                 if atoms_left == 0
#                     create_frag = true
#                 end
                
#             elseif startswith(line, "CONECT")

#                 # this is an array having the following format:
#                 #   [pivot, at1, at1, ..., atN]
#                 # in which at1:atN are all connected to pivot.
#                 line = strip(line) # remove trailing spaces
#                 ids = map(s->parse(Int, s), [line[n:n+4] for n=7:5:length(line)])
                
#                 for id2 in ids[2:end]
                    
#                     # id1 is the pivot (above) and id2 is 
#                     # connected to id1. To ensure that only a single
#                     # link is ever created, one requires that id1 < id2
#                     id1 = ids[1]
#                     if id1 > id2
#                         id1, id2 = id2, id1
#                     end

#                     # so far, id1 and id2 are global indices. To find
#                     # the l-fragments they belong to, one uses the
#                     # atmid2frag map. Once the l-fragment is found, then
#                     # one can convert from global to local (local to the
#                     # Fragment) by subtracting the offset.
#                     lf1 = atmid2frag[id1]
#                     lf2 = atmid2frag[id2]
#                     id1 -= lf1.offset
#                     id2 -= lf2.offset
                    
#                     # Intermolecular links are not allowed. If any
#                     # 2 molecules were to be connected, they would
#                     # form a single molecule, and would no longer be
#                     # 2 separate molecules!!
#                     if lfrag2mol[lf1] !== lfrag2mol[lf2]
#                         error("Intermolecular links are not allowed")
#                     end

#                     # if this is a intramoleculer link, then check if
#                     # it is already defined in the Fragment internal
#                     # connectivity table. If not, this is a custom bond
#                     # that should be defined by a Link.
#                     if (lf1===lf2) && (id1 in lf2.source.bonds[id2])
#                         continue
#                     end

#                     # to ensure a single link is created for any pair
#                     # of linked l-fragments, one requires that the id
#                     # of l-fragment1 is less than that of l-fragment2.
#                     if lf1.id > lf2.id
#                         continue
#                     end

#                     println("$(id1+lf1.offset).$(lf1.id) <-> $(id2+lf2.offset).$(lf2.id)")
                    
#                     # generate the pair key for l-fragments 1 and 2
#                     # this is a simple way to generate a unique key
#                     # for the pair (others exist!).
#                     pairkey = hash(lf1, hash(lf2))
                    
#                     if !haskey(pair2link, pairkey)
#                         # if this pair does not exist, then create a new link,
#                         # add it to the map, molecule, and owner l-fragments,
#                         link = Link(lf1, lf2, ConnectGraph())
#                         pair2link[pairkey] = link
#                         push!(lfrag2mol[lf1], link)
#                         push!(lf1, link)
#                         push!(lf2, link)
#                     else
#                         # otherwise, simply retrive the link from the map
#                         println("fragment $(lf1.id) already linked to $(lf2.id)")
#                         link = pair2link[pairkey]
#                     end

#                     # given the link, check if id1 is already connected
#                     # to id2. If not, connect them.
#                     bonds = get!(link.bonds, id1, [])
#                     if id2 ∉ bonds
#                         push!(bonds, id2)
#                     end

#                 end # for id2 in indices

#             elseif startswith(line, "TER")
#                 # molecules are delimited by TER tags. Whenever such
#                 # a tag is encountered, one should indicate that a new
#                 # molecule ought to be created.
#                 create_mol = true
#             end

#         end # end for line
#     end # end length

#     molecules, state
# end
