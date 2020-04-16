using Printf

const PDB = Val{1}
const XYZ = Val{2}

Base.read(filename::AbstractString, ::Type{PDB}) = begin
    open(filename) do fin
        read(fin, PDB)
    end
end


Base.read(io::IO, ::Type{PDB}) = begin
    
    seekstart(io)
    natoms = mapreduce(l->startswith(l, "ATOM")||startswith(l, "HETATM"), +, eachline(io); init=0)
    
    mol = Molecule()    # returned object
    res = Residue()     # orphan residue
    state = State{Float64}(natoms)

    serial2lid = Dict{Int,Int}()     # atom serial number to atom ID
    serial2res  = Dict{Int,Residue}() # atom serial number to residue
    bonds = Set{Pair{Int,Int}}()
    lid = gid = 0    # local (intra-residue) and global atom ID
    
    seekstart(io)
    for line in eachline(io)
        
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            
            resname = string(strip(line[18:20]))
            resid = parse(Int, line[23:26])
            if res.id != resid || res.name != resname
                res = Residue()
                push!(mol, res)
                res.name = resname
                res.id = resid
                lid = 0
            end

            atom = Atom()
            #lid = parse(Int, line[7:11])
            #atom.id = parse(Int, line[7:11])
            
            # atom.metadata = Dict(
            #     #:x => parse(Float64, line[31:38]),
            #     #:y => parse(Float64, line[39:46]),
            #     #:z => parse(Float64, line[47:54]),
            #     :symbol => length(line)>77 ? string(strip(line[77:78])) : "?"
            # )
            atom.symbol = length(line)>77 ? string(strip(line[77:78])) : "?"
            atom.name = string(strip(line[13:16]))
            atom.id = (lid += 1)
            push!(res,atom)

            gid += 1
            state.coords[1,gid] = parse(Float64, line[31:38])
            state.coords[2,gid] = parse(Float64, line[39:46])
            state.coords[3,gid] = parse(Float64, line[47:54])

            serial = parse(Int, line[7:11])
            serial2lid[serial] = atom.id
            serial2res[serial] = res
        elseif startswith(line, "CONECT")
            idxs = map(s->parse(Int, s), [line[n:n+4] for n=7:5:length(line)])
            for i in idxs[2:end]
                idxs[1] < i ? push!(bonds, idxs[1]=>i) : push!(bonds, i=>idxs[1])
            end
        end
    end
    
    if !isempty(bonds)
        links = Dict{Tuple{Residue,Residue},Link}()
        for bond in bonds
            r1 = serial2res[bond.first]
            r2 = serial2res[bond.second]
            id1 = serial2lid[bond.first]
            id2 = serial2lid[bond.second]

            # intra-residue bonds
            if r1 === r2
                push!(get!(r1.bonds, id1, []), id2)
                push!(get!(r2.bonds, id2, []), id1)
            
            # inter-residue bonds
            else
                if r1.id > r2.id
                    r1, r2 = r2, r1
                    id1, id2 = id2, id1
                end

                key = (r1,r2)
                if !haskey(links, key)
                    links[key] = link((a,b)->Link(a, b, id1=>id2), r1, r2)
                else
                    push!(links[key], id1=>id2)
                end
            end
        end
    end
    # identify the root residue as the one with the
    # lowest ID
    if length(mol) > 0
        i = argmin(map(r->r.id, mol.residues))
        mol.root = mol.residues[i]
    end
    update!(mol), state
end


Base.write(io::IO, mol::Molecule, state::State, ::Type{PDB}) = begin
    println(io, "MODEL")
    xyz = state.coords

    for atom in eachatom(mol)
        id = atom.id + atom.residue.offset
        s = @sprintf("ATOM %6d %4s %-4sA %3d    %8.3f%8.3f%8.3f%24s",
            id, atom.name,
            atom.residue.name, atom.residue.id,
            xyz[1,id],
            xyz[2,id],
            xyz[3,id],
            # atom.metadata[:x],
            # atom.metadata[:y],
            # atom.metadata[:z],
            atom.symbol)
        println(io, s)
    end
    for atom in eachatom(mol)
        id = atom.id + atom.residue.offset
        # @sprintf("CONECT%5d%s", id,
        # jo mol.bonds[id]
        # )
        # join()
        print(io, @sprintf("CONECT%5d",id))
        foreach(i->print(io, @sprintf("%5d",i)), mol.bonds[id])
        println(io,"")
    end
    println(io, "ENDMDL")
end

# ATOM     32  N  AARG A  -3      11.281  86.699  94.383  0.50 35.88           N  
# ATOM     33  N  BARG A  -3      11.296  86.721  94.521  0.50 35.60           N
# ATOM     34  CA AARG A  -3      12.353  85.696  94.456  0.50 36.67           C
# ATOM     35  CA BARG A  -3      12.333  85.862  95.041  0.50 36.42           C
# ATOM     11 HD21 ASN     1       4.408   6.734   2.315  0.00  0.00           H
