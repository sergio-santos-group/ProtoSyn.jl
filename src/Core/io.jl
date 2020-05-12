using Printf: @sprintf


const PDB = Val{1}


Base.read(::Type{T}, filename::AbstractString, ::Type{PDB}) where {T<:AbstractFloat} = begin
    pose = open(filename) do fin
        read(T, fin, PDB)
    end
    pose.graph.name = basename(filename)
    pose
end
Base.read(filename::AbstractString, ::Type{PDB}) = read(Float64, filename, PDB)

Base.read(::Type{T}, io::IO, ::Type{PDB}) where {T<:AbstractFloat} = begin
    
    top = Topology("?", -1)
    seg =  Segment("?", -1)     # orphan segment
    res =  Residue("?", -1)     # orphan residue
    
    seekstart(io)
    natoms = mapreduce(l->startswith(l, "ATOM")||startswith(l, "HETATM"), +, eachline(io); init=0)
    
    id2atom = Dict{Int,Atom}()
    
    state = State{T}(natoms)
    
    segid = atmindex = 1
    
    seekstart(io)
    for line in eachline(io)
        
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            
            resname = string(strip(line[18:20]))
            segname = string(strip(string(line[22])))
            resid = parse(Int, line[23:26])

            if seg.name != segname
                seg = Segment!(top, segname, segid)
                segid += 1
            end

            if res.id != resid || res.name != resname
                res = Residue!(seg, resname, resid)
            end

            atsymbol = length(line)>77 ? string(strip(line[77:78])) : "?"
            atname = string(strip(line[13:16]))
            atid = parse(Int, line[7:11])

            atom = Atom!(res, atname, atid, atmindex, atsymbol)
            id2atom[atid] = atom
            
            s = state[atmindex]
            s.t[1] = parse(T, line[31:38])
            s.t[2] = parse(T, line[39:46])
            s.t[3] = parse(T, line[47:54])
            atmindex += 1

        elseif startswith(line, "CONECT")
            idxs = map(s->parse(Int, s), [line[n:n+4] for n=7:5:length(line)])
            pivot = id2atom[idxs[1]]
            for i in idxs[2:end]
                bond(pivot, id2atom[i])
            end
        end
    end
    top.id = state.id = genid()
    
    # request conversion from cartesian to internal
    # request_c2i(state; all=true)
    
    Pose(top,state)
end


# Base.write(io::IO, top::Topology, state::State) = begin
Base.write(io::IO, top::AbstractContainer, state::State) = begin
    
    println(io, "MODEL")
    for atom in eachatom(top)
        sti = state[atom.index]
        # s = @sprintf("ATOM %6d %4s %-4sA %3d    %8.3f%8.3f%8.3f%24s",
        s = @sprintf("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%24s",
            atom.index, atom.name,
            atom.container.name, atom.container.container.name,
            atom.container.id,
            sti.t[1], sti.t[2], sti.t[3],
            atom.symbol)
        println(io, s)
    end

    for atom in eachatom(top)
       print(io, @sprintf("CONECT%5d", atom.index))
       foreach(n->print(io, @sprintf("%5d",n.index)), atom.children)
       println(io,"")
    end
    println(io, "ENDMDL")
end
