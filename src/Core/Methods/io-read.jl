
# --- YML ----------------------------------------------------------------------

load(::Type{T}, io::IO, ::Type{YML}; alternative_location::String = "A", kwargs...) where {T<:AbstractFloat} = begin
    
    yml = YAML.load(io)
    natoms = length(yml["atoms"])
    
    state = State{T}(natoms)
    top = Topology(yml["name"], 1)
    seg = Segment!(top, top.name, 1)
    res = Residue!(seg, top.name, 1)
    
    # add atoms
    for (index, pivot) in enumerate(yml["atoms"])
        atom = Atom!(res, pivot["name"], pivot["id"], index, pivot["symbol"])
        s = state[index]

        s.θ = ProtoSyn.Units.tonumber(T, pivot["theta"])
        s.ϕ = ProtoSyn.Units.tonumber(T, pivot["phi"])
        s.b = ProtoSyn.Units.tonumber(T, pivot["b"])

        # load charge
        if "c" in keys(pivot)
            s.δ = ProtoSyn.Units.tonumber(T, pivot["c"])
        end
    end

    # add bonds
    for (pivot, others) in yml["bonds"]
        atom = res[pivot]
        foreach(other -> bond(atom, res[other]), others)
    end

    # bond graph
    graph = yml["graph"]
    for (pivot, others) in graph["adjacency"]
        atom = res[pivot]
        foreach(other -> setparent!(res[other], atom), others)
    end

    setparent!(
        res[graph["root"]],
        ProtoSyn.root(top)
    )
    
    for atom in eachatom(top)
        atom.ascendents = ascendents(atom, 4)
    end

    request_i2c!(state; all=true)
    top.id = state.id = genid()
    sync!(Pose(top, state))
end

# --- PDB ----------------------------------------------------------------------

load(::Type{T}, io::IO, ::Type{PDB}; alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) where {T<:AbstractFloat} = begin
    
    top  = Topology("UNK", -1)
    seg  = Segment("", -1)     # orphan segment
    res  = Residue("", -1)     # orphan residue
    
    id2atom = Dict{Int, Atom}()
    
    state = State{T}() # empty state
    
    segid = atmindex = 1 # ! segment and atom index are overwritten by default 

    er = r"^\w+\s+(?<aid>\d+)\s+(?<an>\w{1,4}?)\s*(?<al>(?:\w{1}(?=\w{3}\s+\w{1}\s))|\s)(?<rn>\S{2,4})\s*(?<sn>[a-zA-Z]|[?]{0,1})\s*(?<rid>\d+)\s+(?<x>-*\d+\.\d+)\s+(?<y>-*\d+\.\d+)\s+(?<z>-*\d+\.\d+)\s+(?:\w)*\s+(?:\d+\.\d+)*\s+(?:\d+\.\d+)*\s*(?|(?:\w*\s+(\w)[0-9]*[+-]*\s*)|(?:(\w)[0-9]*[+-]*\s*)|(?<as>\w+[+-]*)(?=\s*$))$"

    aid = 0
    ignored_atoms     = Vector{Int}()
    uncommon_residues = Vector{Int}()
    conect_records_present = false
    for line in eachline(io)
        
        # if startswith(line, "TITLE")
        #     top.name = string(strip(line[11:end]))

        if startswith(line, "ATOM") || startswith(line, "HETATM")
            
            atom = match(er, line)

            if in(atom["rn"], ignore_residues) || in(atom["sn"], ignore_chains)
                push!(ignored_atoms, parse(Int, atom["aid"]))
                continue
            end

            # println(atom)

            # * Choose alternative locations
            al = string(atom["al"])
            if al !== " " && al !== alternative_location
                continue
            end

            segname   = atom["sn"] in [" ", ""] ? "?" : string(atom["sn"])             # * Default
            atom_symb = atom["as"] == ""  ? string(atom["an"][1]) : string(atom["as"]) # * Default

            if seg.name != segname
                seg = Segment!(top, segname, segid)
                seg.code = segname[1]
                segid += 1
            end

            resid = parse(Int, atom["rid"])
            resname = string(atom["rn"])
            if length(resname) > 3
                if !(resid in uncommon_residues)
                    @warn "Uncommon residue name found in residue $resid: $resname.\nProtoSyn.jl will use $(resname[1:3]) as the residue's name. Check if this is the desired behaviour."
                    push!(uncommon_residues, resid)
                end

                resname = resname[1:3]
            end

            # New residue
            if res.id != resid || res.name.content != resname
                res = Residue!(seg, resname, resid)
                setparent!(res, ProtoSyn.root(top).container)
                aid = 1
            end

            # Deal with repeated atom names
            an = string(atom["an"])
            while an in keys(res.itemsbyname)
                @warn "Atom named $an already found in residue $res. Adding atom identifier $aid."
                aid += 1
                an = an * string(aid)
            end

            new_atom = Atom!(res,
                an,
                parse(Int, atom["aid"]),
                atmindex,
                atom_symb)

            id2atom[parse(Int, atom["aid"])] = new_atom
            
            s = AtomState()
            s.t[1] = parse(T, atom["x"])
            s.t[2] = parse(T, atom["y"])
            s.t[3] = parse(T, atom["z"])
            push!(state, s)

        elseif startswith(line, "CONECT")
            conect_records_present = true
            idxs = map(s -> parse(Int, s), split(line)[2:end])

            # Consider ignored atoms in previous steps
            if idxs[1] in ignored_atoms
                continue
            end

            if !(idxs[1] in keys(id2atom))
                @warn "Found CONECT record $(idxs[1]) but not the corresponding atom."
                continue
            end
            pivot = id2atom[idxs[1]]
            for i in idxs[2:end]
                # Case this connect doesn't have a corresponding atom
                if !(i in keys(id2atom))
                    @warn "Found CONECT record $i but not the corresponding atom."
                    continue
                end
                other_atom = id2atom[i]
                bond(pivot, other_atom)
            end
        end
    end

    if !conect_records_present
        @warn "It seems the loaded pose does not have CONECT records present. Consider setting the `bonds_by_distance` flag to `true`."
    end

    top.id = state.id = genid()
   
    pose = Pose(top, state)
    reindex(pose.graph, set_ascendents = false)
    reindex(pose.state)
    return pose
end

# --- PQR ----------------------------------------------------------------------
load(::Type{T}, io::IO, ::Type{PQR}; alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) where {T<:AbstractFloat} = begin
    
    top  = Topology("UNK", -1)
    seg  = Segment("", -1)     # orphan segment
    res  = Residue("", -1)     # orphan residue
    
    id2atom = Dict{Int, Atom}()
    
    state = State{T}() # empty state
    
    segid = atmindex = 1 # ! segment and atom index are overwritten by default 

    er = r"\w+\s+(?<aid>\d+)\s+(?|((?:(?<an>\w{1,4})(?<al>\w))(?=\w{3}\s)(?<rn>\w{3}))|((\w+)\s+(\w?)(\w{3}))|((\w+)\s(\w*)\s(\w*)))\s+(?<sn>\D{1})\s*(?<rid>\d+)\s+(?<x>-*\d+\.\d+)\s+(?<y>-*\d+\.\d+)\s+(?<z>-*\d+\.\d+)\s+(?<c>-*\d+\.\d+)\s+(?:-*\d+\.\d+)"
    
    aid = 0
    ignored_atoms = Vector{Int}()
    for line in eachline(io)
        
        # if startswith(line, "TITLE")
        #     top.name = string(strip(line[11:end]))

        if startswith(line, "ATOM") || startswith(line, "HETATM")
            
            atom = match(er, line)

            if in(atom["rn"], ignore_residues) || in(atom["sn"], ignore_chains)
                push!(ignored_atoms, parse(Int, atom["aid"]))
                continue
            end

            # * Choose alternative locations
            al = string(atom["al"])
            if al !== "" && al !== alternative_location
                continue
            end

            segname = atom["sn"] == " " ? "?" : string(atom["sn"]) # * Default

            if seg.name != segname
                seg = Segment!(top, segname, segid)
                seg.code = isempty(segname) ? '-' : segname[1]
                segid += 1
            end

            resid = parse(Int, atom["rid"])
            resname = string(atom["rn"])

            # New residue
            if res.id != resid || res.name.content != resname
                res = Residue!(seg, resname, resid)
                setparent!(res, ProtoSyn.root(top).container)
                aid = 1
            end

            # Deal with repeated atom names
            an = string(atom["an"])
            while an in keys(res.itemsbyname)
                @warn "Atom named $an already found in residue $res. Adding atom identifier $aid."
                aid += 1
                an = an * string(aid)
            end

            new_atom = Atom!(res,
                an,
                parse(Int, atom["aid"]),
                atmindex,
                string(an[1]))

            id2atom[parse(Int, atom["aid"])] = new_atom
            
            s = AtomState()
            s.t[1] = parse(T, atom["x"])
            s.t[2] = parse(T, atom["y"])
            s.t[3] = parse(T, atom["z"])
            s.δ    = parse(T, atom["c"])
            push!(state, s)

        elseif startswith(line, "CONECT")
            idxs = map(s -> parse(Int, s), split(line)[2:end])

            # Consider ignored atoms in previous steps
            if idxs[1] in ignored_atoms
                continue
            end

            if !(idxs[1] in keys(id2atom))
                @warn "Found CONECT record $(idxs[1]) but not the corresponding atom."
                continue
            end
            pivot = id2atom[idxs[1]]
            for i in idxs[2:end]
                # Case this connect doesn't have a corresponding atom
                if !(i in keys(id2atom))
                    @warn "Found CONECT record $i but not the corresponding atom."
                    continue
                end
                other_atom = id2atom[i]
                bond(pivot, other_atom)
            end
        end
    end

    top.id = state.id = genid()
   
    pose = Pose(top, state)
    reindex(pose.graph, set_ascendents = false)
    reindex(pose.state)
    return pose
end


# --- XYZ ----------------------------------------------------------------------
load(::Type{T}, io::IO, ::Type{XYZ}; alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) where {T<:AbstractFloat} = begin
   
    top  = Topology("UNK", -1)
    seg  = Segment!(top, "UNK", -1)     # orphan segment
    res  = Residue!(seg, "UNK", -1)     # orphan residue
        
    state = State{T}() # empty state
        
    for (atmindex, line) in enumerate(eachline(io))
        elem = split(line)

        new_atom = Atom!(res,
            elem[1] * string(atmindex),
            atmindex,
            atmindex,
            string(elem[1]))
        
        s = AtomState()
        s.t[1] = parse(T, elem[2])
        s.t[2] = parse(T, elem[3])
        s.t[3] = parse(T, elem[4])
        push!(state, s)
    end

    top.id = state.id = genid()
   
    pose = Pose(top, state)
    reindex(pose.graph, set_ascendents = false)
    reindex(pose.state)
    return pose 
end

# --- SDF ----------------------------------------------------------------------
load(::Type{T}, io::IO, ::Type{SDF}; alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) where {T<:AbstractFloat} = begin
   
    id2atom = Dict{Int, Atom}()

    name = readline(io)
    top  = Topology(name, -1)
    seg  = Segment!(top, name, -1)     # orphan segment
    res  = Residue!(seg, name, -1)     # orphan residue
        
    state = State{T}() # empty state
        
    readline(io)
    readline(io)
    data = split(readline(io))
    N_atoms = parse(Int, data[1])
    N_bonds = parse(Int, data[2])

    for atmindex in 1:N_atoms
        line = readline(io)
        elem = split(line)

        new_atom = Atom!(res,
            elem[4] * string(atmindex),
            atmindex,
            atmindex,
            string(elem[4]))

        id2atom[atmindex] = new_atom
        
        s = AtomState()
        s.t[1] = parse(T, elem[1])
        s.t[2] = parse(T, elem[2])
        s.t[3] = parse(T, elem[3])
        push!(state, s)
    end

    for bondindex in 1:N_bonds
        line = readline(io)
        elem = split(line)

        pivot = id2atom[parse(Int, elem[1])]
        other = id2atom[parse(Int, elem[2])]
        bond(pivot, other)
    end

    top.id = state.id = genid()
   
    pose = Pose(top, state)
    reindex(pose.graph, set_ascendents = false)
    reindex(pose.state)
    return pose 
end