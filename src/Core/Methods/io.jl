using Printf: @sprintf
using YAML

const PDB = Val{1}
const YML = Val{2}

"""
    load([::Type{T}], filename::AbstractString) where {T <: AbstractFloat}

Load the given `filename` into a pose, parametrized by T. If this is not
provided, the default ProtoSyn.Units.defaultFloat is used. The file format is
infered from the extension (Supported: .pdb, .yml). If `bonds_by_distance` is
set to true (false, by default), the CONECT records will be complemented with
bonds infered by distance. The threshold distances for each pair of atoms is
defined in ProtoSyn.bond_lengths.
*NOTE:* This function does not infer any data of parenthood or ascedence. To
calculate that information, specific implementations of this function are
provided in other modules (such as `Peptides.load`).

# Examples
```julia-repl
julia> ProtoSyn.load("1ctf.pdb")
...
```
"""
function load(::Type{T}, filename::AbstractString; bonds_by_distance::Bool = false) where {T <: AbstractFloat}
    if endswith(filename, ".pdb") 
        return load(T, filename, PDB, bonds_by_distance = bonds_by_distance)
    elseif endswith(filename, ".yml")
        return load(T, filename, YML, bonds_by_distance = bonds_by_distance)
    else
        error("Unable to load '$filename': unsupported file type")
    end
end

load(filename::AbstractString; bonds_by_distance = false) = begin
    @info "Consider using Peptides.load when dealing with peptide chains."
    load(Float64, filename, bonds_by_distance = bonds_by_distance)
end

load(::Type{T}, filename::AbstractString, ::Type{K}; bonds_by_distance = false) where {T, K} = begin
    
    pose = load(T, open(filename), K)
    name, ext = splitext(basename(filename))
    pose.graph.name = name

    if bonds_by_distance
        dm        = ProtoSyn.Calculators.full_distance_matrix(pose)
        threshold = T(0.1)

        atoms   = collect(eachatom(pose.graph))
        for (i, atom_i) in enumerate(atoms)
            for (j, atom_j) in enumerate(atoms)
                i == j && continue
                atom_j = atoms[j]
                atom_j in atom_i.bonds && continue
                putative_bond = "$(atom_i.symbol)$(atom_j.symbol)"
                !(putative_bond in keys(bond_lengths)) && continue
                d = bond_lengths[putative_bond]
                d += d * threshold
                dm[i, j] < d && ProtoSyn.bond(atom_i, atom_j)
            end
        end
    end

    ProtoSyn.request_c2i(pose.state)
    # println()
    # sync!(pose)
    pose
end

load(filename::AbstractString, ::Type{K}; bonds_by_distance = false) where K = begin
    load(Float64, filename, K, bonds_by_distance = bonds_by_distance)
end

load(::Type{T}, io::IO, ::Type{YML}) where {T<:AbstractFloat} = begin
    
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

    root = origin(top)
    setparent!(
        res[graph["root"]],
        origin(top)
    )
    
    for atom in eachatom(top)
        atom.ascendents = ascendents(atom, 4)
    end

    request_i2c(state; all=true)
    top.id = state.id = genid()
    sync!(Pose(top, state))
end

load(::Type{T}, io::IO, ::Type{PDB}) where {T<:AbstractFloat} = begin
    
    top  = Topology("UNK", -1)
    seg  =  Segment("", -1)     # orphan segment
    res  =  Residue("", -1)     # orphan residue
    
    seekstart(io)
    natoms = mapreduce(l->startswith(l, "ATOM")|| startswith(l, "HETATM"), +, eachline(io); init=0)
    x = zeros(T, 3, natoms)
    
    id2atom = Dict{Int, Atom}()
    
    state = State{T}(natoms)
    
    segid = atmindex = 1
    
    seekstart(io)
    for line in eachline(io)
        
        if startswith(line, "TITLE")
            top.name = string(strip(line[11:end]))

        elseif startswith(line, "ATOM") || startswith(line, "HETATM")
            
            resname = string(strip(line[18:20]))
            segname = string(strip(string(line[22])))
            resid = parse(Int, line[23:26])

            if seg.name != segname
                seg = Segment!(top, segname, segid)
                seg.code = isempty(segname) ? '-' : segname[1]
                segid += 1
            end

            if res.id != resid || res.name != resname
                res = Residue!(seg, resname, resid)
                setparent!(res, ProtoSyn.origin(top).container)
            end

            atsymbol = length(line)>77 ? string(strip(line[77:78])) : "?"
            atname = string(strip(line[13:16]))
            atid = parse(Int, line[7:11])

            atom = Atom!(res, atname, atid, atmindex, atsymbol)
            id2atom[atid] = atom
            
            s = state[atmindex]
            xi = parse(T, line[31:38])
            s.t[1] = xi
            yi = parse(T, line[39:46])
            s.t[2] = yi
            zi = parse(T, line[47:54])
            s.t[3] = zi
            x[:, atmindex] = [xi, yi, zi]
            atmindex += 1

        elseif startswith(line, "CONECT")
            idxs = map(s -> parse(Int, s), split(line)[2:end])
            pivot = id2atom[idxs[1]]
            for i in idxs[2:end]
                other_atom = id2atom[i]
                bond(pivot, other_atom)
            end
        end
    end
    state.x = StateMatrix(state, x)
    top.id = state.id = genid()
    
    # request conversion from cartesian to internal
    Pose(top, state)
end

write(io::IO, top::AbstractContainer, state::State, ::Type{PDB}; model::Int = 1) = begin
    
    @printf(io, "MODEL %8d\n", model)
    for atom in eachatom(top)
        sti = state[atom.index]
        s = @sprintf("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%24s",
            atom.index, atom.name,
            atom.container.name, atom.container.container.code,
            atom.container.id,
            sti.t[1], sti.t[2], sti.t[3],
            atom.symbol)
        println(io, s)
    end

    for atom in eachatom(top)
       print(io, @sprintf("CONECT%5d", atom.index))
       foreach(n->print(io, @sprintf("%5d",n.index)), atom.bonds)
       println(io,"")
    end
    println(io, "ENDMDL")
end

write(io::IO, top::AbstractContainer, state::State, ::Type{YML}) = begin
    println(io, "name: ", top.name)
    println(io, "atoms:")
    
    byatom = eachatom(top)
    for at in byatom
        st = state[at]
        println(io,
            @sprintf("  - {name: %3s, id: %3d, symbol: %2s, b: %10.6f, theta: %10.6f, phi: %10.6f}",
            at.name, at.id, at.symbol, st.b, st.θ, st.ϕ)
        )
    end
    
    println(io, "bonds:")
    for at in byatom
        print(io, "  ", at.name, ": [")
        print(io, Base.join(map(a->a.name, at.bonds), ", "))
        println(io, "]")
    end

    println(io, "graph:")
    println(io, "  root: N")
    println(io, "  adjacency:")
    for at in byatom
        if haschildren(at)
            print(io, "    ", at.name, ": [")
            print(io, Base.join(map(a->a.name, at.children), ", "))
            println(io, "]")
        end
    end

end

"""
    ProtoSyn.write(pose::Pose, filename::String)

Write to file the given `pose`. The structure format is infered from the
`filename` extension (Supported: .pdb, .yml).

# Examples
```julia-repl
julia> ProtoSyn.write(pose, "1ctf_modified.pdb")
```
"""
function write(pose::Pose, filename::String)
    sync!(pose)
    io = open(filename, "w")
    if endswith(filename, ".pdb") 
        write(io, pose.graph, pose.state, PDB)
    elseif endswith(filename, ".yml")
        write(io, pose.graph, pose.state, YML)
    else
        error("Unable to write to '$filename': unsupported file type")
    end
    close(io)
end

"""
    ProtoSyn.append(pose::Pose, filename::String)

Append to file the given `pose` (as a new frame). The structure format is
infered from the `filename` extension (Supported: .pdb, .yml).

# Examples
```julia-repl
julia> ProtoSyn.append(pose, "1ctf_modified.pdb")
```
"""
function append(pose::Pose, filename::String; model::Int = 1)
    sync!(pose)
    io = open(filename, "a")
    if endswith(filename, ".pdb")
        write(io, pose.graph, pose.state, PDB, model = model)
    elseif endswith(filename, ".yml")
        write(io, pose.graph, pose.state, YML)
    else
        error("Unable to write to '$filename': unsupported file type")
    end
    close(io)
end