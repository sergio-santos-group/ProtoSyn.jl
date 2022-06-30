using Printf: @sprintf
using YAML
using Downloads: download

const PDB = Val{1}
const YML = Val{2}
const PQR = Val{3}
const XYZ = Val{4}

const supported_formats = [
    ".PDB => Read: ✓ | Write: ✓",
    ".YML => Read: ✓ | Write: ✓",
    ".PQR => Read: ✗ | Write: ✓",
    ".XYZ => Read: ✗ | Write: ✓"]

"""
    load([::Type{T}], filename::AbstractString, [bonds_by_distance::Bool = false], [alternative_location::String = "A"], [ignore_residues::Vector{String} = Vector{String}()], [ignore_chains::Vector{String} = Vector{String}()]) where {T <: AbstractFloat}

Load the given `filename` into a [`Pose`](@ref), parametrized by `T`. If this is
not provided, the default `ProtoSyn.Units.defaultFloat` is used. The file format
is infered from the extension (See `ProtoSyn.supported_formats` for all
supported formats). If `bonds_by_distance` is set to `true` (`false`, by
default), the CONECT records will be complemented with bonds infered by
distance. The distances for each pair of atoms is defined in
`ProtoSyn.Units.bond_lengths` (in Angstrom Å, with a standard deviation
threshold of 0.1 Å). Return the resulting [`Pose`](@ref) instance. By default,
and when available, ProtoSyn will use `alternative_location` `A`, unless
specified in the flag `alternative_location`. If the input file if of type PDB
and a trajectory, returns a vector of [`Pose`](@ref) instances instead.
Optionally, by setting `ignore_residues` and `ignore_chains`, ProtoSyn will skip
the load of any atom belonging to either the given residues or chains (by name).

# See also
[`distance`](@ref)

!!! ukw "Note:"
    This function tries to infer information of the parenthood and ascendents of each atom, using the CONECT records or infered `bonds_by_distance`. The parents are arbitrarily defined as the first bond found, by order, and any atom without parent is connected to the [`root`](@ref ProtoSyn.root). All [`Residue`](@ref) instances have the [`root`](@ref ProtoSyn.root)`.container` as parent. Note that this infered information may need to be manually corrected.

# Examples
```
julia> ProtoSyn.load("2a3d.pdb")
Pose{Topology}(Topology{/2a3d:6263}, State{Float64}:
 Size: 1140
 i2c: false | c2i: true
 Energy: Dict(:Total => Inf)
)
```
"""
function load(::Type{T}, filename::AbstractString; bonds_by_distance::Bool = false, alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) where {T <: AbstractFloat}
    if endswith(filename, ".pdb")
        # Check if this is a trajectory
        if !is_trajectory(filename)
            return load(T, filename, PDB, bonds_by_distance = bonds_by_distance, alternative_location = alternative_location, ignore_chains = ignore_chains, ignore_residues = ignore_residues)
        else
            return load_trajectory(T, filename, PDB, bonds_by_distance = bonds_by_distance, alternative_location = alternative_location, ignore_chains = ignore_chains, ignore_residues = ignore_residues)
        end
    elseif endswith(filename, ".yml")
        return load(T, filename, YML, bonds_by_distance = bonds_by_distance, infer_parenthood = false)
    elseif endswith(filename, ".xyz")
        return load(T, filename, XYZ, bonds_by_distance = bonds_by_distance, infer_parenthood = false)
    elseif endswith(filename, ".pqr")
        return load(T, filename, PQR, bonds_by_distance = bonds_by_distance, infer_parenthood = true)
    else
        error("Unable to load '$filename': unsupported file type")
    end
end

load(filename::AbstractString; bonds_by_distance = false, alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) = begin
    @info "Consider using Peptides.load when dealing with peptide chains."
    load(ProtoSyn.Units.defaultFloat, filename, bonds_by_distance = bonds_by_distance; alternative_location = alternative_location, ignore_chains = ignore_chains, ignore_residues = ignore_residues)
end


"""
    is_trajectory(filename::String)

Read the given `filename` and check if multiple "MODEL" entries are found. File
must be in PDB format.

# Examples
```jldoctest
julia> ProtoSyn.is_trajectory("teste.pdb")
true
```
"""
function is_trajectory(filename::String)
    
    @assert filename[(end-3):end] == ".pdb" "File $filename must be in PDB format."
    @assert isfile(filename) "File \"$filename\" not found."

    io = open(filename, "r")
    models = 0
    for line in eachline(io)
        if startswith(line, "MODEL")
            models += 1
            models > 1 && return true
        end
    end

    return false
end


"""
    splice_trajectory(filename::String)

Create a new temporary folder holding all the "MODEL" entries in a given
`filename` separated (one per file, input should be in PDB format).

# Examples
```
julia> ProtoSyn.splice_trajectory("teste.pdb")
"teste_spliced"
```
"""
function splice_trajectory(filename::String)

    @assert filename[(end-3):end] == ".pdb" "File $filename must be in PDB format."

    # Create temporary folder
    dirname::String = "$(filename[1:(end-4)])_spliced"
    isdir(dirname) && begin
        @info @info "Found pre-existent $dirname folder. Overwritting."
        rm(dirname, recursive = true)
    end
    mkdir(dirname)

    model::Int = 0
    io_out = nothing
    io = open(filename, "r")

    for line in eachline(io)
        if startswith(line, "MODEL")
            model += 1
            model_name = joinpath(dirname, "$model.pdb")
            io_out !== nothing && close(io_out)
            io_out = open(model_name, "w")
        end

        if io_out !== nothing
            Base.write(io_out, line * "\n")
        end
    end
    close(io_out)

    return dirname
end


"""
    load_trajectory([::Type{T}], filename::AbstractString, ::Type{K}; [bonds_by_distance = false], [alternative_location::String = "A"], [ignore_residues::Vector{String} = Vector{String}()], [ignore_chains::Vector{String} = Vector{String}()]) where {T <: AbstractFloat, K}

Load the given `filename` into a vector of [`Pose`](@ref) instances,
parametrized by `T`, separated by new "MODEL" entries. If `T` is not provided,
the default `ProtoSyn.Units.defaultFloat` is used. The file format is infered
from the extension (Supported: .pdb only). If `bonds_by_distance` is set to
`true` (`false`, by default), the CONECT records will be complemented with bonds
infered by distance. The distances for each pair of atoms is defined in
`ProtoSyn.Units.bond_lengths` (in Angstrom Å, with a standard deviation
threshold of 0.1 Å). Return the resulting vector of [`Pose`](@ref) instances. By
default, and when available, ProtoSyn will use `alternative_location` `A`,
unless specified in the flag `alternative_location`. Optionally, by setting
`ignore_residues` and `ignore_chains`, ProtoSyn will skip the load of any atom
belonging to either the given residues or chains (by name).

# See also
[`is_trajectory`](@ref) [`splice_trajectory`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.load_trajectory("teste.pdb")
2-element Vector{Pose}:
 Pose{Topology}(Topology{/1:5584}, State{Float64}:
 Size: 39
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
 Pose{Topology}(Topology{/2:48484}, State{Float64}:
 Size: 39
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function load_trajectory(::Type{T}, filename::AbstractString, ::Type{K}; bonds_by_distance = false, alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) where {T <: AbstractFloat, K}
    models = splice_trajectory(filename)
    a = readdir(models)
    files = sort([parse(Int, split(x, ".")[1]) for x in a if endswith(x, ".pdb")])
    files = map(x -> joinpath(models, string(x) * ".pdb"), files)

    poses = Vector{Pose}()
    N = length(files)
    for (i, file) in enumerate(files)
        @info "Loading pose $i out of $N"
        pose = load(T, file, K; bonds_by_distance = bonds_by_distance, alternative_location = alternative_location, ignore_chains = ignore_chains, ignore_residues = ignore_residues)
        push!(poses, pose)
    end

    rm(models, recursive = true)

    return poses
end

load_trajectory(filename::AbstractString, ::Type{K}; bonds_by_distance = false, alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) where {T <: AbstractFloat, K} = begin
    load_trajectory(ProtoSyn.Units.defaultFloat, filename, K, bonds_by_distance = bonds_by_distance, alternative_location = alternative_location, ignore_chains = ignore_chains, ignore_residues = ignore_residues)
end

load_trajectory(filename; bonds_by_distance::Bool = false, alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) = begin
    load_trajectory(filename, PDB, bonds_by_distance = bonds_by_distance, alternative_location = alternative_location, ignore_chains = ignore_chains, ignore_residues = ignore_residues)
end

# The function bellow is similar to the first `load`, excepts it already
# receives the ::Type{K}. The only purpose is to perform `bonds by distance` and
# `infer_parenthood!` operations. In future version of ProtoSyn, the two methods
# might be combined in 1.
load(::Type{T}, filename::AbstractString, ::Type{K}; bonds_by_distance = false, alternative_location::String = "A", infer_parenthood::Bool = true, ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) where {T <: AbstractFloat, K} = begin
    
    pose = load(T, open(filename), K; alternative_location = alternative_location, ignore_chains = ignore_chains, ignore_residues = ignore_residues)
    name, _ = splitext(basename(filename))
    pose.graph.name = name

    bonds_by_distance && infer_bonds!(pose, threshold = 0.1)

    # Set parenthood of atoms (infered)
    infer_parenthood && infer_parenthood!(pose.graph)

    ProtoSyn.request_c2i!(pose.state)
    sync!(pose)
    pose
end

load(filename::AbstractString, ::Type{K}; bonds_by_distance = false, alternative_location::String = "A", ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}()) where K = begin
    load(Float64, filename, K, bonds_by_distance = bonds_by_distance, alternative_location = alternative_location, ignore_chains = ignore_chains, ignore_residues = ignore_residues)
end

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

    er = r"\w+\s+(?<aid>\d+)\s+(?|((?:(?<an>\w{1,4})(?<al>\w))(?=\w{3}\s)(?<rn>\w{3}))|((\w+)\s+(\w?)(\w{3}))|((\w+)\s(\w*)\s(\w*)))\s+(?<sn>\D{1})\s*(?<rid>\d+)\s+(?<x>-*\d+\.\d+)\s+(?<y>-*\d+\.\d+)\s+(?<z>-*\d+\.\d+)\s+(?:\w)*\s+(?:\d+\.\d+)*\s+(?:\d+\.\d+)*\s+(?|(?<as>\w+[+-]*)(?=\s*$)|(?:\w*\s+(\w)[0-9]*[+-]*\s*)|(\s*))"

    aid = 0
    ignored_atoms = Vector{Int}()
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
                string(atom["as"]))

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

# --- WRITE --------------------------------------------------------------------
# --- PDB ----------------------------------------------------------------------

write(io::IO, top::AbstractContainer, state::State, ::Type{PDB}; model::Int = 1, B_factors::Vector{Float64} = Float64[]) = begin
    
    B_factors_exist = length(B_factors) > 0
    @printf(io, "MODEL %8d\n", model)
    for (index, segment) in enumerate(eachsegment(top))
        index > 1 && println(io, "TER")
        for (i, atom) in enumerate(eachatom(segment))
            sti = state[atom.index]

            B_factor = B_factors_exist ? B_factors[i] : 0.0

            s = @sprintf("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%12.2f%12s",
                atom.index, atom.name,
                atom.container.name, atom.container.container.code,
                atom.container.id,
                sti.t[1], sti.t[2], sti.t[3], B_factor,
                atom.symbol)
            println(io, s)
        end
    end

    for atom in eachatom(top)
       print(io, @sprintf("CONECT%5d", atom.index))
       foreach(n->print(io, @sprintf("%5d", n.index)), atom.bonds)
       println(io,"")
    end
    println(io, "ENDMDL")
end

# --- PQR ----------------------------------------------------------------------

write(io::IO, top::AbstractContainer, state::State, ::Type{PQR}; model::Int = 1) = begin
    
    @printf(io, "MODEL %8d\n", model)
    for (index, segment) in enumerate(eachsegment(top))
        index > 1 && println(io, "TER")
        for atom in eachatom(segment)
            sti = state[atom.index]

            s = @sprintf("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%8.3f%8.3f",
                atom.index, atom.name,
                atom.container.name, atom.container.container.code,
                atom.container.id,
                sti.t[1], sti.t[2], sti.t[3], sti.δ, 0.0)
            println(io, s)
        end
    end

    for atom in eachatom(top)
       print(io, @sprintf("CONECT%5d", atom.index))
       foreach(n->print(io, @sprintf("%5d", n.index)), atom.bonds)
       println(io,"")
    end
    println(io, "ENDMDL")
end

# --- YML ----------------------------------------------------------------------

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
    println(io, "  root: $(ProtoSyn.root(top).children[1].name)")
    println(io, "  adjacency:")
    for at in byatom
        if !isempty(at.children)
            print(io, "    ", at.name, ": [")
            print(io, Base.join(map(a->a.name, at.children), ", "))
            println(io, "]")
        end
    end

end

# --- XYZ ----------------------------------------------------------------------

write(io::IO, top::AbstractContainer, state::State, ::Type{XYZ}) = begin
    
    for atom in eachatom(top)
        sti = state[atom.index]
        s = @sprintf("%3s %8.3f %8.3f %8.3f", atom.symbol, sti.t[1], sti.t[2],
            sti.t[3])
        println(io, s)
    end
end

# The following two function allow ProtoSyn to write a
# Vector{T <: AbstractFloat} in a semi-XYZ format (atoms are set to be 'X').
# This allows for a quick visualization is common molecular viz tools.
write(io::IO, v::Vector{Vector{T}}) where {T <: AbstractFloat} = begin
    
    for atom in v
        s = @sprintf("%3s %8.3f %8.3f %8.3f", "X", atom[1], atom[2], atom[3])
        println(io, s)
    end
end

write(filename::String, v::Vector{Vector{T}}) where {T <: AbstractFloat} = begin
    if split(filename, ".")[2] != "xyz"
        @error "ProtoSyn doesn't allow Vector{Vector{T}} print to formats other than .xyz!"
    end
    open(filename, "w") do io
        write(io, v)
    end
end


"""
    ProtoSyn.write(pose::Pose, filename::String)

Write to file the given [`Pose`](@ref) `pose`. The file format is infered from
the `filename` extension (See `ProtoSyn.supported_formats` for all supported
formats). The [`Pose`](@ref) `pose` structure is automatically synched (using
the [`sync!`](@ref) method) when writting to file, as only the cartesian
coordinates are used.

# See also
[`append`](@ref)

# Examples
```
julia> ProtoSyn.write(pose, "new_file.pdb")
```
"""
function write(pose::Pose{Topology}, filename::String)
    sync!(pose)
    io = open(filename, "w")
    if endswith(filename, ".pdb")
        write(io, pose.graph, pose.state, PDB)
    elseif endswith(filename, ".yml")
        write(io, pose.graph, pose.state, YML)
    elseif endswith(filename, ".pqr")
        write(io, pose.graph, pose.state, PQR)
    elseif endswith(filename, ".xyz")
        write(io, pose.graph, pose.state, XYZ)
    else
        error("Unable to write to '$filename': unsupported file type")
    end
    close(io)
end

function write(pose::Pose{Segment}, filename::String)
    error("MethodError: no method ProtoSyn.write is available for Fragment instances. Consider generating a Pose from the fragment, using `Pose(frag)`.")
end


"""
    ProtoSyn.append(pose::Pose, filename::String, [model::Int = 1])

Append to file the given [`Pose`](@ref) `pose` (as a new frame, identified by
the model number `model`: default is 1). The file format is infered from the
`filename` extension (See `ProtoSyn.supported_formats` for all supported
formats). The [`Pose`](@ref) `pose` structure is automatically synched (using
the [`sync!`](@ref) method) when writting to file, as only the cartesian
coordinates are used.

# See also
[`write`](@ref)

# Examples
```
julia> ProtoSyn.append(pose, "new_file.pdb")
```
"""
function append(pose::Pose, filename::String; model::Int = 1)
    sync!(pose)
    io = open(filename, "a")
    if endswith(filename, ".pdb")
        write(io, pose.graph, pose.state, PDB, model = model)
    elseif endswith(filename, ".yml")
        write(io, pose.graph, pose.state, YML)
    elseif endswith(filename, ".pqr")
        write(io, pose.graph, pose.state, PQR)
    elseif endswith(filename, ".xyz")
        write(io, pose.graph, pose.state, XYZ)
    else
        error("Unable to write to '$filename': unsupported file type")
    end
    close(io)
end


"""
    ProtoSyn.download([::T], pdb_code::String) where {T <: AbstractFloat}

Download the PDB file (for the given PDB code) from the RCSB
Protein Data Bank into a [`Pose`](@ref). The downloaded file can be found in the
current working directory. If `T` is specified, the downloaded file will be
loaded into a [`Pose`](@ref) parametrized by `T`, otherwise uses the default
`ProtoSyn.Units.defaultFloat`.

# See also
[`load`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.download("2A3D")
```
"""
function ProtoSyn.download(::Type{T}, pdb_code::String; bonds_by_distance::Bool = false) where {T <: AbstractFloat}
    if endswith(pdb_code, ".pdb"); pdb_code = pdb_code[1:(end - 4)]; end
    filename = pdb_code * ".pdb"
    url = "https://files.rcsb.org/download/" * filename
    download(url, filename)
    return load(T, filename, bonds_by_distance = bonds_by_distance)
end

ProtoSyn.download(pdb_code::String; bonds_by_distance::Bool = false) = begin
    ProtoSyn.download(ProtoSyn.Units.defaultFloat, pdb_code; bonds_by_distance = bonds_by_distance)
end