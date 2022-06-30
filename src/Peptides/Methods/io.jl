"""
    load([::Type{T}], filename::AbstractString; [bonds_by_distance::Bool = false], [alternative_location::String = "A"], [include_residues::Vector{String} = Vector{String}()], [ignore_residues::Vector{String} = Vector{String}()], [ignore_chains::Vector{String} = Vector{String}()], [sort_atoms_by_graph::Bool = false]) where {T <: AbstractFloat}

Load the given `filename` into a [`Pose`](@ref), parametrized by `T`. If this is
not provided, the default `ProtoSyn.Units.defaultFloat` is used instead. The
file format is infered from the extension (See `ProtoSyn.supported_formats` for
all supported formats). If `bonds_by_distance` is set to `true` (`false`, by
default), the CONECT records will be complemented with bonds infered by
distance. The threshold distances for each pair of atoms is defined in
`ProtoSyn.bond_lengths`. Infers parenthood and ascedence from bonds
(N-[`Residue`](@ref) instances have the connected C-[`Residue`](@ref) as child).
By default, and when available, ProtoSyn will use `alternative_location` `A`,
unless specified in the flag `alternative_location`. If the input file if of
type PDB and a trajectory, returns a vector of [`Pose`](@ref) instances instead.
If provided, ProtoSyn will ignore any [`Residue`](@ref) and [`Segment`](@ref)
instances in `ignore_residues` and `ignore_chains`, respectivelly (by name).
If a [`Residue`](@ref) instance is not found in
`ProtoSyn.Peptides.available_aminoacids` (for example, a ligand or a
non-canonical aminoacid), ProtoSyn will skip setting inter-residue parenthoods.
This can be overwritten by adding the [`Residue`](@ref) name in
`include_residues`, in which case ProtoSyn will attempt to find and established
all possible inter-residue parenthood relationships. In more complex cases this
should be manually checked. Despite the overall [`Pose`](@ref)'s [Graph](@ref graph-types),
[`Atom`](@ref) order in the encompassing `AbstractContainer` is, by default,
maintained the same as the input `filename`. This can se re-ordered to match the
infered parenthood by setting the `sort_atoms_by_graph` flag to `true` (`false`,
by default), in which case the [`ProtoSyn.sort_atoms_by_graph!`](@ref) method
will be used (employing the `ProtoSyn.Peptides.IUPAC` search algorithm).

!!! ukw "Note:"
    This function is an extension of [`ProtoSyn.load`](@ref).

# Examples
```
julia> ProtoSyn.Peptides.load("1ctf.pdb", bonds_by_distance = true)
Pose{Topology}(Topology{/2a3d:61708}, State{Float64}:
 Size: 1140
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function load(::Type{T}, filename::AbstractString; bonds_by_distance::Bool = false, alternative_location::String = "A", include_residues::Vector{String} = Vector{String}(), ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}(), sort_atoms_by_graph::Bool = false) where {T <: AbstractFloat}

    if !bonds_by_distance
        @info "Flag `bonds_by_distance` is set to False. Make sure the loaded $filename file has connect records."
    end

    poses = ProtoSyn.load(T, filename, bonds_by_distance = bonds_by_distance, alternative_location = alternative_location, ignore_chains = ignore_chains, ignore_residues = ignore_residues)

    available_aminoacids = [ProtoSyn.ProtoSyn.one_2_three[x] for x in keys(ProtoSyn.Peptides.available_aminoacids)]
    
    # Set parenthood of residues
    if isa(poses, Pose) 
        poses = Vector{Pose}([poses])
    end

    for pose in poses

        c = 0; for atom in eachatom(pose.graph)
            c += length(atom.bonds)
        end

        if c === 0
            @warn "Attempted to set-up parenthood in the Atom-level graph, but no bond records were found. Consider loading the pose with the `bonds_by_distance` flag set to `true`."
        else
            for segment in eachsegment(pose.graph)

                if string(segment.code) in ignore_chains
                    @warn "Skipping segment '$(segment.code)' (it's in 'ignore_chains' list)"
                    continue
                end

                n_residues = ProtoSyn.count_residues(segment)

                # Intra-residue parenthood (on residue 1)
                # Note: the pose already has an infered graph from ProtoSyn.load,
                # but the Peptides module makes sure that all intra-residue graphs
                # emerge from the N atom, matching the default Peptides grammar.
                residue = segment[1]

                residue.name.content in ignore_residues && continue

                starting_atom = ProtoSyn.identify_atom_by_bonding_pattern(residue, ["N", "C", "C", "O"])
                if isa(starting_atom, Vector{Atom})
                    if length(starting_atom) === 0
                        @warn "No starting atom found on residue $residue. Will use atom $(collect(eachatom(residue))[1]), check if this is the desired behaviour."
                        starting_atom = collect(eachatom(residue))[1]
                    else
                        @warn "Multiple starting atoms identified on residue $residue. Will use atom $(starting_atom[1]), check if this is the desired behaviour."
                        starting_atom = starting_atom[1]
                    end
                end
                ProtoSyn.infer_parenthood!(residue, overwrite = true, start = starting_atom)
                sort_atoms_by_graph && ProtoSyn.sort_atoms_by_graph!(pose.state, residue, start = starting_atom, search_algorithm = Peptides.IUPAC)

                for residue_index in 2:n_residues
                    residue = segment[residue_index]

                    residue.name.content in ignore_residues && continue

                    # Intra-residue parenthood
                    starting_atom = ProtoSyn.identify_atom_by_bonding_pattern(residue, ["N", "C", "C", "O"])
                    if isa(starting_atom, Vector{Atom})
                        if length(starting_atom) === 0
                            @warn "No starting atom found on residue $residue. Will use atom $(collect(eachatom(residue))[1]) to set intra-residue graph, check if this is the desired behaviour."
                            starting_atom = collect(eachatom(residue))[1]
                        else
                            @warn "Multiple starting atoms identified on residue $residue. Will use atom $(starting_atom[1]) to set intra-residue graph, check if this is the desired behaviour."
                            starting_atom = starting_atom[1]
                        end
                    end
                    ProtoSyn.infer_parenthood!(residue, overwrite = true, start = starting_atom)
                    sort_atoms_by_graph && ProtoSyn.sort_atoms_by_graph!(pose.state, residue, start = starting_atom, search_algorithm = Peptides.IUPAC)

                    # Inter-residue parenthood
                    if !(residue.name in available_aminoacids) && !(residue.name in include_residues)
                        @warn "Found a possible ligand or unknown NCAA at residue $residue. ProtoSyn will skip it when setting up inter residue parenthoods. This behaviour can be changed by adding this residue name to the 'include_residues' vector."
                        continue
                    end
                    
                    # Re-set the inter-residue parenthood to match the default
                    # Peptides grammar
                    popparent!(residue) # Was root before.
                    popparent!(starting_atom) # Was root before.
                    # println("Bonds of $(starting_atom): $(starting_atom.bonds)")
                    C_bond_index = findfirst(x -> x.symbol == "C", starting_atom.bonds)
                    # ! Note: It may happen that there is only the CA connected to
                    # ! an N. In such case, the residue parent will be itself. These
                    # ! cases are adressed further down.
                    @assert C_bond_index !== nothing "No atom \"C\" found in atom $(starting_atom) bonds list ($(starting_atom.bonds))."
                    parent_residue = starting_atom.bonds[C_bond_index].container
                    setparent!(residue, parent_residue)
                    setparent!(starting_atom, starting_atom.bonds[C_bond_index])
                    # println("$(starting_atom.bonds[C_bond_index])-$(starting_atom)")
                end
            end
        end
        
        # Adress cases where the parent's residue is itself (happens when
        # there's missing residues in the chain)
        for residue in eachresidue(pose.graph)
            if residue.parent === residue
                popparent!(residue)
                setparent!(residue, ProtoSyn.root(pose.graph).container)
            end
        end

        reindex(pose.graph)
        ProtoSyn.request_c2i!(pose.state; all = true)
        sync!(pose)
    end

    if length(poses) == 1
        poses = poses[1]
    end

    return poses
end

load(filename::AbstractString; bonds_by_distance::Bool = false, alternative_location::String = "A", include_residues::Vector{String} = Vector{String}(), ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}(), sort_atoms_by_graph::Bool = false) = begin
    Peptides.load(ProtoSyn.Units.defaultFloat, filename, bonds_by_distance = bonds_by_distance, alternative_location = alternative_location, include_residues = include_residues, ignore_residues = ignore_residues, ignore_chains = ignore_chains, sort_atoms_by_graph = sort_atoms_by_graph)
end


"""
    ProtoSyn.Peptides.download([::Type{T}], pdb_code::String; [bonds_by_distance::Bool = false], [include_residues::Vector{String} = Vector{String}()], [ignore_residues::Vector{String} = Vector{String}()], [ignore_chains::Vector{String} = Vector{String}()], [sort_atoms_by_graph::Bool = true]) where {T <: AbstractFloat}

Download the PDB file (for the given PDB code) from the RCSB
Protein Data Bank into a [`Pose`](@ref). The downloaded file can be found in the
current working directory. If `T` is specified, the downloaded file will be
loaded into a [`Pose`](@ref) parametrized by `T`, otherwise uses the default
`ProtoSyn.Units.defaultFloat`. Uses the specific `Peptides.load` method. For
more information on the input arguments, check the
[`load`](@ref ProtoSyn.Peptides.load) method documentation.

# Examples
```jldoctest
julia> ProtoSyn.download("2A3D")
```
"""
function ProtoSyn.Peptides.download(::Type{T}, pdb_code::String; bonds_by_distance::Bool = false, include_residues::Vector{String} = Vector{String}(), ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}(), sort_atoms_by_graph::Bool = false) where {T <: AbstractFloat}
    if endswith(pdb_code, ".pdb"); pdb_code = pdb_code[1:(end - 4)]; end
    filename = pdb_code * ".pdb"
    url = "https://files.rcsb.org/download/" * filename
    Base.download(url, filename)
    return ProtoSyn.Peptides.load(T, filename, bonds_by_distance = bonds_by_distance, include_residues = include_residues, ignore_residues = ignore_residues, ignore_chains = ignore_chains, sort_atoms_by_graph = sort_atoms_by_graph)
end

ProtoSyn.Peptides.download(pdb_code::String; bonds_by_distance::Bool = false, include_residues::Vector{String} = Vector{String}(), ignore_residues::Vector{String} = Vector{String}(), ignore_chains::Vector{String} = Vector{String}(), sort_atoms_by_graph::Bool = true) = begin
    ProtoSyn.Peptides.download(ProtoSyn.Units.defaultFloat, pdb_code, bonds_by_distance = bonds_by_distance, include_residues = include_residues, ignore_residues = ignore_residues, ignore_chains = ignore_chains, sort_atoms_by_graph = sort_atoms_by_graph)
end