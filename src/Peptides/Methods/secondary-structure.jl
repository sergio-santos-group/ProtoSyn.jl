using ProtoSyn.Units

"""
# TODO
# * Add include variation explanation

    setss!(state::State, ss::SecondaryStructureTemplate, residues::Vector{Residue})

Set the phi `ϕ`, psi `ψ` and omega `ω` backbone angles of all [`Residue`](@ref)
instances in the given `residues` vector to match the provided
[`SecondaryStructureTemplate`](@ref). This function acts on the internal
coordinates and does not update cartesian coordinates (using the
[`sync!`](@ref ProtoSyn.sync!) method), although a request for conversion is
made (by calling the [`request_i2c!`](@ref ProtoSyn.request_i2c!)).

    setss!(pose::Pose, ss::SecondaryStructureTemplate)
    setss!(pose::Pose, ss::SecondaryStructureTemplate, sele::ProtoSyn.AbstractSelection)
    setss!(pose::Pose, ss::SecondaryStructureTemplate, residue::Residue)

Set the phi `ϕ`, psi `ψ` and omega `ω` backbone angles of all [`Residue`](@ref)
instances in the given [`Pose`](@ref) `pose` to match the provided
[`SecondaryStructureTemplate`](@ref). If an optional `AbstractSelection` `sele`
is provided, apply the transformation only to the selected [`Residue`](@ref)
instances.

# See also
[`SecondaryStructure`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.setss!(pose, ProtoSyn.Peptides.SecondaryStructure[:helix])
State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
```
"""
function setss!(state::State, ss::SecondaryStructureTemplate{T}, residues::Vector{Residue}; include_variation::Bool = false, min_prob::T = 0.0) where {T <: AbstractFloat}
    for r in residues
        if r.name == "PRO"
            # Proline is restricted to TRANS conformation
            # This conformation is most abundant (~95%) in globular proteins.
            # https://pubs.acs.org/doi/10.1021/jacs.0c02263
            # PHI: -75° | PSI: 145° | OMEGA: 180°
            ProtoSyn.setdihedral!(state, Dihedral.phi(r), T(-75°))
            # Last residues of chain might not have a psi angle.
            Dihedral.psi(r) !== nothing && ProtoSyn.setdihedral!(state, Dihedral.psi(r),  T(145°))
            ProtoSyn.setdihedral!(state, Dihedral.omega(r), T(180°))
            continue
        end
        if !include_variation
            ProtoSyn.setdihedral!(state, Dihedral.phi(r), ss.ϕ.angle)
            # Last residues of chain might not have a psi angle.
            Dihedral.psi(r) !== nothing && ProtoSyn.setdihedral!(state, Dihedral.psi(r),  ss.ψ.angle)
            ProtoSyn.setdihedral!(state, Dihedral.omega(r), ss.ω.angle)
        else
            _ϕ = ss.ϕ.sampler === nothing ? ss.ϕ.angle : Peptides.sample(ss.ϕ.sampler, min_prob = min_prob)
            ProtoSyn.setdihedral!(state, Dihedral.phi(r), _ϕ)
            # Last residues of chain might not have a psi angle.
            _ψ = ss.ψ.sampler === nothing ? ss.ψ.angle : Peptides.sample(ss.ψ.sampler, min_prob = min_prob)
            Dihedral.psi(r) !== nothing && ProtoSyn.setdihedral!(state, Dihedral.psi(r), _ψ)
            _ω = ss.ω.sampler === nothing ? ss.ω.angle : Peptides.sample(ss.ω.sampler, min_prob = min_prob)
            ProtoSyn.setdihedral!(state, Dihedral.omega(r), _ω)
        end
    end

    return state
end

setss!(pose::Pose, ss::SecondaryStructureTemplate{T}; include_variation::Bool = false, min_prob::T = 0.0) where {T <: AbstractFloat} = begin
    residues::Vector{Residue} = collect(eachresidue(pose.graph))
    Peptides.setss!(pose.state, ss, residues, include_variation = include_variation, min_prob = min_prob)
end

setss!(pose::Pose, ss::SecondaryStructureTemplate{T}, sele::ProtoSyn.AbstractSelection; include_variation::Bool = false, min_prob::T = 0.0) where {T <: AbstractFloat} = begin
    residues = ProtoSyn.PromoteSelection(sele, Residue, any)(pose.graph, gather = true)
    Peptides.setss!(pose.state, ss, residues, include_variation = include_variation, min_prob = min_prob)
end

setss!(pose::Pose, ss::SecondaryStructureTemplate{T}, residue::Residue; include_variation::Bool = false, min_prob::T = 0.0) where {T <: AbstractFloat} = begin
    Peptides.setss!(pose.state, ss, [residue], include_variation = include_variation, min_prob = min_prob)
end


# ---

"""
# TODO
"""
function read_ss_map(pose::Pose, filename::String, ss_type::String)

    sequence = ProtoSyn.Peptides.sequence(pose)
    N = length(sequence)
    i = 0
    selection = !TrueSelection{Residue}()
    
    ss_map = open(filename, "r")

    for line in eachline(ss_map)
        startswith(line, "#") && continue
        elems = split(line)

        
        if (i+1) > N
            @info "Secondary structure map continues after last aminoacid on sequence (> $N): skipping this aminoacid on the secondary structure map."
            continue
        elseif string(elems[2]) !== string(sequence[i + 1])
            @info "Aminoacid type on secondary structure map ($(elems[2])) does not match the given sequence on the same position ($(sequence[i + 1])): skiping this aminoacid on the secondary structure map."
            continue
        end
        
        i += 1

        if string(elems[3]) === ss_type
            id = parse(Int, elems[1])
            this_aminoacid = ProtoSyn.SerialSelection{Residue}(id, :id)
            selection = selection | this_aminoacid
        end
    end
    
    close(ss_map)
    
    return selection
end


function read_ss_map(filename::String, ss_type::String)
    
    selection = !TrueSelection{Residue}()

    open(filename, "r") do ss_map
        for line in eachline(ss_map)
            startswith(line, "#") && continue
            elems = split(line)

            if string(elems[3]) === ss_type
                id = parse(Int, elems[1])
                this_aminoacid = ProtoSyn.SerialSelection{Residue}(id, :id)
                selection = selection | this_aminoacid
            end
        end
    end

    return selection
end