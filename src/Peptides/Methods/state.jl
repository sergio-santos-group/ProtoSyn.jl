using ProtoSyn.Units

"""
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
function setss!(state::State, ss::SecondaryStructureTemplate{T}, residues::Vector{Residue}) where {T <: AbstractFloat}
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
        ProtoSyn.setdihedral!(state, Dihedral.phi(r), ss.ϕ)
        # Last residues of chain might not have a psi angle.
        Dihedral.psi(r) !== nothing && ProtoSyn.setdihedral!(state, Dihedral.psi(r),  ss.ψ)
        ProtoSyn.setdihedral!(state, Dihedral.omega(r), ss.ω)
    end

    return state
end

setss!(pose::Pose, ss::SecondaryStructureTemplate{T}) where {T <: AbstractFloat} = begin
    residues::Vector{Residue} = collect(eachresidue(pose.graph))
    Peptides.setss!(pose.state, ss, residues)
end

setss!(pose::Pose, ss::SecondaryStructureTemplate{T}, sele::ProtoSyn.AbstractSelection) where {T <: AbstractFloat} = begin
    residues = ProtoSyn.PromoteSelection(sele, Residue, any)(pose.graph, gather = true)
    Peptides.setss!(pose.state, ss, residues)
end

setss!(pose::Pose, ss::SecondaryStructureTemplate{T}, residue::Residue) where {T <: AbstractFloat} = begin
    Peptides.setss!(pose.state, ss, [residue])
end