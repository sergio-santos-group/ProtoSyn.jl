using ProtoSyn
using ProtoSyn.Units
using CurveFit: Polynomial

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


"""
# TODO: Documentation
# Compare with DSSP
"""
function catergorize_ss_from_dihedral_angles(pose::Pose, selection::Opt{ProtoSyn.AbstractSelection})
    T = eltype(pose.state)

    if selection === nothing
        sele = TrueSelection{Residue}()
    else
        sele = promote(selection, Residue)
    end

    residues = sele(pose, gather = true)

    N = length(residues)
    points = zeros((3, N))

    ss_assignment = Dict(1 => "H", 2 => "E", 3 => "C")

    # 1. Get geometrical criteria ----------------------------------------------

    phi_potentials = Dict{String, Polynomial}([
        "1L" => ProtoSyn.Peptides.phi_α_L_potential,
        "1R" => ProtoSyn.Peptides.phi_α_R_potential,
        "2" => ProtoSyn.Peptides.phi_β_potential,
        "3" => ProtoSyn.Peptides.phi_coil_potential])

    phi_values = Dict{String, T}([
        "1L" => 0.0, # Helix
        "1R" => 0.0, # Helix
        "2" => 0.0,  # Beta-sheet
        "3" => 0.0]) # Coil

    psi_potentials = Dict{String, Polynomial}([
        "1L" => ProtoSyn.Peptides.psi_α_L_potential,
        "1R" => ProtoSyn.Peptides.psi_α_R_potential,
        "2" => ProtoSyn.Peptides.psi_β_potential,
        "3" => ProtoSyn.Peptides.psi_coil_potential])

    psi_values = Dict{String, T}([
        "1L" => 0.0, # Helix
        "1R" => 0.0, # Helix
        "2" => 0.0,  # Beta-sheet
        "3" => 0.0]) # Coil

    for (residue_index, residue) in enumerate(residues)
        phi_atom = ProtoSyn.Peptides.Dihedral.phi(residue)
        phi_atom === nothing && continue
        phi = ProtoSyn.getdihedral(pose.state, phi_atom)
        
        for (ss, potential) in phi_potentials
            value = potential(phi)
            phi_values[ss] = value
        end
        min_phi_values = argmin(collect(values(phi_values)))
        phi_ss = parse(Int, collect(keys(phi_values))[min_phi_values][1])
        
        psi_atom = ProtoSyn.Peptides.Dihedral.psi(residue)
        psi_atom === nothing && continue
        psi = ProtoSyn.getdihedral(pose.state, psi_atom)
        for (ss, potential) in psi_potentials
            value = potential(psi)
            psi_values[ss] = value
        end
        min_psi_values = argmin(collect(values(psi_values)))
        psi_ss = parse(Int, collect(keys(psi_values))[min_psi_values][1])

        if phi_ss === psi_ss
            points[phi_ss, residue_index] += 2.0
            points[psi_ss, residue_index] += 2.0
        else
            points[phi_ss, residue_index] += 1.0
            points[psi_ss, residue_index] += 1.0
            points[3, residue_index] += 0.5
        end
    end

    # 2. Get hydrogen bonding pattern criteria ---------------------------------

    hbn = ProtoSyn.Calculators.HydrogenBonds.generate_hydrogen_bond_network(
        pose, sele & !SidechainSelection())
    _, _, hb_list = ProtoSyn.Calculators.HydrogenBonds.calc_hydrogen_bond_network(
        pose, nothing, false; hydrogen_bond_network = hbn)

    involved_in_hb = falses(N)
    for (acceptor, donor) in hb_list
        di = donor.container.index
        ai = acceptor.container.index
        involved_in_hb[di] = true
        involved_in_hb[ai] = true
        if (di === ai + 4) | (di === ai + 3)
            points[1, di] += 4.0
            points[1, ai] += 4.0
        elseif di === ai + 1
            points[3, di] += 4.0
            points[3, ai] += 4.0
        else
            points[2, di] += 2.0
            points[2, ai] += 2.0
            points[3, di] += 1.0
            points[3, ai] += 1.0
        end
    end

    for (residue_index, _involved_in_hb) in enumerate(involved_in_hb)
        if !_involved_in_hb
            points[1, residue_index] -= 2.0
            points[2, residue_index] -= 2.0
            points[3, residue_index] += 2.0
        end
    end

    # 3. Flanking residue SS criteria ------------------------------------------
    _points = copy(points)
    blur_amount = 0.5
    points[:, 2] += _points[:, 1] .* blur_amount
    for i in 2:(N-1)
        value = _points[:, i] .* blur_amount
        points[:, (i-1)] += value
        points[:, (i+1)] += value
    end
    points[:, (N-1)] += _points[:, N] .* blur_amount

    # 4. Terminals criteria ----------------------------------------------------
    points[3, 1] += 5.0
    points[3, N] += 5.0

    max_points = map((x) -> trunc(Int, ceil(x[1])), argmax(points, dims=1))
    converted_points = map((x) -> convert(Float64,x), max_points[1, :])
    return Base.join(map((x) -> ss_assignment[x], converted_points))
end

catergorize_ss_from_dihedral_angles(pose::Pose) = begin
    catergorize_ss_from_dihedral_angles(pose, nothing)
end