using ProtoSyn
using ProtoSyn.Units
using CurveFit: Polynomial

"""
    setss!(state::State, ss::SecondaryStructureTemplate{T}, residues::Vector{Residue}; [include_variation::Bool = false], [min_prob::T = 0.0]) where {T <: AbstractFloat}

Set the phi `ϕ`, psi `ψ` and omega `ω` backbone angles of all [`Residue`](@ref)
instances in the given `residues` vector to match the provided
[`SecondaryStructureTemplate`](@ref). This function acts on the internal
coordinates and does not update cartesian coordinates (using the
[`sync!`](@ref ProtoSyn.sync!) method), although a request for conversion is
made (by calling the [`request_i2c!`](@ref ProtoSyn.request_i2c!)). If
`include_variation` is set to `true` (`false, by default`), the phi `ϕ`, psi `ψ`
and omega `ω` backbone angles are sampled from each distribution, instead of
using the ideal angle. `min_prob` defines the minimum probability of the sampled
angle (from 0.0 to 1.0). Using higher `min_prob` values results in variations
conformationally closer to the ideal dihedral angle values.

    setss!(pose::Pose, ss::SecondaryStructureTemplate; [include_variation::Bool = false], [min_prob::T = 0.0])
    setss!(pose::Pose, ss::SecondaryStructureTemplate, sele::ProtoSyn.AbstractSelection; [include_variation::Bool = false], [min_prob::T = 0.0])
    setss!(pose::Pose, ss::SecondaryStructureTemplate, residue::Residue; [include_variation::Bool = false], [min_prob::T = 0.0])

Set the phi `ϕ`, psi `ψ` and omega `ω` backbone angles of all [`Residue`](@ref)
instances in the given [`Pose`](@ref) `pose` to match the provided
[`SecondaryStructureTemplate`](@ref). If an optional `AbstractSelection` `sele`
is provided, apply the transformation only to the selected [`Residue`](@ref)
instances. Optionally, a single [`Residue`](@ref) instance can also be provided:
any changes will only apply to the selected [`Residue`](@ref) instance.

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
            ProtoSyn.setdihedral!(state, Peptides.phi(r), T(-75°))
            # Last residues of chain might not have a psi angle.
            Peptides.psi(r) !== nothing && ProtoSyn.setdihedral!(state, Peptides.psi(r),  T(145°))
            ProtoSyn.setdihedral!(state, Peptides.omega(r), T(180°))
            continue
        end
        if !include_variation
            ProtoSyn.setdihedral!(state, Peptides.phi(r), ss.ϕ.angle)
            # Last residues of chain might not have a psi angle.
            Peptides.psi(r) !== nothing && ProtoSyn.setdihedral!(state, Peptides.psi(r),  ss.ψ.angle)
            ProtoSyn.setdihedral!(state, Peptides.omega(r), ss.ω.angle)
        else
            _ϕ = ss.ϕ.sampler === nothing ? ss.ϕ.angle : Peptides.sample_ramachandran(ss.ϕ.sampler, min_prob = min_prob)
            ProtoSyn.setdihedral!(state, Peptides.phi(r), _ϕ)
            # Last residues of chain might not have a psi angle.
            _ψ = ss.ψ.sampler === nothing ? ss.ψ.angle : Peptides.sample_ramachandran(ss.ψ.sampler, min_prob = min_prob)
            Peptides.psi(r) !== nothing && ProtoSyn.setdihedral!(state, Peptides.psi(r), _ψ)
            _ω = ss.ω.sampler === nothing ? ss.ω.angle : Peptides.sample_ramachandran(ss.ω.sampler, min_prob = min_prob)
            ProtoSyn.setdihedral!(state, Peptides.omega(r), _ω)
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
    read_ss_map(pose::Pose, filename::String, ss_type::String)

Reads a [`SecondaryStructure`](@ref) map from `filename` into an
`AbstractSelection` of [`Residue`](@ref) instances marked with the desired
`ss_type` in the map. The map is checked against the given [`Pose`](@ref) for
compatibility (number of [`Residue`](@ref) instances, etc). The expected map
format follows the DeepConCNF_SS3 format (3-mode categorization), as given by
[RaptorX prediction server](http://raptorx.uchicago.edu/StructurePropertyPred/predict/).

# Examples
```
julia> ProtoSyn.Peptides.read_ss_map("ss_map.txt", "H")
 (...)
```
"""
function read_ss_map(pose::Pose, filename::String, ss_type::String)

    sequence = ProtoSyn.sequence(pose)
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
    categorize_ss_from_dihedral_angles(pose::Pose, [selection::Opt{ProtoSyn.AbstractSelection} = nothing]; [blur_amount::T = 0.5]) where {T <: AbstractFloat}

This method attempts to categorize the given [`Pose`](@ref) `pose`
[`Residue`](@ref) instances in a 3-mode categorization format: "H" for helix,
"E" for beta sheets and "C" for coils. The categorization takes into account:

 * Geometrical criteria: accordance of the current phi `ϕ`, psi `ψ` and omega `ω` backbone angles to known potentials. Make sure the [`Pose`](@ref) `pose` is synched (see [`sync!`](@ref ProtoSyn.sync!));
 * Hydrogen bonding pattern: counts hydrogen bonds using [`generate_hydrogen_bond_network`](@ref ProtoSyn.Calculators.HydrogenBonds.generate_hydrogen_bond_network). Make sure the [`Pose`](@ref) `pose` has charges (see [`assign_default_charges!`](@ref ProtoSyn.Peptides.Calculators.Electrostatics.assign_default_charges!));
 * Flanking [`Residue`](@ref) secondary structure: flanking [`Residue`](@ref)'s secondary structure "spills" over and influences neighboring [`Residue`](@ref) instances. The influence amount can be set with `blur_amount` (0.5, by default);
 * Terminal [`Residue`](@ref) instances: Terminal [`Residue`](@ref) instances tend to adopt coil conformations. 

!!! ukw "Note:"
    In ProtoSyn 1.1, this method offers notoriously sub-par predictions. Using external tools such as DSSP or RaptorX secondary structure prediction server is recommended.

# Examples
```
julia> ProtoSyn.Peptides.categorize_ss_from_dihedral_angles(pose)
"CHHHHHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHHHHHHHCCCHHHHHHHHHHHHHHHHHHHHHHHHC"
```
"""
function categorize_ss_from_dihedral_angles(pose::Pose, selection::Opt{ProtoSyn.AbstractSelection} = nothing; blur_amount::T = 0.5) where {T <: AbstractFloat}

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

    phi_potentials = Vector{Polynomial}([
        ProtoSyn.Peptides.phi_α_L_potential,
        ProtoSyn.Peptides.phi_α_R_potential,
        ProtoSyn.Peptides.phi_β_potential,
        ProtoSyn.Peptides.phi_coil_potential])

    phi_keys   = [  1,   1,   2,   3]
    phi_values = [0.0, 0.0, 0.0, 0.0]

    psi_potentials = Vector{Polynomial}([
        ProtoSyn.Peptides.psi_α_L_potential,
        ProtoSyn.Peptides.psi_α_R_potential,
        ProtoSyn.Peptides.psi_β_potential,
        ProtoSyn.Peptides.psi_coil_potential])

    psi_keys   = [  1,   1,   2,   3]
    psi_values = [0.0, 0.0, 0.0, 0.0]

    for (residue_index, residue) in enumerate(residues)
        phi_atom = Peptides.phi(residue)
        phi_atom === nothing && continue
        phi = ProtoSyn.getdihedral(pose.state, phi_atom)
        
        for (ss_index, potential) in enumerate(phi_potentials)
            value = potential(phi)
            phi_values[ss_index] = value
        end
        @info "Phi:\n values = $phi_values\n keys = $phi_keys\n Min value: $(argmin(phi_values))\n Chosen key: $(phi_keys[argmin(phi_values)])"
        phi_ss = phi_keys[argmin(phi_values)]
        
        psi_atom = Peptides.psi(residue)
        psi_atom === nothing && continue
        psi = ProtoSyn.getdihedral(pose.state, psi_atom)
        for (ss_index, potential) in enumerate(psi_potentials)
            value = potential(psi)
            psi_values[ss_index] = value
        end
        @info "Psi:\n values = $psi_values\n keys = $psi_keys\n Min value: $(argmin(psi_values))\n Chosen key: $(psi_keys[argmin(psi_values)])"
        psi_ss = psi_keys[argmin(psi_values)]

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
        pose, nothing, false; hydrogen_bond_network = hbn,
        potential = ProtoSyn.Calculators.get_bump_potential_charges(c = 3.0, r = 1.5))

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
            points[1, residue_index] -= 12.0
            points[2, residue_index] -= 12.0
            points[3, residue_index] += 12.0
        end
    end

    # 3. Flanking residue SS criteria ------------------------------------------
    _points     = copy(points)
    points[:, 2] += _points[:, 1] .* blur_amount
    for i in 2:(N-1)
        value = _points[:, i] .* blur_amount
        points[:, (i-1)] += value
        points[:, (i+1)] += value
    end
    points[:, (N-1)] += _points[:, N] .* blur_amount

    # 4. Terminals criteria ----------------------------------------------------
    points[3, 1] += 12.0
    points[3, N] += 12.0

    max_points = map((x) -> trunc(Int, ceil(x[1])), argmax(points, dims=1))
    converted_points = map((x) -> convert(Float64,x), max_points[1, :])
    return Base.join(map((x) -> ss_assignment[x], converted_points))
end