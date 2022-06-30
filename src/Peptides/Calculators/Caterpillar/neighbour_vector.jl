"""
    ProtoSyn.Peptides.Calculators.Caterpillar.neighbour_vector([::A], pose::Pose, selection::Opt{AbstractSelection}, [update_forces::Bool = false]; [identification_curve::Function = null_identification_curve], [hydrophobicity_weight::Function = null_hydrophobicity_weight], [rmax::T = 9.0], [sc::T = 1.0], [Ω::Union{Int, T} = 750.0], [hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity]) where {A, T <: AbstractFloat}
    
Calculate the given [`Pose`](@ref) `pose` caterpillar solvation energy using the
Neighbour Vector (NV) algorithm (see [this article](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0020853)).
If an `AbstractSelection` `selection` is provided, consider only the selected
[`Atom`](@ref) instances (any given `selection` will be promoted to be of
[`Atom`](@ref) type, see [`ProtoSyn.promote`](@ref)). In this model, vectors
`ωi` are defined between an [`Atom`](@ref) and all neighbouring [`Atom`](@ref)
instances (within the defined `rmax` cut-off (in Angstrom Å)). Note that only
the selected atoms by the `selection` are considered. A resulting vector `Ωi` is
calculated by suming all `ωi` vectors, multiplied by a `w1` weight, provided by
the `identification_curve` `Function`. This `Function` receives the distance
between each pair of neighbouring atoms (as a float), the `rmax` value and,
optionally, a slope control `sc` value, and return a weight `w1` (as a float).
The `identification_curve` signature is as follows:

```
identification_curve(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}
```

In order to use pre-defined `identification_curve` `Function` instances defined
in ProtoSyn, check [`linear`](@ref), [`sigmoid`](@ref) and
[`sigmoid_normalized`](@ref).

Note that:
 - Shorter `rmax` value identify buried residues in the local environment (i.e.:
 in the scale of the secondary structure, recommended) while a larger `rmax`
 value identifies buried residues in the global scale (i.e.: in comparison with
 the whole structure).
 - The slope control `sc` value only has effect in sigmoid
 `identification_curve` `Function` instances. A smaller value augments the
 prevalence of distance information in the `w1` weight calculation, while a
 larger value defines a more strict cut-off (recommended);

An [`Atom`](@ref) is, therefore, considered buried if the magnitude of the
resulting vector from the sum of the `ωi` (multiplied by `w1`) is within a
defined cut-off value `Ω`. Buried hydrophobic aminoacids receive an energetic
reward, while exposed hydrophobic [`Residue`](@ref) instances receive a penalty
(and vice-versa for hydrophylic aminoacids), defined in the provided
`hydrophobicity_map` (hydrophobicity map examples can be found in
`Peptides.constants.jl`) and multiplied by `w2`, calculated by the
`hydrophobicity_weight` `Function`. This `Function` receives the vector
magnitude `Ωi`, the `hydrophobicity_map_value` and the cut-off value `Ω`,
returning a `w2` (as a float).  The `hydrophobicity_weight` signature is as
follows:

```
hydrophobicity_weight(Ωi::Union{Int, T}; hydrophobicity_map_value::T = 0.0, Ω::Union{Int, T} = 0.0) where {T <: AbstractFloat}
```

In order to use pre-defined `hydrophobicity_weight` `Function` instances defined
in ProtoSyn, check [`nv_scalling_exposed_only`](@ref),
[`nv_non_scalling_exposed_only`](@ref), [`nv_scalling_all_contributions`](@ref)
(recommended) and [`nv_non_scalling_all_contributions`](@ref).

The optional `A` parameter defines the acceleration
mode used (SISD_0, SIMD_1 or CUDA_2). If left undefined the default
`ProtoSyn.acceleration.active` mode will be used. This function does not
calculate forces (not applicable), and therefore the `update_forces` flag serves
solely for uniformization with other energy-calculating functions.

# See also
[`neighbour_count`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.Caterpillar.neighbour_vector(pose, false)
(0.0, nothing)
```
"""
function neighbour_vector(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; identification_curve::Function = null_identification_curve, hydrophobicity_weight::Function = null_hydrophobicity_weight, rmax::T = 9.0, sc::T = 1.0, Ω::Union{Int, T} = 750.0, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
    
    if selection === nothing
        selection = TrueSelection{Atom}()
    else
        selection = ProtoSyn.promote(selection, Atom)
    end

    dm    = ProtoSyn.Calculators.full_distance_matrix(A, pose, selection)
    if A === ProtoSyn.CUDA_2
        dm = collect(dm)
    end

    atoms = selection(pose, gather = true)
    
    # if length(atoms) != length(eachresidue(pose.graph))
    #     @warn "The number of selected residues doesn't match the number of residues in the pose ($(length(atoms)) ≠ $(length(eachresidue(pose.graph))))"
    #     return 0.0, nothing
    # end
    
    Ωis   = Vector{T}() # !
    esols = Vector{T}() # !
    esol = T(0.0)
    for i in 1:size(dm)[1]
        atom_i = atoms[i]
        xi = pose.state[atom_i].t
        Ωi = zeros(T, 3)
        for j in 1:size(dm)[2]
            i === j && continue
            w = identification_curve(dm[i, j]; rmax = rmax, sc = sc)
            v  = (pose.state[atoms[j]].t .- xi) .* w
            Ωi += v
        end

        DHI = hydrophobicity_map[atom_i.container.name]
        _Ωi = norm(Ωi)
        esol_i = hydrophobicity_weight(_Ωi, hydrophobicity_map_value = DHI, Ω = Ω)
        esol += esol_i

        push!(Ωis, _Ωi) # !
        push!(esols, esol_i) # !
    end

    return esol, nothing, Ωis, esols
    # return esol, nothing
end


neighbour_vector(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; identification_curve::Function = null_identification_curve, hydrophobicity_weight::Function = null_hydrophobicity_weight, rmax::T = 9.0, sc::T = 1.0, Ω::Union{Int, T} = 750.0, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity, kwargs...) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat} = begin
    neighbour_vector(ProtoSyn.acceleration.active, pose, selection, update_forces,
        identification_curve = identification_curve,
        hydrophobicity_weight = hydrophobicity_weight,
        rmax = rmax,
        sc = sc, Ω = Ω,
        hydrophobicity_map = hydrophobicity_map)
end