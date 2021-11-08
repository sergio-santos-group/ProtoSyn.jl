"""
    ProtoSyn.Peptides.Calculators.Caterpillar.neighbour_vector([::A], pose::Pose, update_forces::Bool = false; selection::AbstractSelection = ProtoSyn.TrueSelection{Atom}(), identification_curve::Function = () -> nothing, hydrophobicity_weight::Function = () -> nothing, rmax::T = 9.0, sc::T = 1.0, Ω::Union{Int, T} = 750.0, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity) where {A, T <: AbstractFloat}
    
Calculate the given [`Pose`](@ref) `pose` caterpillar solvation energy using the
Neighbour Vector (NV) algorithm (see [this article](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0020853)).
In this model, vectors `ωi` are defined between an [`Atom`](@ref) and all
neighbouring [`Atom`](@ref) instances (within the defined `rmax` cut-off (in
Angstrom Å)). Note that only the selected atoms by the `selection` are
considered. A resulting vector `Ωi` is calculated by suming all `ωi` vectors,
multiplied by a `w1` weight, provided by the `identification_curve` `Function`.
This `Function` receives the distance between each pair of neighbouring atoms
(as a float), the `rmax` value and, optionally, a slope control `sc` value, and
return a weight `w1` (as a float). The `identification_curve` signature is as
follows:

```
linear(distance::T; rmax::T = 9.0, sc::T = 1.0) where {T <: AbstractFloat}
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
`hydrophobicity_map`. Other map examples can be found in
`Peptides.constants.jl`. The optional `A` parameter defines the acceleration
mode used (SISD_0, SIMD_1 or CUDA_2). If left undefined the default
`ProtoSyn.acceleration.active` mode will be used. This function does not
calculate forces (not applicable), and therefore the `update_forces` flag serves
solely for uniformization with other energy-calculating functions.

# See also
[`neighbour_count`](@ref)

# Examples
```
julia> ProtoSyn.Peptides.Calculators.Caterpillar.calc_solvation_energy(pose)
(925.5142248612556, nothing)
```
"""
function neighbour_vector(::Type{A}, pose::Pose, update_forces::Bool; selection::AbstractSelection = ProtoSyn.TrueSelection{Atom}(), identification_curve::Function = () -> nothing, hydrophobicity_weight::Function = () -> nothing, rmax::T = 9.0, sc::T = 1.0, Ω::Union{Int, T} = 750.0, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
    
    dm    = ProtoSyn.Calculators.full_distance_matrix(A, pose, selection)
    if A === ProtoSyn.CUDA_2
        dm = collect(dm)
    end
    atoms = selection(pose, gather = true)
    
    # Ωis   = Vector{T}()
    # esols = Vector{T}()
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
        # push!(Ωis, _Ωi)
        # push!(esols, esol_i)
    end

    # return Ωis, esols
    return esol, nothing
end

neighbour_vector(pose::Pose, update_forces::Bool; selection::AbstractSelection = ProtoSyn.TrueSelection{Atom}(), identification_curve::Function = () -> nothing, hydrophobicity_weight::Function = () -> nothing, rmax::T = 9.0, sc::T = 1.0, Ω::Union{Int, T} = 750.0, hydrophobicity_map::Dict{String, T} = ProtoSyn.Peptides.doolitle_hydrophobicity, kwargs...) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat} = begin
    neighbour_vector(ProtoSyn.acceleration.active, pose, update_forces;
        selection = selection,
        identification_curve = identification_curve,
        hydrophobicity_weight = hydrophobicity_weight,
        rmax = rmax,
        sc = sc, Ω = Ω,
        hydrophobicity_map = hydrophobicity_map)
end