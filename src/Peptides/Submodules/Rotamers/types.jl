using ProtoSyn: State, setdihedral!, Residue
using Random
using Printf

abstract type RotamerLibrary end

"""
    Rotamer{T <: AbstractFloat}(name::String, chis::Dict{AbstractSelection, Tuple{T, T}})

A [`Rotamer`](@ref) holds information regarding a single conformation for all
chi dihedral angles of a sidechain belonging to an aminoacid identified by the 
given `name`. The `chis` list is, therefore, a dictionary, where the key is the
`AbstractSelection` (chi1, chi2, chi3 or chi4) and the value is a `Tuple{T, T}`,
where the first entry is the average dihedral angle and the second entry is the
standard deviation expected for that dihedral angle.
    
# See also
[`get_rotamer`](@ref) [`apply!`](@ref)

# Examples
```
julia> rot_lib["LYS"][35°, -35°][1]
Rotamer{Float64}: LYS | Chi1:   -67.5° ±  6.9 | Chi2:  -179.6° ±  9.7 | Chi3:  -179.7° ± 11.8 | Chi4:   178.5° ± 12.0
```
"""
mutable struct Rotamer{T <: AbstractFloat}
    name::String
    chis::Dict{AbstractSelection, Tuple{Opt{T}, T}}
end

function Base.show(io::IO, r::Rotamer{T}) where {T <: AbstractFloat}
    print(io, "Rotamer{$T}: $(as_string(r))\n")
end

function Base.getindex(r::Rotamer{T}, i::AbstractSelection) where {T <: AbstractFloat}
    k = collect(keys(r.chis))
    return r.chis[k[findfirst((c) -> c.n === i.n, k)]]
end

function Base.setindex!(r::Rotamer{T}, v::Tuple{T, T}, i::AbstractSelection) where {T <: AbstractFloat}
    k = collect(keys(r.chis))
    r.chis[k[findfirst((c) -> c.n === i.n, k)]] = v
end

function as_string(r::Rotamer{T}) where {T <: AbstractFloat}
    chis = ""
    for index in 1:4
        if index <= length(r.chis)
            name = "chi$index"
            value, sd = r[ChiSelection(index)]
            if value === nothing
                chis *= @sprintf " | %s: %11s ± %4.1f" titlecase(name) "Not defined" rad2deg(sd)
            else
                chis *= @sprintf " | %s: %7.1f° ± %4.1f" titlecase(name) rad2deg(value) rad2deg(sd)
            end
        else
            name = "chi$index"
            chis *= @sprintf " | %s: %-15s" titlecase(name) "--"
        end
    end
    return "$(r.name)$chis"
end

# ---

"""
    BBI_RotamerLibrary{T <: AbstractFloat}(rotamers::Vector{Rotamer{T}}, weights::Weights{T, T, Array{T, 1}})

A [`BBI_RotamerLibrary`](@ref) is a backbone-independent list of
[`Rotamer`](@ref) instances, ordered based on the natural probability of
occurrence, where each position in the `rotamers` list has a corresponding
position in the `weights` list.

# See also
[`sample_rotamer`](@ref)

# Examples
```
julia> rot_lib["VAL"][35°, -35°]

+---------------------------------------------------------------------------------------------------------------------------+
| Index | Probability | Rotamer description                                                                                 |
+---------------------------------------------------------------------------------------------------------------------------+
| 1     |      82.07% | VAL | Chi1:   175.8° ±  5.6 | Chi2: --              | Chi3: --              | Chi4: --              |
| 2     |      13.74% | VAL | Chi1:    64.7° ±  7.5 | Chi2: --              | Chi3: --              | Chi4: --              |
| 3     |       4.18% | VAL | Chi1:   -61.6° ±  6.4 | Chi2: --              | Chi3: --              | Chi4: --              |
+---------------------------------------------------------------------------------------------------------------------------+
```
"""
mutable struct BBI_RotamerLibrary{T <: AbstractFloat} <: RotamerLibrary
    rotamers::Vector{Rotamer{T}}
    weights::Weights{T, T, Array{T, 1}}
end

BBI_RotamerLibrary(::Type{T}) where {T <: AbstractFloat} = begin
    BBI_RotamerLibrary(Vector{Rotamer{T}}(), Weights(Vector{T}()))
end

BBI_RotamerLibrary() = BBI_RotamerLibrary(ProtoSyn.Units.defaultFloat)

function Base.show(io::IO, r::BBI_RotamerLibrary)
    println(io, "\n+"*repeat("-", 123)*"+")
    @printf(io, "| %-5s | %10s | %-99s |\n", "Index", "Probability", "Rotamer description")
    println(io, "+"*repeat("-", 123)*"+")
    for (index, rotamer) in enumerate(r.rotamers)
        @printf(io, "| %-5d | %10.2f%% | %-80s |\n", index, r.weights[index]*100, as_string(rotamer))
    end
    println(io, "+"*repeat("-", 123)*"+")
end

# ---

"""
    BBD_RotamerLibrary{T <: AbstractFloat}(name::String, phis::Vector{T}, psis::Vector{T}, rotamer_stacks::Matrix{BBI_RotamerLibrary}) where {T <: AbstractFloat}

A [`BBD_RotamerLibrary`](@ref) is a 2D backbone dependent matrix of
[`BBI_RotamerLibrary`](@ref) instances (`rotamer_stacks`), indexed by both the list of
backbone phi dihedrals (`phis`) and the list of backbone psi dihedrals (`psis`).
Each [`BBD_RotamerLibrary`](@ref) only hold information regarding one type of
aminoacid, identified by the `name`. Each entry in the `rotamer_stacks` matrix
only holds information of the list of rotamers available/plausible for a
specific set of phi and psi backbone dihedrals.

# See also
[`load_dunbrack`](@ref)

# Examples
```
julia> rot_lib["VAL"]
Name: VAL | Shape: (37, 37)
```
"""
mutable struct BBD_RotamerLibrary{T <: AbstractFloat} <: RotamerLibrary
    name::String
    phis::Vector{T}
    psis::Vector{T}
    rotamer_stacks::Matrix{BBI_RotamerLibrary}
end

function Base.show(io::IO, r::BBD_RotamerLibrary)
    println(io, "Name: $(r.name) | Shape: $(size(r.rotamer_stacks))")
end

Base.copy(rot_lib::BBD_RotamerLibrary{T}) where {T <: AbstractFloat} = begin
    return BBD_RotamerLibrary(
        rot_lib.name,
        copy(rot_lib.phis),
        copy(rot_lib.psis),
        copy(rot_lib.rotamer_stacks)
    )
end
