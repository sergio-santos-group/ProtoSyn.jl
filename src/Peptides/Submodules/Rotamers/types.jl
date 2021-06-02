using ProtoSyn: State, setdihedral!, Residue
using ProtoSyn.Peptides.Dihedral: DihedralType
using Random
using StatsBase
using Printf

abstract type RotamerLibrary end

"""
    Rotamer{T <: AbstractFloat}(name::String, chis::Dict{DihedralType, Tuple{T, T}})

A [`Rotamer`](@ref) holds information regarding a single conformation for all
chi dihedral angles of a sidechain belonging to an aminoacid identified by the 
given `name`. The `chis` list is, therefore, a dictionary, where the key is the
`DihedralType` (chi1, chi2, chi3 or chi4) and the value is a `Tuple{T, T}`,
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
    chis::Dict{DihedralType, Tuple{T, T}}
end

function Base.show(io::IO, r::Rotamer{T}) where {T <: AbstractFloat}
    print(io, "Rotamer{$T}: $(as_string(r))\n")
end

function as_string(r::Rotamer{T}) where {T <: AbstractFloat}
    chis = ""
    for index in 1:4
        if index <= length(r.chis)
            name = "chi$index"
            value, sd = r.chis[getfield(ProtoSyn.Peptides.Dihedral, Symbol(name))]
            chis *= @sprintf " | %s: %7.1f° ± %4.1f" titlecase(name) rad2deg(value) rad2deg(sd)
        else
            name = "chi$index"
            chis *= @sprintf " | %s: %-15s" titlecase(name) "--"
        end
    end
    return "$(r.name)$chis"
end

# ---

"""
    RotamerStack{T <: AbstractFloat}(rotamers::Vector{Rotamer{T}}, weights::Weights{T, T, Array{T, 1}})

A [`RotamerStack`](@ref) is a smart list of [`Rotamer`](@ref) instances, ordered
based on the natural probability of occurrence, where each position in the
`rotamers` list has a corresponding position in the `weights` list.

# See also
[`sample`](@ref)

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
mutable struct RotamerStack{T <: AbstractFloat}
    rotamers::Vector{Rotamer{T}}
    weights::Weights{T, T, Array{T, 1}}
end

RotamerStack(::Type{T}) where {T <: AbstractFloat} = begin
    RotamerStack(Vector{Rotamer{T}}(), Weights(Vector{T}()))
end

RotamerStack() = RotamerStack(ProtoSyn.Units.defaultFloat)

function Base.show(io::IO, r::RotamerStack)
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
    BBD_RotamerLibrary{T <: AbstractFloat}(name::String, phis::Vector{T}, psis::Vector{T}, rotamer_stacks::Matrix{RotamerStack}) where {T <: AbstractFloat}

A [`BBD_RotamerLibrary`](@ref) is a 2D backbone dependent matrix of
[`RotamerStack`](@ref) instances (`rotamer_stacks`), indexed by both the list of
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
    rotamer_stacks::Matrix{RotamerStack}
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