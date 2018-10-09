module Dihedral

using ..Common
using ..Aux

@doc raw"""
    ConfigParameters(p_mut::Float64 = 0.1)

Define the runtime parameters for Dihedral movements.

# Arguments
- `p_mut::Float64`: Probability of mutation of this dihedral (Default: 0.1).

# Examples
```julia-repl
julia> Mutators.Diehdral.ConfigParameters(0.2)
ConfigParameters(p_mut=0.2)

julia> Mutators.Diehdral.ConfigParameters()
ConfigParameters(p_mut=0.1)
```
"""
struct ConfigParameters

    pmut::Float64

    ConfigParameters(; p_mut::Float64 = 0.1) = new(p_mut)
end
Base.show(io::IO, b::ConfigParameters) = print(io, "ConfigParameters(p_mut=$(b.p_mut))")

# ----------------------------------------------------------------------------------------------------------

@doc raw"""
    NewDihedral(a1::Int64, a2::Int64, a3::Int64, a4::Int64, movable::Array{Int64, 1}, residue::Union{Common.Residue, Int64}, dtype::String)

Define a dihedral.

# Arguments
- `a1::Int64, a2::Int64, a3::Int64, a4::Int64`: *global* atom indices.
- `movable::Array{Int64, 1}`: List of *global* atom indices that will be moved during the dihedral movement in *this* residue.
- `residue::Union{Common.Residue, Int64}`: Residue that this dihedral belongs to. Should be a [`Common.Residue`](@ref) object.
- `dtype::String`: Dihedral type (i.e. "PHI", "PSI", ...)

# Examples
```julia-repl
julia> Mutators.Diehdral.NewDihedral(2, 3, 5, 7, [5, 6], Common.Residue([1, 2, 3, 4, 5, 6], (...), "A"), "PSI")
Dihedral(a1=2, a2=3, a3=5, a4=7, movable=[5, 6], residue=Common.Residue(atoms=[1, 2, 3, 4, 5, 6], next=V, name=A), type="PSI")
```
See also: [`load_topology`](@ref)
"""
mutable struct NewDihedral
    a1::Int64
    a2::Int64
    a3::Int64
    a4::Int64
    movable::Array{Int64, 1}
    residue::Union{Common.Residue, Int64}
    dtype::String
end
Base.show(io::IO, b::NewDihedral) = print(io, "Dihedral(a1=$(b.a1), a2=$(b.a2), a3=$(b.a3), a4=$(b.a4), movable=$(b.movable), residue=$(b.residue), type=$(b.dtype))")

# ------------------------------------------------------------------------------------------------------------

@doc raw"""
    load_topology(p::Dict{String, Any})

Parse a dictionary containing the dihedral and residue topology. Return a [`NewDihedral`](@ref) array
and a [`Common.Residue`](@ref) array.

# Examples
```julia-repl
julia> Mutators.Diehdral.load_topology(p)
(ProtoSyn.Mutators.Dihedral.NewDihedral[...], ProtoSyn.Common.Residue[...])
```
See also: [`Aux.read_JSON`](@ref)
"""
function load_topology(p::Dict{String, Any})

    dihedrals = Array{NewDihedral, 1}()
    residues = Array{Common.Residue, 1}()

    for content in p["dihedrals"]
        push!(dihedrals, NewDihedral(
            content["a1"],
            content["a2"],
            content["a3"],
            content["a4"],
            content["movable"],
            content["parent"],
            content["type"]
        ))
    end
    
    for content in p["residues"]

        #Create residue
        new_residue = Common.Residue(
            content["atoms"],
            content["next"],
            content["type"]
        )
        
        #Set parent of this residue dihedrals
        for dihedral in dihedrals
            if dihedral.residue == content["n"]
                dihedral.residue = new_residue
            end
        end
        
        push!(residues, new_residue)
    end

    #Set correct references for dihedrals previous and next
    for residue in residues
        try
            residue.next = residues[residue.next]
        catch LoadError
            residue.next = nothing
        end
    end

    return dihedrals, residues
end

# -----------------------------------------------------------------------------------------------------------

@doc raw"""
    run!(state::Common.State, dihedrals::Array{NewDihedral, 1}, params::Config.Parameters, angle_sampler::Function[, ostream::IO = stdout])

Iterate over a list of [`NewDihedral`](@ref) (`dihedrals`) and perform dihedral movements on the current
[`Common.State`](@ref). The probability of each dihedral undergo movements is defined in the
[`ConfigParameters`](@ref) (`params.pmut`). The new angle is obtained from `angle_sampler`, who should
return a `Float64` in radians.
After movement, the [`Common.State`](@ref) is updated with the new conformation.
Any logging is written to `ostream` (Default: `stdout`).

# Examples
```julia-repl
julia> Mutators.Diehdral.run!(state, dihedrals, params, () -> randn())
```
See also: [`rotate_dihedral!`](@ref)
"""
function run!(
    state::Common.State,
    dihedrals::Array{NewDihedral, 1},
    params::ConfigParameters,
    angle_sampler::Function;
    ostream::IO = stdout)
    
    for dihedral in dihedrals
        if rand() < params.pmut
            angle::Float64 = angle_sampler()
            rotate_dihedral!(state.xyz, dihedral, angle)
        end
    end
end

# -----------------------------------------------------------------------------------------------------------

@doc raw"""
    rotate_dihedral!(xyz::Array{Float64, 2}, dihedral::NewDihedral, angle::Float64)

Perform a dihedral movement, adding the provided `angle` (in radians). If the `dihedral.dtype` is "PHI"
or "PSI" the `dihedral.residue.next` is also rotated and this is propagated recursively until the end of
the molecule. 

# Examples
```julia-repl
julia> Mutators.Diehdral.rotate_dihedral(state.xyz, dihedral, π/2)
```
See also: [`run!`](@ref) [`Aux.rotation_matrix_from_axis_angle`](@ref)
"""
function rotate_dihedral!(
    xyz::Array{Float64, 2},
    dihedral::NewDihedral,
    angle::Float64)

    pivot = xyz[dihedral.a2,:]'
    axis = xyz[dihedral.a3,:]' - pivot

    #Define the rotation matrix based on the rotation axis and angle
    rmat = Aux.rotation_matrix_from_axis_angle(axis', angle)

    #Rotate movable atoms pertaining to this dihedral
    xyz[dihedral.movable, :] = (rmat * (xyz[dihedral.movable, :] .- pivot)')' .+ pivot

    #Rotate all downstream residues
    if dihedral.dtype in ["PHI", "PSI"]
        residue = dihedral.residue
        while residue.next != nothing
            residue = residue.next
            xyz[residue.atoms, :] = (rmat * (xyz[residue.atoms, :] .- pivot)')' .+ pivot
        end
    end
end

end