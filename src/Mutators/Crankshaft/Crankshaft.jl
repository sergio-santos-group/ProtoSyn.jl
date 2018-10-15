module Crankshaft

using ..Drivers
using ..Common
using ..Aux

@doc raw"""
    ConfigParameters(p_mut::Float64 = 0.1)

Define the runtime parameters for Crankshaft movements.

# Arguments
- `p_mut::Float64`: Probability of mutation on this pair of residues (Default: 0.1).

# Examples
```julia-repl
julia> Mutators.Crankshaft.ConfigParameters(0.2)
ConfigParameters(p_mut=0.2)

julia> Mutators.Crankshaft.ConfigParameters()
ConfigParameters(p_mut=0.1)
```
"""
mutable struct ConfigParameters

    p_mut::Float64
    step_size::Float64

    ConfigParameters(; p_mut::Float64 = 0.1, step_size::Float64 = 0.1) = new(p_mut, step_size)
end
Base.show(io::IO, b::ConfigParameters) = print(io, "ConfigParameters(p_mut=$(b.p_mut), step_size=$(b.step_size))")

# ----------------------------------------------------------------------------------------------------
function run!(
    state::Common.State,
    dihedrals::Array{Common.Dihedral, 1},
    params::ConfigParameters,
    angle_sampler::Function;
    ostream::IO = stdout)

    #Choose dihedrals involved in the crankshaft movement
    l = length(dihedrals)
    crank_count::Int64 = 0
    for i in 1:(l-1)
        for j in (i + 1):l
            if rand() < params.p_mut
                #Find angle
                angle::Float64 = angle_sampler()

                #Rotate crankshaft
                rotate_crankshaft!(state.xyz, dihedrals[i], dihedrals[j], dihedrals[j - 1], angle)
                crank_count += 1
            end
        end
    end
    write(ostream, "(CS) Performed $crank_count crankshaft movements (step_size: $(params.step_size)).\n")
end


# function rotate_crankshaft!(
#     xyz::Array{Float64, 2},
#     dihedral_1::Common.Dihedral,
#     dihedral_2::Common.Dihedral,
#     dihedral_2_prev::Common.Dihedral,
#     angle::Float64)

#     #Verify that both dihedrals are PHI
#     if dihedral_1.dtype != "PHI" || dihedral_2.dtype != "PHI"
#         error("Both dihedrals involved in crankshaft movements need to be of dtype 'PHI'")
#     end

#     #Get carbon alphas and define rotation axis
#     pivot = xyz[dihedral_1.a3, :]'
#     axis = xyz[dihedral_2.a3, :]' - pivot

#     #Define the rotation matrix based on the rotation axis and angle
#     rmat = Aux.rotation_matrix_from_axis_angle(axis', angle)
    
#     #Disconect
#     dihedral_2_prev.residue.next = nothing
    
#     #Rotate movable atoms pertaining to this dihedral
#     xyz[dihedral_1.movable, :] = (rmat * (xyz[dihedral_1.movable, :] .- pivot)')' .+ pivot

#     #Rotate all downstream residues
#     residue = dihedral_1.residue
#     while residue.next != nothing
#         residue = residue.next
#         xyz[residue.atoms, :] = (rmat * (xyz[residue.atoms, :] .- pivot)')' .+ pivot
#     end

#     #Rotate first 3 atoms in the last residue
#     f3 = setdiff(dihedral_2.residue.atoms, dihedral_2.movable)
#     xyz[f3, :] = (rmat * (xyz[f3, :] .- pivot)')' .+ pivot

#     #Reconnect
#     dihedral_2_prev.residue.next = dihedral_2.residue
# end




function rotate_crankshaft!(
    xyz::Array{Float64, 2},
    dihedral1::Dihedral,
    dihedral2::Dihedral,
    angle::Float64)

    next = dihedral2.residue.next 
    dihedral2.residue.next = nothing
    rotate_dihedral!(xyz,
        dihedral1.a3, dihedral2.a3,
        dihedral1.dtype, angle,
        dihedral1.movable,
        dihedral1.residue)
    dihedral2.residue.next = next
    
    #movable = setdiff(dihedral2.residue.atoms, dihedral2.movable)
    rotate_dihedral!(xyz,
        dihedral1.a3, dihedral2.a3,
        dihedral1.dtype, -angle,
        dihedral2.movable)
    
end


end