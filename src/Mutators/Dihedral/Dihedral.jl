module Dihedral

using ..Common
using ..Aux


struct ConfigParameters

    pmut::Float64

end


mutable struct New_Dihedral
    a1::Int64
    a2::Int64
    a3::Int64
    a4::Int64
    movable::Array{Int64, 1}
    residue::Union{Common.Residue, Int64}
    dtype::String
end


function load_topology(p::Dict{String, Any})

    dihedrals = Array{New_Dihedral, 1}()
    residues = Array{Common.Residue, 1}()

    for content in p["dihedrals"]
        push!(dihedrals, New_Dihedral(
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


function run!(
    state::Common.State,
    dihedrals::Array{New_Dihedral, 1},
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


function rotate_dihedral!(
    xyz::Array{Float64, 2},
    dihedral::New_Dihedral,
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