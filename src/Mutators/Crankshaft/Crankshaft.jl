module Crankshaft

using ..Common
using ..Aux

# ----------------------------------------------------------------------------------------------------

#TODO: Document structure
mutable struct CrankshaftMutator
    dihedrals::Vector{Common.Dihedral}
    pmut::Float64
    angle_sampler::Function
    stepsize::Float64
end

# ----------------------------------------------------------------------------------------------------

#TODO: Document function
@inline function run!(state::Common.State, mutator::CrankshaftMutator)

    l = length(dihedrals)
    for i in 1:(l-1)
        for j in (i + 1):l
            if rand() < params.p_mut
                rotate_crankshaft!(state.xyz, dihedrals[i], dihedrals[j], mutator.angle_sampler())
            end
        end
    end
end

# ----------------------------------------------------------------------------------------------------

#TODO: Document function
function rotate_crankshaft!(
    xyz::Array{Float64, 2},
    dihedral1::Common.Dihedral,
    dihedral2::Common.Dihedral,
    angle::Float64)

    next = dihedral2.residue.next 
    dihedral2.residue.next = nothing
    Common.rotate_dihedral!(xyz, dihedral1.a3, dihedral2.a3, dihedral1.dtype, angle, dihedral1.movable, dihedral1.residue)
    dihedral2.residue.next = next
    
    Common.rotate_dihedral!(xyz, dihedral1.a3, dihedral2.a3, dihedral1.dtype, -angle, dihedral2.movable)
end

end