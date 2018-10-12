using ProtoSyn

# using ProtoSyn.Aux
# using ProtoSyn.Common

# @enum DIHEDRALTYPE begin
#     phi   = 0
#     psi   = 1
#     omega = 2
#     chi1  = 3
#     chi2  = 4
#     chi3  = 5
#     chi4  = 6
#     chi5  = 7
# end


# mutable struct MutableDihedral
#     a1::Int64
#     a2::Int64
#     a3::Int64
#     a4::Int64
#     movable::Vector{Int64}
#     residue::Union{Common.Residue, Int64}
#     dtype::DIHEDRALTYPE
# end


# function load_topology(p::Dict{String, Any})

#     dihedrals = Array{MutableDihedral, 1}()
#     residues = Array{Common.Residue, 1}()
    
#     str2enum = Dict(string(s) => s for s in instances(DIHEDRALTYPE))
#     residues = Dict(d["n"]=>Common.Residue(d["atoms"],d["next"],d["type"]) for d in p["residues"])
    
#     dihedrals = [
#         MutableDihedral(d["a1"], d["a2"], d["a3"], d["a4"],
#             d["movable"], residues[d["parent"]], str2enum[lowercase(d["type"])])
#         for d in p["dihedrals"]
#     ]
#     # for d in p["dihedrals"]
#     #     push!(dihedrals, MutableDihedral(
#     #         d["a1"],
#     #         d["a2"],
#     #         d["a3"],
#     #         d["a4"],
#     #         d["movable"],
#     #         residues[d["parent"]],
#     #         str2enum[lowercase(d["type"])]
#     #     ))
#     # end
    

#     # for content in p["residues"]

#     #     #Create residue
#     #     new_residue = Common.Residue(
#     #         content["atoms"],
#     #         content["next"],
#     #         content["type"]
#     #     )
        
#     #     #Set parent of this residue dihedrals
#     #     for dihedral in dihedrals
#     #         if dihedral.residue == content["n"]
#     #             dihedral.residue = new_residue
#     #         end
#     #     end
        
#     #     push!(residues, new_residue)
#     # end

#     #Set correct references for dihedrals previous and next
#     for residue in values(residues)
#         try
#             residue.next = residues[residue.next]
#         catch LoadError
#             residue.next = nothing
#         end
#     end

#     return dihedrals, residues
# end


# function rotate_dihedral!(
#     xyz::Array{Float64,2},
#     dihedral::MutableDihedral,
#     angle::Float64)

#     pivot = xyz[dihedral.a2, :]
#     axis  = xyz[dihedral.a3, :] - pivot
#     pivot = pivot'
    
#     # Define the rotation matrix based on the rotation axis and angle
#     rmat = Aux.rotation_matrix_from_axis_angle(axis, angle)

#     #Rotate movable atoms pertaining to this dihedral
#     xyz[dihedral.movable, :] = (rmat * (xyz[dihedral.movable, :] .- pivot)')' .+ pivot

    
#     # Rotate all downstream residues
#     if dihedral.dtype < omega
#         idxs = Vector{Int64}()
#         residue = dihedral.residue
#         while residue.next != nothing
#             residue = residue.next
#             append!(idxs, residue.atoms)
#         end
#         xyz[idxs, :] = (rmat * (xyz[idxs, :] .- pivot)')' .+ pivot
#     end
# end


# mutable struct DihedralMutator
#     dihedrals::Vector{MutableDihedral}
#     pmut::Float64
#     angle_sampler::Function
# end




# function run!(state::Common.State, mutator::DihedralMutator)
#     for dihedral in mutator.dihedrals
#         if rand() < mutator.pmut
#             rotate_dihedral!(state.xyz, dihedral, mutator.angle_sampler())
#         end
#     end
# end

#-------------------------------------------

dihd_json = Aux.read_JSON("data/1ctf_mc_top.json")
dihedrals, residues = Mutators.Dihedral.load_topology(dihd_json)

state = Common.load_from_pdb("data/mol.pdb")

dihedral_mutator = Mutators.Dihedral.DihedralMutator(
    dihedrals,          # list of dihedrals
    0.1,                # single dihedral mutation probability 
    () -> randn()       # angle sampler
)

function do_work(p::Float64, nsteps::Int64)
    # fout = open("out/dihd.xyz", "w")
    dihedral_mutator.pmut = p
    # Print.as_xyz(state, ostream=fout, title="Step 0")
    for n = 1:nsteps
        Mutators.Dihedral.run!(state, dihedral_mutator)
        # Print.as_xyz(state, ostream=fout, title="Step $n")
    end
    # close(fout)
end

do_work(0.0, 1)

using Profile
@time do_work(0.01, 1000)
