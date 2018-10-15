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



using Printf

# mutable struct CallbackObject
#     freq::Int64
#     callback::Function
# end


# mutable struct MonteCarloDriver
#     sampler! :: Function
#     evaluator! :: Function
#     temperature :: Float64
#     nsteps::Int64
#     # callback::Union{Function, Nothing}
#     # callback_freq::Int64
# end


# macro cbcall(cb, step, args...)
#     ex = quote
#         if (($cb != nothing) && (getproperty($cb, :freq)>0) && ($step%getproperty($cb, :freq)==0))
#             getproperty($cb, :callback)($step, $(args...))
#         end
#     end
#     println(esc(ex))
#     return esc(ex)
# end


# function run!(state::Common.State, driver::MonteCarloDriver, callback::Union{CallbackObject, Nothing}=nothing)
    
#     step = 0
#     xyz0 = copy(state.xyz)
#     ene0 = driver.evaluator!(state, false)
#     acceptance_count = 0

#     while step < driver.nsteps
#         step += 1
#         driver.sampler!(state)
#         ene1 = driver.evaluator!(state, false)

#         if (ene1 < ene0) || (rand() < exp(-(ene1 - ene0) / driver.temperature))
#             ene0 = ene1
#             xyz0[:] = state.xyz
#             acceptance_count += 1
#         else
#             state.xyz[:] = xyz0
#         end
#         # write(stdout, @sprintf "(MC) Step: %4d | Energy: %9.4f\n" step ene1)

#         # @cbcall driver.callback driver.callback.freq step state
#         @cbcall callback step state driver acceptance_count/step

#     end
# end




# -------------------------------------
# LOAD
# -------------------------------------
dihd_json = Aux.read_JSON("data/1ctf_mc_top.json")
dihedrals, residues = Mutators.Dihedral.load_topology(dihd_json)
state = Common.load_from_pdb("data/mol.pdb")
state.energy = Common.Energy()


# -------------------------------------
# SAMPLER
# -------------------------------------
dihedral_mutator = Mutators.Dihedral.DihedralMutator(
    dihedrals   # list of dihedrals
    , 0.1       # single dihedral mutation probability
    , randn     # angle sampler
    , 1.0       # stepsize
)
my_sampler(st::Common.State) = Mutators.Dihedral.run!(st, dihedral_mutator)


# -------------------------------------
# EVALUATOR
# -------------------------------------
function my_evaluator(st::Common.State, do_forces::Bool)
    #st.energy.eTotal = rand()
    #return st.energy.eTotal
    st.energy.eTotal
end


# -------------------------------------
# DRIVER
# -------------------------------------
mc_driver = Drivers.MonteCarlo.MonteCarloDriver(
    my_sampler      # sampler
    , my_evaluator  # evaluator
    , 1.0           # temperature
    , 1000          # nsteps
)
    

# -------------------------------------
# CALLBACK
# -------------------------------------
function callback(step::Int64, st::Common.State, dr::Drivers.MonteCarlo.MonteCarloDriver, ac_ratio::Float64)
    write(stdout, @sprintf "(MC) Step: %4d | Energy: %9.4f\n" step state.energy.eTotal)
    println("$ac_ratio -> $(dihedral_mutator.stepsize)")
    dihedral_mutator.stepsize *= 0.95

end

my_callback = Common.CallbackObject(
    10          # calling frequency
    , callback  # the actual callback function
)



Drivers.MonteCarlo.run!(state, mc_driver, my_callback)








# fout = open("out/dihd.xyz", "w")

# function my_cb(step::Int64, st::Common.State, dr::MonteCarloDriver, ac_ratio::Float64)
#     write(stdout, @sprintf "(MC) Step: %4d | Energy: %9.4f\n" step state.energy.eTotal)
#     println(dr.temperature)
#     #Print.as_xyz(st, ostream=fout, title="Step $step")
# end

# cb = CallbackObject(10, my_cb)

# mc_driver = MonteCarloDriver(
#     dihedral_sampler    # sampler!
#     , (a,b) -> 1.0      # evaluator
#     , 1.0               # temperature
#     , 1000              # nsteps
#     # , cb        # callback object
#     # , my_cb     # callback
#     # , 10        # callback_freq
# )

# # mc_driver.nsteps = 1000
# # mc_driver.callback = eq_callback
# # run!(state, mc_driver)

# # mc_driver.nsteps = 10000
# # mc_driver.callback = prod_callback
# # run!(state, mc_driver)

# run!(state, mc_driver, cb)



# function do_work(p::Float64, nsteps::Int64)
#     # fout = open("out/dihd.xyz", "w")
#     #dihedral_mutator.pmut = p
#     # Print.as_xyz(state, ostream=fout, title="Step 0")
#     #for n = 1:nsteps
#     #    Mutators.Dihedral.run!(state, dihedral_mutator)
#     #    # Print.as_xyz(state, ostream=fout, title="Step $n")
#     #end
#     # close(fout)
    
#     dihedral_mutator.pmut = p
#     mc_driver.nsteps = nsteps
#     run!(state, mc_driver)
# end

# do_work(0.0, 1)

# using Profile
# @time do_work(0.01, 100)

# close(fout)