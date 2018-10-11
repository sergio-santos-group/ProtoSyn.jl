using ProtoSyn

#Load state from PDB file
state = Common.load_from_pdb("data/EVK.pdb")
state.energy = Forcefield.Energy()

#Load topologies from JSON files
topol = Forcefield.load_from_json("data/EVK.json")

# function do_work(n::Int64)
#     for i = 1:n
#         fill!(state.forces, 0.0)
#         # Forcefield.evaluate!(topol.bonds, state, do_forces=true)
#         # Forcefield.evaluate!(topol.angles, state, do_forces=true)
#         Forcefield.evaluate!(topol.dihedralsCos, state, do_forces=true)
#         # Forcefield.evaluate!(topol.atoms, state, do_forces=true, cut_off=Inf)
#     end
# end

# # force JIT compilation
# do_work(2)

# @time do_work(10000)


sd_parm = Drivers.SteepestDescent.ConfigParameters(
    max_step = 1e-2,
    log_freq = 100,
    callback_freq = 500,
    f_tol = 0.1,
)

x = zeros(size(state.xyz))
x[:] = state.xyz[:]

function my_evaluator!(st::Common.State, do_forces::Bool)
    # if do_forces
    #     fill!(st.forces, 0.0)
    # end
    # state.xyz[:] = x[:]
    #println(st.xyz)
    energy = 0.0
    energy += Forcefield.evaluate!(topol.bonds,        st, do_forces=do_forces)
    energy += Forcefield.evaluate!(topol.angles,       st, do_forces=do_forces)
    energy += Forcefield.evaluate!(topol.dihedralsCos, st, do_forces=do_forces)
    energy += Forcefield.evaluate!(topol.atoms,        st, do_forces=do_forces, cut_off=Inf)
    # println(st.forces)
    # println(energy)
    return energy
end


function do_work(n::Int64)
    fout = open("out/teste_sd.xyz", "w")
    sd_parm.n_steps = n
    Drivers.SteepestDescent.run!(
        state,
        my_evaluator!,
        sd_parm,
        #callback = (st::Common.State, step::Int) -> Print.as_xyz(st, ostream=fout, title="Step $step")
    )
    close(fout)
end

# state.xyz *= 2.0
do_work(0)
state.xyz[:] = x[:]
@time do_work(1000000)