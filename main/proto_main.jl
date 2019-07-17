using ProtoSyn
using Printf

include("src/aux.jl")
include("src/conf.jl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   M A I N    P R O T O S Y N    A L G O R I T H M
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

state, metadata = Common.load_from_pdb(input_pdb)
Common.apply_ss!(state, metadata, ss)
nb_dhs          = filter(x -> x.residue.ss == Common.SS.COIL, metadata.dihedrals)
nb_phi_dhs      = filter(x -> x.dtype == Common.DIHEDRAL.phi, nb_dhs)

include("src/energy_components.jl")
include("src/sampling_stage.jl")
include("src/refinement_stage.jl")

include("src/check_native_str.jl")

if ""!=PROGRAM_FILE && realpath(@__FILE__) == realpath(PROGRAM_FILE)
    sampling_bests   = open("out/sampling.pdb", "w")
    refinement_bests = open("out/refinement.pdb", "w")
    log_energy       = open("out/energy.log", "w")
    debug_file       = open("out/debug_file.pdb", "w")
    
    write(log_energy, print_energy_keys("", 0)*"\n")
    check_native_str()
    
    for repeat in 1:1000
        debug_file       = open("out/debug_file.pdb", "w")
        state, metadata = Common.load_from_pdb(input_pdb)
        Common.apply_ss!(state, metadata, ss)
        printstyled("\n (CYCLE $repeat) Sampling stage:\n", color = :bold)
        @time Drivers.run!(state, smpl_ilsrr_driver)
        printstyled("\n (CYCLE $repeat) Refinement stage:\n", color = :bold)
        @time Drivers.run!(state, rfnm_sa_driver)
        Print.as_pdb(refinement_bests, state, metadata)
        flush(refinement_bests)
    end
end