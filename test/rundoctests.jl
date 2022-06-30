# * This testing script evaluates all `jldoctest` entries in the ProtoSyn
# * documentation, testing to see if the displayed output matches the expected
# * return value. The necessary objects and variables are set using
# * `DocMeta.setdocmeta!`. In order to avoid errors, run this script after the
# * regular test script `runtests.jl`. 

using Test, Documenter, ProtoSyn
ProtoSyn.set_logger_to_error()

DocMeta.setdocmeta!(ProtoSyn, :DocTestSetup, :(begin
    using ProtoSyn;
    using ProtoSyn.Units;
    using Random; Random.seed!(1);
    res_lib  = ProtoSyn.Peptides.grammar;
    frag     = fragment(res_lib, seq"AAA");
    pose     = build(res_lib, seq"SESEAEFKQRLAAIKTRLQAL"); sync!(pose);
    pose_mod = copy(pose);
    ProtoSyn.setdihedral!(pose_mod.state, pose_mod.graph[1][2]["N"], deg2rad(90));
    sync!(pose_mod);
    rrbm = ProtoSyn.Mutators.RotationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 0.4, an"CA" & an"CB")
    trbm = ProtoSyn.Mutators.TranslationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, 1.0, nothing)
    vl = ProtoSyn.Calculators.VerletList(pose.state.size); vl.cutoff = 4.0
    ProtoSyn.Calculators.update!(ProtoSyn.SIMD_1, vl, pose)
    driver_state = ProtoSyn.Drivers.MonteCarloState{Float64}()
    dihedral_mutator = ProtoSyn.Mutators.DihedralMutator(randn, 0.01, 0.5, an"C|N"r)
    energy_function = ProtoSyn.Common.default_energy_function()
    monte_carlo = ProtoSyn.Drivers.MonteCarlo(energy_function, dihedral_mutator, nothing, 10, ProtoSyn.Drivers.get_linear_quench(1.0, 10))
    cb = ProtoSyn.Common.default_energy_step_frame_callback(10, "teste.pdb")
end); recursive=true)

doctest(ProtoSyn)