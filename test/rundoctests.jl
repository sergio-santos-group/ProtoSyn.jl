# * This testing script evaluates all `jldoctest` entries in the ProtoSyn
# * documentation, testing to see if the displayed output matches the expected
# * return value. The necessary objects and variables are set using
# * `DocMeta.setdocmeta!`. In order to avoid errors, run this script after the
# * regular test script `runtests.jl`. 

using Test, Documenter, ProtoSyn
DocMeta.setdocmeta!(ProtoSyn, :DocTestSetup, :(begin
    using ProtoSyn;
    using Random; Random.seed!(1);
    res_lib  = ProtoSyn.Peptides.grammar(Float64, verbose = false);
    frag     = fragment(res_lib, seq"AAA");
    pose     = build(res_lib, seq"SESEAEFKQRLAAIKTRLQAL"); sync!(pose);
    pose_mod = copy(pose);
    ProtoSyn.setdihedral!(pose_mod.state, pose_mod.graph[1][2]["N"], deg2rad(90));
    sync!(pose_mod);
    rrbm = ProtoSyn.Mutators.RotationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, randn, ProtoSyn.center_of_mass, 0.4, an"CA" & an"CB")
    trbm = ProtoSyn.Mutators.TranslationRigidBodyMutator(ProtoSyn.rand_vector_in_sphere, 1.0, nothing)
end); recursive=true)

doctest(ProtoSyn)