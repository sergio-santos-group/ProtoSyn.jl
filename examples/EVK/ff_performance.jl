using ProtoSyn

#Load state from PDB file
state = Common.load_from_pdb("data/EVK.pdb")
state.energy = Forcefield.Energy()

#Load topologies from JSON files
topol = Forcefield.load_from_json("data/EVK.json")

function do_work(n::Int64)
    for i = 1:n
        fill!(state.forces, 0.0)
        Forcefield.evaluate!(topol.bonds, state, do_forces=true)
        Forcefield.evaluate!(topol.angles, state, do_forces=true)
        Forcefield.evaluate!(topol.dihedralsCos, state, do_forces=true)
        Forcefield.evaluate!(topol.atoms, state, do_forces=true, cut_off=Inf)
    end
end

# force JIT compilation
do_work(2)

@time do_work(10000)