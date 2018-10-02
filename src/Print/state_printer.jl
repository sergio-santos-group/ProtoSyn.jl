#TO DO:
# 1) DOCUMENT THE FUNCTION
function as_xyz(state::Common.State;
    ostream::IO = stdout, title::String = "mol")

    atom_count = size(state.xyz, 1)
    write(ostream, "$atom_count\n $title\n")
    for atom_index in 1:atom_count
        write(ostream, " $(@sprintf("%-4s", state.atnames[atom_index]))")
        write(ostream, " $(@sprintf("%9.4f", state.xyz[atom_index, :][1]*10))")
        write(ostream, " $(@sprintf("%9.4f", state.xyz[atom_index, :][2]*10))") #Angstrom
        write(ostream, " $(@sprintf("%9.4f", state.xyz[atom_index, :][3]*10))\n")
    end
end