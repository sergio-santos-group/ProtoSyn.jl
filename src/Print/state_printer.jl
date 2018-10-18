@doc raw"""
    as_xyz(state::Common.State[, ostream::IO = stdout, title::String = "mol"])

Print the current [`Common.State`](@ref) as a .xyz file.

# Examples
```julia-repl
julia> Drivers.MonteCarlo.load_parameters(state, title = "molecule")
2
 molecule
 N      -0.0040    0.2990    0.0000
 H1      0.1200    1.3010    0.0000
```
"""
function as_xyz(io::IO, state::Common.State, title::String="mol")
    xyz = 10*state.xyz
    write(io, "$(state.size)\n$title\n")
    for i in 1:state.size
        write(io, "$(@sprintf("%-4s %9.4f %9.4f %9.4f\n", state.atnames[i][1], xyz[i,1], xyz[i,2], xyz[i,3]))")
        # write(io, " $(@sprintf("%9.4f", xyz[i,1]))")
        # write(io, " $(@sprintf("%9.4f", xyz[i,2]))") #Angstrom
        # write(io, " $(@sprintf("%9.4f", xyz[i,3]))\n")
    end
end

function as_xyz(state::Common.State, title::String = "mol")
    iobuffer = IOBuffer()
    as_xyz(iobuffer, state, title)
    return String(take!(iobuffer))
end


#TODO: Function is deprecated. Either return it to LIVE (and document it) or DELETE.
function as_pdb(io::IO, state::Common.State, title::String="mol")

    write(io, "TITLE $title\nMODEL\n")
    for atom_index in 1:state.size
        write(io, "ATOM")
        write(io, "$(@sprintf("%7d ", atom_index))")
        write(io, "$(@sprintf("%-5s", state.atnames[atom_index]))")
        write(io, "$(@sprintf("%-3s A", "ALA"))")
        write(io, "$(@sprintf("  %-3d   ", 1))")
        # write(io, "$(@sprintf("%-3s A", state.residues[atom_index][2]))")
        # write(io, "$(@sprintf("  %-3d   ", state.residues[atom_index][1]))")
        write(io, "$(@sprintf("%8.3f", state.xyz[atom_index, :][1]*10))")
        write(io, "$(@sprintf("%8.3f", state.xyz[atom_index, :][2]*10))") #Angstrom
        write(io, "$(@sprintf("%8.3f  1.00  0.00", state.xyz[atom_index, :][3]*10))\n")
    end
    # write(ostream, "TER")
    # for list in state.conects
    #     write(ostream, "\nCONECT")
    #     for at in list
    #         write(ostream, "$(@sprintf("%5d", at))")  
    #     end
    # end
    write(io, "ENDMDL\n")
end

function as_pdb(state::Common.State, title::String = "mol")
    iobuffer = IOBuffer()
    as_pdb(iobuffer, state, title)
    return String(take!(iobuffer))
end