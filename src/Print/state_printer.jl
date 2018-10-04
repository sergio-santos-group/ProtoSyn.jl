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