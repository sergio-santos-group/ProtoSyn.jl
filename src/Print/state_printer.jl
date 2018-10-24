@doc raw"""
    as_xyz(io:IO, state::Common.State[, title::String = "mol"])

Print the current [`Common.State`](@ref) as a .xyz file to the output `io`.

# Examples
```julia-repl
julia> Print.as_xyz(file_xyz, state, title = "molecule")
```
"""
function as_xyz(io::IO, state::Common.State, title::String="mol")
    write(io, "$(state.size)\n$title\n")
    for (xyz, metadata) in state
        write(io, @sprintf("%-4s %9.4f %9.4f %9.4f\n", metadata.elem, xyz[1]*10, xyz[2]*10, xyz[3]*10))
    end
end


@doc raw"""
    as_xyz(state::Common.State[, title::String = "mol"])::String

Print the current [`Common.State`](@ref) in .xyz format and returns a String.

# Examples
```julia-repl
julia> Print.as_xyz(state, title = "molecule")
2
 molecule
 N      -0.0040    0.2990    0.0000
 H1      0.1200    1.3010    0.0000
...
```
"""
function as_xyz(state::Common.State, title::String = "mol")::String
    iobuffer = IOBuffer()
    as_xyz(iobuffer, state, title)
    return String(take!(iobuffer))
end


@doc raw"""
    as_pdb(io:IO, state::Common.State[, title::String = "mol", step::Int64 = 1])

Print the current [`Common.State`](@ref) as a .pdb file to the output `io`.

# Examples
```julia-repl
julia> Print.as_pdb(file_xyz, state, title = "molecule", step=2)
```
"""
function as_pdb(io::IO, state::Common.State; title::String="mol", step::Int64=1)
    write(io, "TITLE $title\nMODEL $step\n")
    if length(state.metadata.ss) > 0
        n_sheets = length(filter(ss -> ss.ss_type == Common.SS.SHEET, state.metadata.ss))
        for (index, ss) in enumerate(state.metadata.ss)
            write(io, "$(@sprintf("%-6s", string(ss.ss_type)))")
            write(io, "$(@sprintf("%4d", index))")
            write(io, "$(@sprintf("%4s", ss.name))")
            ss.ss_type == Common.SS.SHEET ? write(io, "$(@sprintf("%2d", n_sheets))") : nothing
            write(io, "$(@sprintf("%4s A", Aux.conv_aa_123(ss.i_res_name)))")
            ss.ss_type == Common.SS.SHEET ? write(io, "$(@sprintf("%4d ", ss.i_res_num))") : write(io, "$(@sprintf("%5d ", ss.i_res_num))")
            write(io, "$(@sprintf("%4s A", Aux.conv_aa_123(ss.f_res_name)))")
            ss.ss_type == Common.SS.SHEET ? write(io, "$(@sprintf("%4d ", ss.f_res_num))") : write(io, "$(@sprintf("%5d ", ss.f_res_num))")
            write(io, "$(@sprintf("%2d", ss.conf))")
            write(io, "$(@sprintf("%36d\n", ss.f_res_num - ss.i_res_num))")
        end
    end
    for (index, (xyz, metadata)) in enumerate(state)
        write(io, "$(@sprintf("%-6s", "ATOM"))")
        write(io, "$(@sprintf("%5d  ", index))")
        write(io, "$(@sprintf("%-3s ", metadata.name[1:min(end, 3)]))")
        write(io, "$(@sprintf("%-3s A", metadata.res_name))")
        write(io, "$(@sprintf("%4d    ", metadata.res_num))")
        write(io, "$(@sprintf("%8.3f", xyz[1]*10))")
        write(io, "$(@sprintf("%8.3f", xyz[2]*10))") #Angstrom
        write(io, "$(@sprintf("%8.3f  1.00  0.00", xyz[3]*10))\n")
    end
    write(io, "TER")
    for (index, (xyz, metadata)) in enumerate(state)
        if metadata.connects != nothing
            write(io, "\nCONECT$(@sprintf("%5d", index))")
            for at in metadata.connects
                write(io, "$(@sprintf("%5d", at))")  
            end
        end
    end
    write(io, "\nENDMDL\n")
end


@doc raw"""
    as_pdb(state::Common.State[, title::String = "mol"])::String

Print the current [`Common.State`](@ref) in .pdb format and returns a String.

# Examples
```julia-repl
julia> Print.as_pdb(state, title = "molecule", step = 2)
TITLE molecule
MODEL 2
SHEET    1  BA 3 ASP A   2  ALA A   8  1                                   6
HELIX    2  HA VAL A   13  GLY A   24  1                                  11
ATOM      1  N   GLU A   1      -0.004   0.299   0.000  1.00  0.00
ATOM      2  H1  GLU A   1       0.120   1.301   0.000  1.00  0.00
...
"""
function as_pdb(state::Common.State, title::String = "mol")
    iobuffer = IOBuffer()
    as_pdb(iobuffer, state, title)
    return String(take!(iobuffer))
end