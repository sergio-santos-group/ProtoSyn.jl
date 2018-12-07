@doc raw"""
    apply_ss!(state::State, dihedrals::Vector{Dihedral}, ss::String)

Apply predefined angles to all `PHI` and `PSI` dihedrals defined in `dihedrals`, based on the provided `ss`, changing the State.xyz
to apply the secondary structure. The applied angles (in degrees) are the following:

Beta sheet ("E"):
PHI = -139.0 | PSI = 135.0
Alpha helix ("H"):
PHI = -57.0  | PSI = -47.0

# Examples
```julia-repl
julia> Common.apply_ss!(state, dihedrals, "CCCHHHHCCEEEECCC")
```
See also: [`compile_ss`](@ref)
"""
function apply_ss!(state::State, metadata::Metadata, ss::String)

    # Save secondary structure as metadata
    metadata.ss     = compile_ss(metadata.dihedrals, ss)
    metadata.blocks = compile_blocks(metadata.residues, ss)

    index::Int64 = 1
    for dihedral in metadata.dihedrals
        #Verify input so that only the phi and psi dihedrals are iterated over, and COILS remain "free"
        if !(dihedral.dtype >= DIHEDRAL.omega || dihedral.residue.ss == Common.SS.COIL)
            rotate_dihedral_to!(state.xyz, dihedral, Common.ss2bbd[dihedral.residue.ss][dihedral.dtype])
        end
    end
    println("(  PRE) ▲ Applied secondary structure angles to $(length(metadata.blocks)) blocks")
end


@doc raw"""
    apply_dihedrals_from_file!(state::State, bb_dihedrals::Vector{Dihedral}, file_i::String)

Read the input file `file_i`, extract the dihedrals angles and apply to the backbone dihedrals `bb_dihedrals` of `state`.

# Examples
```julia-repl
julia> Common.apply_initial_conf(state, bb_dihedrals, "native_conf.pdb")
```
"""
function apply_dihedrals_from_file!(state::State, bb_dihedrals::Vector{Dihedral}, file_i::String)

    #Extract backbone coordinates
    backbone = Vector{Vector{Float64}}()
    open(file_i, "r") do f
        for line in eachline(f)
            if startswith(line, "ATOM") && string(strip(line[13:16])) in ["N", "CA", "C"]
                push!(backbone, map(x -> 0.1*parse(Float64, x), [line[31:38], line[39:46], line[47:54]]))
            end
        end
    end

    #Calculate native dihedrals and rotate on state
    count::Int64 = 0
    cd = 1
    for i in 1:(length(backbone) - 3)
        if i == (2 + 3 * count)
            count += 1
            continue
        end
        current_angle = Aux.calc_dih_angle(state.xyz[bb_dihedrals[cd].a1, :], state.xyz[bb_dihedrals[cd].a2, :], state.xyz[bb_dihedrals[cd].a3, :], state.xyz[bb_dihedrals[cd].a4, :])
        displacement = Aux.calc_dih_angle(backbone[i], backbone[i + 1], backbone[i + 2], backbone[i + 3]) - current_angle
        rotate_dihedral!(state.xyz, bb_dihedrals[cd], displacement)
        cd += 1
    end

end


@doc raw"""
    stretch_conformation!(state::State, dihedrals::Vector{Dihedral})

Iterate over the `dihedrals` list and rotate all PHI and PSI to 180º.

# Examples
```julia-repl
julia> Common.stretch_conformation!(state, dihedrals)
```
"""
function stretch_conformation!(state::State, dihedrals::Vector{Dihedral})

    for dihedral in dihedrals
        dihedral.dtype > omega ? continue : nothing
        rotate_dihedral_to!(state.xyz, dihedral, 3.1416)
    end
end


@doc raw"""
    fix_proline!(state::State, dihedrals::Vector{Dihedral})

( TEMPORARY ) Find all Prolines in the `dihedrals` list and solve the two following issues:
1. Add the first two atoms to the movables list.
2. Rotate the dihedral to be at 0º instead of 180º.

# Examples
```julia-repl
julia> Common.fix_proline(state, dihedrals)
```
"""
function fix_proline!(state::State, dihedrals::Vector{Dihedral})

    for dihedral in dihedrals
        if Aux.conv123(dihedral.residue.name) == "PRO" && dihedral.dtype == Common.DIHEDRAL.omega
            rotate_dihedral!(state.xyz, dihedral.a3, dihedral.a2, deg2rad(180), Common.DIHEDRAL.omega, Vector{Int64}(), dihedral.residue)
        end
    end
end