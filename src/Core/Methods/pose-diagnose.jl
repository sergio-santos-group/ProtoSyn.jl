"""
    print_diagnose_results(title::String, issues::Vector{String}, init_level_code::LevelCode, level::Int)

Prints a list of `issues` (under the scop of `title` grouping) in a stylized way
(uses `init_level_code` and adds inner `LevelCode` instances at `level`). 

!!! ukw "Note:"
    This function is intended to be employed from [`diagnose`](@ref), and not as a standalone.

# Examples
```
julia> ProtoSyn.print_diagnose_results("Charges", ["Missing charges", "Charges sum is not 0.0"], LevelCode(), 4)
 |
 └──  • Charges (2 issues identified)
      ├── Missing charges
      └── Charges sum is not 0.0
```
"""
function print_diagnose_results(title::String, issues::Vector{String}, init_level_code::LevelCode, level::Int)
    level_code        = vcat(init_level_code, level)
    lead              = ProtoSyn.get_lead(level_code)
    N                 = length(issues)
    s                 = N > 1 ? "s" : ""
    result = N == 0 ? "OK" : "($N issue$s identified)"
    issue_color       = N == 0 ? :green : :yellow
    println(" |")
    print(lead*" • $title ")
    printstyled(result*"\n", color = issue_color)
    if N > 0
        for issue in issues[1:(end-1)]
            inner_level_code = vcat(level_code, 3)
            lead             = ProtoSyn.get_lead(inner_level_code)
            inner_lead       = ProtoSyn.get_inner_lead(inner_level_code)
            issue_parts      = split(issue, "\n")
            println(lead, issue_parts[1])
            if length(issue_parts) > 1
                for issue_part in issue_parts[2:end]
                    println(inner_lead*issue_part)
                end
            end
        end

        inner_level_code = vcat(level_code, 4)
        lead             = ProtoSyn.get_lead(inner_level_code)
        inner_lead       = ProtoSyn.get_inner_lead(inner_level_code)
        issue_parts      = split(issues[end], "\n")
        println(lead, issue_parts[1])
        if length(issue_parts) > 1
            for issue_part in issue_parts[2:end]
                println(inner_lead*issue_part)
            end
        end
    end
end


"""
    diagnose(pose::Pose; [return_issues::Bool = false], [atom_order_search_algorithm::F = ProtoSyn.BFS]) where {F <: SearchAlgorithm}

Measure several agreement criteria on the given [`Pose`](@ref) `pose`:

 1. Checks residue-level graph for any [`Residue`](@ref) instance without parent
 2. Checks atom-level graph for any [`Atom`](@ref) instance without parent
 3. Checks if atom-level graph travels all [`Atom`](@ref) instances in the given [`Pose`](@ref) `pose`
 4. Checks if atom-level graph and list of [`Atom`](@ref) instances in the given [`Pose`](@ref) `pose` have the same order
 5. Checks if any internal to cartesian coordinate conversion (or vice-versa) is pending
 6. Check if the [`Pose`](@ref) `pose` indexation matches the order of atoms (both in the :id and :index fields)

Any [Graph](@ref graph-types) travel is done using the [`ProtoSyn.travel_graph`](@ref)
method, employing the given `atom_order_search_algorithm` (`ProtoSyn.BFS`, by
default). If `return_issues` is set to `true` (`false`, by default) doesn't
print results to `stdout`, returns them as a `Vector{String}` instead.

# Examples
```
julia> ProtoSyn.diagnose(pose)
 ⬤   Diagnosing pose 4J88 ...
 |
 ├──  • Residue-level graph OK
 |
 ├──  • Atom-level graph (1 issue identified)
 |    └── Travelling from the first atom on the pose list, not all atoms were visited (From graph: 466 | From list: 1805).
 |        Check for breaks in the parenthood relationships. Suggested fix: consider using the infer_parenthood! function.
 |
 ├──  • Pose synchronization status OK
 |
 └──  • Pose indexation (1 issue identified)
      └── Pose indexation doesn't match the current order of atoms in each AbstractContainer (on the :id fields).
          Check atoms 63 and 64. Suggested fix: Consider using the reindex function.
```
"""
function diagnose(pose::Pose; return_issues::Bool = false, atom_order_search_algorithm::F = ProtoSyn.BFS) where {F <: SearchAlgorithm}

    println(" ⬤   Diagnosing pose $(pose.graph.name) ...")
    init_level_code = ProtoSyn.LevelCode()

    # 1. Check residue-level parenthoods for missing parents
    residue_level_graph_issues = Vector{String}()
    for residue in eachresidue(pose.graph)
        if !hasparent(residue)
            push!(residue_level_graph_issues, "No parent on residue $residue") 
        end
    end

    !return_issues && print_diagnose_results("Residue-level graph", residue_level_graph_issues, init_level_code, 3)

    # 2. Check atom-level parenthoods for missing parents
    atom_level_graph_issues = Vector{String}()
    for atom in eachatom(pose.graph)
        if !hasparent(atom)
            push!(atom_level_graph_issues, "No parent on atom $atom") 
        end
    end

    # 3. Check if atom-level graph matches atom list (number and order)
    for segment in eachsegment(pose.graph)
        atoms_from_graph = ProtoSyn.travel_graph(segment[1][1], search_algorithm = atom_order_search_algorithm)
        atoms_from_list  = collect(eachatom(segment))
        NAG              = length(atoms_from_graph)
        NAL              = length(atoms_from_list)
        if !(NAG === NAL)
            push!(atom_level_graph_issues, "Using the $atom_order_search_algorithm graph search algorithm and travelling from the first atom on $segment, not all atoms were visited (From graph: $NAG | From list: $NAL).\nCheck for breaks in the parenthood relationships. Suggested fix: consider using the infer_parenthood! function.")
        else
            N = findall((x) -> x === false, atoms_from_list .=== atoms_from_graph)
            if length(N) > 0
                push!(atom_level_graph_issues, "Using the $atom_order_search_algorithm graph search algorithm and travelling from the first atom on $segment:\nThe atom order doesn't match the established list of atoms in each AbstractContainer.\nTotal number of atoms in different order: $(length(N))/$(length(atoms_from_list)).\nSuggested fix: Consider using the sort_atoms_by_graph! function (with the desired search algorithm convention) or diagnose using a different graph search algorithm.")
            end
        end
    end

    !return_issues && print_diagnose_results("Atom-level graph", atom_level_graph_issues, init_level_code, 3)

    # 4. Pose synchronization status
    pose_synchronization_issues = Vector{String}()
    if pose.state.i2c
        push!(pose_synchronization_issues, "Pose is requesting internal to cartesian coordinates synchronization. Consider using the sync! function.")
    end
    if pose.state.c2i
        push!(pose_synchronization_issues, "Pose is requesting cartesian to internal coordinates synchronization. Consider using the sync! function.")
    end

    !return_issues && print_diagnose_results("Pose synchronization status", pose_synchronization_issues, init_level_code, 3)

    # 5. Indexation
    indexation_issues = Vector{String}()

    atoms_from_list = collect(eachatom(pose.graph))
    ids = [a.id for a in atoms_from_list]
    N   = findall((x) -> x === false, ids .=== sort(ids))
    if length(N) > 0
        push!(indexation_issues, "Pose indexation doesn't match the current order of atoms in each AbstractContainer (on the :id fields).\nCheck atoms $(Base.join(N, ", ", " and ")). Suggested fix: Consider using the reindex function.")
    end

    indexes = [a.index for a in atoms_from_list]
    N       = findall((x) -> x === false, indexes .=== sort(indexes))
    if length(N) > 0
        push!(indexation_issues, "Pose indexation doesn't match the current order of atoms in each AbstractContainer (on the :index fields).\nCheck atoms $(Base.join(N, ", ", " and ")). Suggested fix: Consider using the reindex function.")
    end

    !return_issues && print_diagnose_results("Pose indexation", indexation_issues, init_level_code, 3)

    # 6. Charges
    charge_issues = Vector{String}()
    c = [atom.δ for atom in pose.state.items]
    if all((x) -> x === 0.0, c)
        push!(charge_issues, "Pose doesn't seem to have charges assigned.\nSuggested fix: Consider using ProtoSyn.Calculators.Electrostatics.assign_default_charges!")
    end
    if !isapprox(sum(c), 0.0, atol = 1e-10)
        push!(charge_issues, @sprintf "Sum of charges in the pose is not 0.0 (Δδ = %.3f)." sum(c))
    end

    !return_issues && print_diagnose_results("Pose charges", charge_issues, init_level_code, 4)


    if return_issues
        return Dict{String, Vector{String}}(
            "Residue-level graph"         => residue_level_graph_issues,
            "Atom-level graph"            => atom_level_graph_issues,
            "Pose synchronization status" => pose_synchronization_issues,
            "Pose indexation"             => indexation_issues,
            "Pose charges"                => charge_issues,
            )
    else
        return nothing
    end
end


"""
    set_parenthood_as_forces!(pose::Pose, [selection::Opt{AbstractSelection} = nothing])

Set a vector between any given [`Atom`](@ref) instance in the provided
[`Pose`](@ref) `pose` and its `.parent` as that [`Atom`](@ref) instance's force
(in the `pose.state`). If an `AbstractSelection` `selection` is given, only loop
over the selected atoms. Useful for debuging purposes.

# See also
[`write_forces`](@ref)

# Examples
```
julia> ProtoSyn.set_parenthood_as_forces!(pose)

julia> pose.state.f
3×7 Matrix{Float64}:
  0.512436     -1.4          -0.7          -1.43         0.7           1.4           0.7
  0.512436      0.0           1.21244      -1.46422e-7   1.21244      -4.88498e-15  -1.21244
 -6.27553e-17   1.21234e-16  -1.31353e-16   8.09416e-7  -2.52586e-16  -1.21234e-16   1.31353e-16
```
"""
function set_parenthood_as_forces!(pose::Pose, selection::Opt{AbstractSelection} = nothing)

    if selection === nothing
        sele = TrueSelection{Atom}()
    else
        sele = ProtoSyn.promote(selection, Atom)
    end

    for atom in sele(pose, gather = true)
        s = collect(pose.state[atom].t)
        pose.state.f[:, atom.index] .= collect(pose.state[atom.parent].t) .- s
    end
end


"""
    write_forces(pose::Pose, filename::String, α::T = 1.0) where {T <: AbstractFloat}

Write the [`Pose`](@ref) `pose` forces to `filename` in a specific format to be
read by the companion Python script "cgo_arrow.py" (see
[https://pymolwiki.org/index.php/Cgo_arrow](https://pymolwiki.org/index.php/Cgo_arrow)).
`α` sets a multiplying factor to make the resulting force vectors longer/shorter
(for visualization purposes only).

# See also
[`set_parenthood_as_forces!`](@ref)

# Examples
```
julia> ProtoSyn.write_forces(pose, "forces.dat")
```
"""
function write_forces(pose::Pose, filename::String, α::T = 1.0) where {T <: AbstractFloat}
    open(filename, "w") do file_out
        for (i, atom) in enumerate(eachatom(pose.graph))
            x  = pose.state[atom].t[1]
            y  = pose.state[atom].t[2]
            z  = pose.state[atom].t[3]
            fx = x + (pose.state.f[1, i] * α)
            fy = y + (pose.state.f[2, i] * α)
            fz = z + (pose.state.f[3, i] * α)
            
            s  = @sprintf("%5d %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n", atom.id, x, y, z, fx, fy, fz)
            Base.write(file_out, s)
        end
    end
end