using ProtoSyn.Peptides

"""
    diagnose(pose::Pose; [return_issues::Bool = false], [atom_order_search_algorithm::F = Peptides.IUPAC], [res_lib::LGrammar = Peptides.grammar]) where {F <: ProtoSyn.SearchAlgorithm}

Measure several agreement criteria on the given [`Pose`](@ref) `pose`:

 1. All diagnostics from [`ProtoSyn.diagnose`](@ref), except using the given `atom_order_search_algorithm` (`Peptides.IUPAC`, by default)
 2. Check if [`Atom`](@ref) and [`Residue`](@ref) names are in accordance to the templates in the provided [`LGrammar`] `res_lib`
 3. Check for missing hydrogens (in comparison with the templates in the provided [`LGrammar`] `res_lib`)
 4. Check for missing caps.
   
If `return_issues` is set to `true` (`false`, by default) doesn't print results
to `stdout`, returns them as a `Vector{String}` instead.

# Examples
```
julia> ProtoSyn.Peptides.diagnose(pose)
 ⬤   Diagnosing pose 2a3d ...
 |
 ├──  • Residue-level graph OK
 |
 ├──  • Pose indexation OK
 |
 ├──  • Atom-level graph OK
 |
 ├──  • Pose synchronization status (1 issue identified)
 |    └── Pose is requesting internal to cartesian coordinates synchronization. Consider using the sync! function.
 |
 ├──  • Pose charges (1 issue identified)
 |    └── Sum of charges in the pose is not 0.0.
 |
 ├──  • (Peptides only) Residue and atom naming OK
 |
 └──  • (Peptides only) Missing atoms (1 issue identified)
      └── Pose seems to be missing hydrogens. Suggested fix: Consider using ProtoSyn.add_hydrogens!
```
"""
function diagnose(pose::Pose; return_issues::Bool = false, atom_order_search_algorithm::F = Peptides.IUPAC, res_lib::LGrammar = Peptides.grammar) where {F <: ProtoSyn.SearchAlgorithm}
    
    default_issues = ProtoSyn.diagnose(pose, return_issues = true,
        atom_order_search_algorithm = atom_order_search_algorithm)

    init_level_code = ProtoSyn.LevelCode()
    for (title, issues) in default_issues
        !return_issues && ProtoSyn.print_diagnose_results(title, issues, init_level_code, 3)
    end # for

    # 1. Scan for unknown residue and atom names
    graph_naming_issues = Vector{String}()
    ukws_res    = Vector{Int}()
    ukws_res_n  = Vector{String}()
    ukw_atoms   = Vector{Int}()
    # Get terminal atoms for cap checking
    N = UpstreamTerminalSelection{Residue}()(pose, gather = true)[1]["N"]
    C = DownstreamTerminalSelection{Residue}()(pose, gather = true)[1]["C"]
    for res in eachresidue(pose.graph)
        ukw_res = !(res.name in keys(ProtoSyn.three_2_one))
        if ukw_res
            push!(ukws_res, res.id)
            push!(ukws_res_n, res.name.content)
            continue
        end # if
        code = ProtoSyn.three_2_one[res.name]
        ukw_res = !(code in keys(Peptides.available_aminoacids))
        if ukw_res
            push!(ukws_res, res.id)
            push!(ukws_res_n, res.name.content)
            continue
        end # if
        
        template = res_lib.variables[string(code)]
        if isa(template, Tautomer)
            tautomer = ProtoSyn.find_tautomer(template, res)
            template = tautomer === nothing ? template.list[1] : tautomer
        end # if
        atom_names = [atom.name for atom in template.graph[1].items]
        for atom in eachatom(res)
            ukw_atm = !(atom.name in atom_names)
            # Check if the ukw atoms are caps. A cap is defined as any H or O
            # atom connected to the N-terminal or C-terminal, respectivelly.
            if ukw_atm
                is_n_cap = (atom.symbol === "H") & (N in atom.bonds)
                is_c_cap = (atom.symbol === "O") & (C in atom.bonds)
                if !is_n_cap & !is_c_cap
                    push!(ukw_atoms, atom.id)
                end # if
            end # if
        end # for
    end # for
    if length(ukws_res) > 0
        s = length(ukws_res) > 1 ? "s" : ""
        a_info = ""
        if length(ukws_res) < 10
            a_info = "\nCheck residue$s $(Base.join(["$i-$n" for (i, n) in zip(ukws_res, ukws_res_n)], ", ", " and "))."
        end
        push!(graph_naming_issues, "Unknown residue$s name$s (Not found in Peptides.available_aminoacids or ProtoSyn.three_2_one).\nThis can be due to the presence of ligands or NCAAs.$a_info")
    end # if
    ai = "" # additional information
    if length(ukw_atoms) < 10
        ai = "\nAdditional information: Check atoms $(Base.join(ukw_atoms, ", ", " and ")) (by :id)."
    end # if
    if length(ukw_atoms) > 0
        s = length(ukw_atoms) > 1 ? "s" : ""
        push!(graph_naming_issues, "Unknown atom$s name$s (Not found in the given LGrammar templates).\nTotal number of atoms with unknown names: $(length(ukw_atoms))/$(ProtoSyn.count_atoms(pose.graph)).\nSuggested fix: Consider using ProtoSyn.Peptides.assign_default_atom_names!$ai")
    end # if

    !return_issues && ProtoSyn.print_diagnose_results("(Peptides only) Residue and atom naming", graph_naming_issues, init_level_code, 3)

    # 2. Has hydrogens
    missing_atoms_residues = Vector{String}()
    h = [a for a in eachatom(pose.graph) if a.symbol === "H"]
    if length(h) <= ProtoSyn.count_residues(pose.graph)
        push!(missing_atoms_residues, "Pose seems to be missing hydrogens. Suggested fix: Consider using ProtoSyn.add_hydrogens!")
    end # if
    sc = (SidechainSelection() & !as"H")(pose, gather = true)
    if length(sc) <= ProtoSyn.count_residues(pose.graph)
        push!(missing_atoms_residues, "Pose seems to be missing 1 or more sidechains. Suggested fix: Consider using ProtoSyn.Peptides.add_sidechains!")
    end # if

    !return_issues && ProtoSyn.print_diagnose_results("(Peptides only) Missing atoms", missing_atoms_residues, init_level_code, 3)

    # 3. Check caps
    missing_caps = Vector{String}()
    for segment in eachsegment(pose.graph)
        N_terminal = UpstreamTerminalSelection{Residue}()(segment, gather = true)[1]
        N = ProtoSyn.identify_atom_by_bonding_pattern(N_terminal, ["N", "C", "C", "O"])
        if !isa(N, Vector{Atom}) && length(N.bonds) <= 2
            push!(missing_caps, "Pose seems to be missing the N-terminal NH₃ cap. Suggested fix: Consider using ProtoSyn.Peptides.cap!")
        end

        C_terminal = DownstreamTerminalSelection{Residue}()(segment, gather = true)[1]
        C = identify_c_terminal(segment, supress_warn = true)
        if !isa(C, Vector{Atom}) && length(C.bonds) <= 2
            push!(missing_caps, "Pose seems to be missing the C-terminal CO₂ cap. Suggested fix: Consider using ProtoSyn.Peptides.cap!")
        end
    end

    !return_issues && ProtoSyn.print_diagnose_results("(Peptides only) Missing caps", missing_caps, init_level_code, 4)
end # function