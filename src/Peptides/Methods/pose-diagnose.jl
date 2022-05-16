using ProtoSyn.Peptides

"""
# TODO
"""
function diagnose(pose::Pose; return_issues::Bool = false, atom_order_search_algorithm::F = Peptides.IUPAC, res_lib::LGrammar = Peptides.grammar) where {F <: ProtoSyn.SearchAlgorithm}
    
    default_issues = ProtoSyn.diagnose(pose, return_issues = true,
        atom_order_search_algorithm = atom_order_search_algorithm)

    init_level_code = ProtoSyn.LevelCode()
    for (title, issues) in default_issues
        !return_issues && ProtoSyn.print_diagnose_results(title, issues, init_level_code, 3)
    end

    # 1. Scan for unknown residue and atom names
    graph_naming_issues = Vector{String}()
    ukws_res  = Vector{Int}()
    ukw_atoms = Vector{Int}()
    # Get terminal atoms for cap checking
    N = NTerminalSelection()(pose, gather = true)[1]["N"]
    C = CTerminalSelection()(pose, gather = true)[1]["C"]
    for res in eachresidue(pose.graph)
        ukw_res = !(res.name in keys(ProtoSyn.three_2_one))
        if ukw_res
            push!(ukws_res, res.id)
            continue
        end
        code = ProtoSyn.three_2_one[res.name]
        ukw_res = !(code in keys(Peptides.available_aminoacids))
        if ukw_res
            push!(ukws_res, res.id)
            continue
        end
        
        template = res_lib.variables[string(code)]
        if isa(template, Tautomer)
            tautomer = ProtoSyn.find_tautomer(template, res)
            template = tautomer === nothing ? template.list[1] : tautomer
        end
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
                end
            end
        end
    end
    if length(ukws_res) > 0
        s = length(ukws_res) > 1 ? "s" : ""
        push!(graph_naming_issues, "Unknown residue$s name$s (Not found in Peptides.available_aminoacids or ProtoSyn.three_2_one).\nCheck residue$s $(Base.join(ukws_res, ", ", " and ")).")
    end
    ai = "" # additional information
    if length(ukw_atoms) < 10
        ai = "\nAdditional information: Check atoms $(Base.join(ukw_atoms, ", ", " and ")) (by :id)."
    end
    if length(ukw_atoms) > 0
        s = length(ukw_atoms) > 1 ? "s" : ""
        push!(graph_naming_issues, "Unknown atom$s name$s (Not found in the given LGrammar templates).\nTotal number of atoms with unknown names: $(length(ukw_atoms))/$(ProtoSyn.count_atoms(pose.graph)).\nSuggested fix: Consider using ProtoSyn.Peptides.assign_default_atom_names!$ai")
    end

    !return_issues && ProtoSyn.print_diagnose_results("(Peptides only) Residue and atom naming", graph_naming_issues, init_level_code, 3)

    # 2. Has hydrogens
    missing_atoms_residues = Vector{String}()
    h = [a for a in eachatom(pose.graph) if a.symbol === "H"]
    if length(h) <= ProtoSyn.count_residues(pose.graph)
        push!(missing_atoms_residues, "Pose seems to be missing hydrogens. Suggested fix: Consider using ProtoSyn.add_hydrogens!")
    end
    sc = (SidechainSelection() & !as"H")(pose, gather = true)
    if length(sc) <= ProtoSyn.count_residues(pose.graph)
        push!(missing_atoms_residues, "Pose seems to be missing 1 or more sidechains. Suggested fix: Consider using ProtoSyn.Peptides.add_sidechains!")
    end

    !return_issues && ProtoSyn.print_diagnose_results("(Peptides only) Missing atoms", missing_atoms_residues, init_level_code, 4)

end