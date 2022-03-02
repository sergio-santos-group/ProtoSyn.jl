"""
    unbond!(pose::Pose, residue_1::Residue, residue_2::Residue; [keep_downstream_position::Bool = true])

Unbond the two provided residues by severing the peptidic bond between the C
atom of `residue_1` and the N atom of `residue_2`. In order to do this, perform
the following steps: unset parent/children, unbond neighbours, remove from
[Graph](@ref graph-types), remove from [State](@ref state-types), update the
containers `itemsbyname` field. If `keep_downstream_position` is set to `true`
(is, by default), the downstream [`Residue`](@ref) position is maintained
(by calling [`request_c2i!`](@ref ProtoSyn.request_c2i!) and
[`sync!`](@ref ProtoSyn.sync!) methods).

# Examples
```jldoctest
julia> ProtoSyn.Peptides.unbond!(pose, pose.graph[1][2], pose.graph[1][3]; keep_downstream_position = true)
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function unbond!(pose::Pose, residue_1::Residue, residue_2::Residue; keep_downstream_position::Bool = true)

    if isparent(residue_1, residue_2)
        return ProtoSyn._unbond!(pose, residue_1["C"], residue_2["N"], keep_downstream_position = keep_downstream_position)
    end
    if isparent(residue_2, residue_1)
        return ProtoSyn._unbond!(pose, residue_2["C"], residue_1["N"], keep_downstream_position = keep_downstream_position)
    end
end


"""
# TODO
# SHOULD USE FIND_TAUTOMER?
"""
function assign_default_atom_names!(pose::Pose, grammar::LGrammar = Peptides.grammar)
    
    pose_sequence = [string(residue) for residue in ProtoSyn.sequence(pose)]
    template = ProtoSyn.build(grammar, pose_sequence)

    for segment in eachsegment(pose.graph)
        # Check caps
        # N-terminal
        N = ProtoSyn.identify_atom_by_bonding_pattern(segment[1], ["N", "C", "C", "O"])
        if length(N.bonds) > 2 
            Peptides.cap!(template, rid"1")
        end

        last_res_id = segment[end].id
        # Identify C
        Cs = ProtoSyn.identify_atom_by_bonding_pattern(segment[end], ["C", "C", "N", "C"])
        if isa(Cs, Atom)
            C = Cs
        else
            C = [c for c in Cs if (length(c.bonds) <= 3)]
            if length(C) > 1
                C = [c for c in C if c.container["CA"] in c.bonds]
            end
            C = C[1]
        end
        if length(C.bonds) > 2 
            Peptides.cap!(template, SerialSelection{Residue}(last_res_id, :id))
        end
    end
    
    t_atoms = ProtoSyn.travel_graph(template.graph[1, 1, 1], search_algorithm = ProtoSyn.BFS)
    p_atoms = ProtoSyn.travel_graph(pose.graph[1, 1, 1], search_algorithm = ProtoSyn.BFS)
    
    if all([a.symbol for a in t_atoms] .=== [a.symbol for a in p_atoms])
        for (t_atom, p_atom) in zip(t_atoms, p_atoms)
            ProtoSyn.rename!(p_atom, t_atom.name)
        end
    end
    
    # In case the direct approach wasn't successful
    if ProtoSyn.verbose.mode
        @warn "Possible tautomers identified, assigning default names residue by residue ..."
    end
    pose_temp_res = zip(eachresidue(pose.graph), eachresidue(template.graph))
    for (pose_residue, template_residue) in pose_temp_res
        # println("\n$pose_residue - $template_residue")

        pose_res_name = string(ProtoSyn.three_2_one[pose_residue.name.content])
        ind_pose_res  = copy(pose_residue) # Removes inter-residue connections
        pose_N        = ProtoSyn.identify_atom_by_bonding_pattern(ind_pose_res, ["N", "C", "C", "O"])
        pose_atoms    = ProtoSyn.travel_graph(pose_N, search_algorithm = ProtoSyn.BFS)
        pose_graph    = [a.symbol for a in pose_atoms]
        # println("Residue $pose_residue (Graph: $pose_graph)")

        if isa(grammar.variables[pose_res_name], Tautomer)
            found_match = false
            for tautomer in grammar.variables[pose_res_name].list
                temp_N     = ProtoSyn.identify_atom_by_bonding_pattern(tautomer.graph[1], ["N", "C", "C", "O"])
                temp_atoms = ProtoSyn.travel_graph(temp_N, search_algorithm = ProtoSyn.BFS)
                temp_graph = [a.symbol for a in temp_atoms]

                if all(pose_graph .=== temp_graph)
                    # println("Match found.")
                    for (t_atom, p_atom) in zip(temp_atoms, pose_atoms)
                        pri = pose_residue.items
                        actual_pose_a = findfirst((a)->a.index == p_atom.index, pri)
                        ProtoSyn.rename!(pri[actual_pose_a], t_atom.name)
                    end
                    found_match = true
                end
            end
            !found_match && @warn "No available tautomer matched the residue $ind_pose_res !"
        else
            temp_res   = copy(template_residue)
            temp_N     = ProtoSyn.identify_atom_by_bonding_pattern(temp_res, ["N", "C", "C", "O"])
            temp_atoms = ProtoSyn.travel_graph(temp_N, search_algorithm = ProtoSyn.BFS)
            temp_graph = [a.symbol for a in temp_atoms]

            if all(pose_graph .=== temp_graph)
                # println("Match found.")
                for (t_atom, p_atom) in zip(temp_atoms, pose_atoms)
                    pri = pose_residue.items
                    actual_pose_a = findfirst((a)->a.index == p_atom.index, pri)
                    # println("renaming atom $(pri[actual_pose_a]) to $(t_atom.name)")
                    ProtoSyn.rename!(pri[actual_pose_a], t_atom.name)
                end
            else
                @warn "Available template for $pose_residue ($(template_residue)) doesn't match the graph."
            end
        end
    end

    return pose
end