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
    assign_default_atom_names!(pose::Pose, [selection::Opt{AbstractSelection} = nothing], [grammar::LGrammar = Peptides.grammar]; [force_rename::Bool = false])

Assigns the default [`Atom`](@ref) names to the given [`Pose`](@ref) `pose`. If
an `AbstractSelection` `selection` is provided, only rename the selected
[`Residue`](@ref) instances (any given `selection` will be promoted to be of
[`Residue`](@ref) level using [`ProtoSyn.promote`](@ref)). The [`Atom`](@ref)
default names are retrieved from the given [`LGrammar`](@ref) `grammar`
(`Peptides.grammar`, by default). Both the given [`Pose`](@ref) and the built
template from the `grammar`'s [Graph](@ref graph-types) are travelled to generate a 1-to-1
correspondence between [`Atom`](@ref) instances (the [`Pose`](@ref)
[`Atom`](@ref) names are then renamed to match the template, uses the
[`ProtoSyn.travel_graph`](@ref) method with the `ProtoSyn.Peptides.IUPAC` search
algorithm). This approach may sometimes fail, for example, when tautomers are
present, in wich case this method attempts to find tautomer candidates in a
[`Residue`](@ref) by [`Residue`](@ref) case. This method also attempts to verify
if terminal [`Residue`](@ref) instances are capped, in which case the correct
naming attribution is automatically taken into consideration.

!!! ukw "Note:"
    Some methods (such as, for example, the [`assign_default_charges!`](@ref ProtoSyn.Peptides.Calculators.Electrostatics.assign_default_charges!) method) expect default atom names. Consider using [`ProtoSyn.Peptides.diagnose`](@ref) to check whether the current [`Atom`](@ref) names are known.

# See also
[`rename!`](@ref ProtoSyn.rename!)

# Examples
```
julia> ProtoSyn.Peptides.assign_default_atom_names!(pose)
Pose{Topology}(Topology{/2a3d:42429}, State{Float64}:
 Size: 1140
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function assign_default_atom_names!(pose::Pose, selection::Opt{AbstractSelection} = nothing, grammar::LGrammar = Peptides.grammar; force_rename::Bool = false)

    if selection === nothing
        sele = TrueSelection{Residue}()
    else
        sele = ProtoSyn.promote(selection, Residue)
    end

    for segment in eachsegment(pose.graph)

        residues = sele(segment, gather = true)
        
        if length(residues) === 0
            @info "Skipping $segment : No selected residues in this segment!"
            continue
        end

        pose_sequence = [string(ProtoSyn.three_2_one[residue.name.content]) for residue in residues]
        template = ProtoSyn.build(grammar, pose_sequence)

        # Check cysteines
        cysteines = (sele & rn"CYS" & as"S")(segment, gather = true)
        if length(cysteines) > 0
            template_cysteines = (sele & rn"CYS" & as"S")(template, gather = true)
            for (S_pose, S_template) in zip(cysteines, template_cysteines)
                if !("H" in [a.symbol for a in S_pose.bonds])
                    H = [a for a in S_template.bonds if a.symbol === "H"][1]
                    ProtoSyn.pop_atom!(template, H)
                end
            end
        end

        # Check caps
        # N-terminal
        N_terminal_pose     = UpstreamTerminalSelection{Residue}()(segment, gather = true)[1]
        N_terminal_template = UpstreamTerminalSelection{Residue}()(template, gather = true)[1]
        if !(N_terminal_pose in residues)
            @info "N-terminal not identified as being selected: skipping cap check on this residue."
        else
            N = ProtoSyn.identify_atom_by_bonding_pattern(N_terminal_pose, ["N", "C", "C", "O"])
            if isa(N, Vector{Atom})
                if length(N) === 0
                    @warn "Terminal N was not found on $segment"
                elseif length(N) > 1
                    @warn "Multiple candidates identified for terminal N on segment $segment"
                end
            else
                if length(N.bonds) > 2
                    Peptides.cap!(template, SerialSelection{Residue}(N_terminal_template.id, :id), skip_C_terminal = true)
                else
                    # Case already capped
                end
            end
        end

        # C-terminal
        C_terminal_pose     = DownstreamTerminalSelection{Residue}()(segment, gather = true)[1]
        C_terminal_template = DownstreamTerminalSelection{Residue}()(template, gather = true)[1]
        if !(C_terminal_pose in residues)
            @info "C-terminal not identified as being selected: skipping cap check on this residue."
        else
            C = identify_c_terminal(segment, supress_warn = true) # Uses multiple criteria to identify
            if isa(C, Vector{Atom})
                if length(C) === 0
                    @warn "Terminal C was not found on $segment"
                elseif length(C) > 1
                    @warn "Multiple candidates identified for terminal C on segment $segment"
                end
            else
                if length(C.bonds) > 2
                    Peptides.cap!(template, SerialSelection{Residue}(C_terminal_template.id, :id), skip_N_terminal = true)
                else
                    # Case already capped
                end
            end
        end
        
        t_atoms = ProtoSyn.travel_graph(N_terminal_template[1], search_algorithm = ProtoSyn.Peptides.IUPAC)
        p_atoms = ProtoSyn.travel_graph(N_terminal_pose[1], search_algorithm = ProtoSyn.Peptides.IUPAC)
        p_selected_atoms = reduce(vcat, [r.items for r in residues])
        p_atoms = p_atoms[p_atoms .∈ [p_selected_atoms]]

        @assert length(t_atoms) === length(p_selected_atoms) "Template atoms ($(length(t_atoms))) and pose atoms ($(length(p_selected_atoms))) don't match in number."
        if all([a.symbol for a in t_atoms] .=== [a.symbol for a in p_selected_atoms])
            for (t_atom, p_atom) in zip(t_atoms, p_atoms)
                ProtoSyn.rename!(p_atom, t_atom.name, force_rename = force_rename)
            end
        end
        
        # In case the direct approach wasn't successful
        @info "Possible tautomers identified, assigning default names residue by residue ..."
        pose_temp_res = zip(residues, eachresidue(template.graph))
        for (pose_residue, template_residue) in pose_temp_res
            # println("\n$pose_residue - $template_residue")
            
            pose_res_name = string(ProtoSyn.three_2_one[pose_residue.name.content])
            ind_pose_res  = copy(pose_residue) # Removes inter-residue connections
            pose_N        = ProtoSyn.identify_atom_by_bonding_pattern(ind_pose_res, ["N", "C", "C", "O"])
            pose_atoms    = ProtoSyn.travel_graph(pose_N, search_algorithm = ProtoSyn.BFS)
            pose_graph    = [a.symbol for a in pose_atoms]
            # _pose_graph    = [a.name for a in pose_atoms]
            # println("Residue $pose_residue (Graph: $pose_graph)")
            
            if isa(grammar.variables[pose_res_name], Tautomer)
                found_match = false
                for tautomer in grammar.variables[pose_res_name].list
                    temp_N     = ProtoSyn.identify_atom_by_bonding_pattern(tautomer.graph[1], ["N", "C", "C", "O"])
                    temp_atoms = ProtoSyn.travel_graph(temp_N, search_algorithm = ProtoSyn.BFS)
                    temp_graph = [a.symbol for a in temp_atoms]
                    # _temp_graph = [a.name for a in temp_atoms]
                    
                    # println(_pose_graph)
                    # println(_temp_graph)
                    if all(pose_graph .=== temp_graph)
                        # println("Match found.")
                        for (t_atom, p_atom) in zip(temp_atoms, pose_atoms)
                            pri = pose_residue.items
                            actual_pose_a = findfirst((a)->a.index == p_atom.index, pri)
                            ProtoSyn.rename!(pri[actual_pose_a], t_atom.name, force_rename = force_rename)
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
                
                @assert length(temp_atoms) === length(pose_atoms) "Number of atoms in $ind_pose_res and corresponding template don't match: $(length(pose_atoms)) in Pose <-> $(length(temp_atoms)) in template.\n$(repeat(" ", 23))Does the Pose residue have any missing or added atoms? Consider ignoring this residue by setting an AbstractSelection."
                if all(pose_graph .=== temp_graph)
                    for (t_atom, p_atom) in zip(temp_atoms, pose_atoms)
                        pri = pose_residue.items
                        actual_pose_a = findfirst((a)->a.index == p_atom.index, pri)
                        # println("renaming atom $(pri[actual_pose_a]) to $(t_atom.name)")
                        try
                            ProtoSyn.rename!(pri[actual_pose_a], t_atom.name, force_rename = force_rename)
                        catch AssertionError
                            error("Multiple Atom instances with the same name identified. Suggested fix: run `assign_default_atom_names!` with `force_rename` set to `true`.")
                        end
                    end
                else
                    println("    Graph: $([x.name for x in pose_atoms])")
                    println(" Template: $([x.name for x in temp_atoms])")
                    @warn "Available template for $pose_residue doesn't match the graph.\n    Graph: $pose_graph\n Template: $temp_graph"
                end
            end
        end
    end
        
    return pose
end