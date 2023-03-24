"""
    find_tautomer(tautomer::Tautomer, target::Residue)

Given a `target` [`Residue`](@ref), search the provided Tautomer `tautomer` list
for the corresponding template [`Residue`](@ref), based on the [Graph](@ref graph-types)
(employs the [`travel_graph`](@ref ProtoSyn.travel_graph) method).

# Examples
```
julia> tautomer = Peptides.grammar.variables["H"]
Fragment(Segment{/HIE:22535}, State{Float64}:
 Size: 17
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)(And 1 other tautomer(s) available.)

julia> ProtoSyn.find_tautomer(tautomer, pose.graph[1][72])
Fragment(Segment{/HID:3247}, State{Float64}:
 Size: 17
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function find_tautomer(tautomer::Tautomer, target::Residue)
    target_res   = copy(target)
    target_atoms = ProtoSyn.travel_graph(target_res[1], search_algorithm = ProtoSyn.BFS)
    target_graph = [a.symbol for a in target_atoms]
    _target_graph = [a.name for a in target_atoms]
    # println("Target")
    # println(_target_graph)
    for template_residue in tautomer.list
        template_res   = copy(template_residue.graph[1])
        template_atoms = ProtoSyn.travel_graph(template_res[1], search_algorithm = ProtoSyn.BFS)
        template_graph = [a.symbol for a in template_atoms]
        _template_graph = [a.name for a in template_atoms]

        length(target_graph) !== length(template_graph) && continue

        # println(template_residue.graph[1].name)
        # println(_template_graph)

        match = all(target_graph .=== template_graph)
        match && return template_residue
    end

    return nothing
end