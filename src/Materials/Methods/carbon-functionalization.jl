using Random
using StatsBase

const edge           = BondCountSelection(3, <)
const edge_or_center = BondCountSelection(4, <)

known_functional_groups = Dict{String, AbstractSelection}(
    "ine" => edge,
    "nyl" => edge,
    "eth" => edge
)

"""
# TODO: Documentation
"""
function functionalize(pose::Pose, functional_groups::Dict{Fragment, T}; normalize_frequencies::Bool = false) where {T <: AbstractFloat}

    # Initial verifications
    for (frag, perc) in functional_groups
        @assert 0.0 <= perc <= 1.0 "Relative frequency (%) of functional group $(frag.graph.name) must be between 0.0 and 1.0."
    end

    S = sum(collect(values(functional_groups)))
    if !normalize_frequencies
        @assert S <= 1.0 "The total requested frequency (%) of functional groups exceeds 1.0. Consider setting the `normalize_frequencies` flag to `true` if you wish to saturate the given pose with all the requested functional groups."
    else
        for (frag, perc) in functional_groups
            functional_groups[frag] = perc / S
        end
    end

    init_N_carbons = count((ChargeSelection(0.0) & BondCountSelection(3, <=))(pose))
    functional_groups_n = Dict{Fragment, Int}()
    for (frag, perc) in functional_groups
        functional_groups_n[frag] = floor(Int, perc * init_N_carbons)
    end

    return functionalize(pose, functional_groups_n)
end


"""
# TODO: Documentation
Warning: Based on charge
"""
function functionalize(pose::Pose, functional_groups::Dict{Fragment, Int})

    init_N_carbons = count((ChargeSelection(0.0) & BondCountSelection(3, <=))(pose))
    for (frag, n) in functional_groups
        @info "$(frag.graph.name) : $n/$init_N_carbons"
    end

    # Convert N functional groups into a randomized list
    functional_group_id = Dict{Int, Fragment}()
    functional_group_list = Vector{Int}()
    for (i, (frag, n)) in enumerate(functional_groups)
        functional_group_id[i] = frag
        n_functional_groups    = repeat([i], outer = n)
        functional_group_list  = vcat(functional_group_list, n_functional_groups)
    end
    shuffle!(functional_group_list)

    # Consume list and add functional groups to pose
    available = ChargeSelection(0.0)
    while length(functional_group_list) > 0
        fcn_id = pop!(functional_group_list)
        fcn = functional_group_id[fcn_id]
        ProtoSyn.write(pose, "current.pdb")
        if !(fcn.graph.name in keys(known_functional_groups))
            @warn "Tried to add unknown function group: $fcn"
        else
            fcn_sele    = known_functional_groups[fcn.graph.name]
            random_atom = sample((fcn_sele & available)(pose, gather = true))
            random_atom = pose.graph[1, 1, 7]
            println("\nAdding $(fcn.graph.name) to $random_atom")

            # Replace the parenthood on the randomly selected atom
            random_parent = random_atom.bonds[1]
            if !isparent(random_parent, random_atom)
                ProtoSyn.popparent!(random_atom)
                ProtoSyn.setparent!(random_atom, random_parent)
                reindex(pose.graph, set_ascendents = true)
                reindex(pose.state)
                random_atom.ascendents = (random_atom.index, random_parent.index, random_parent.parent.index, random_parent.parent.parent.index)
                sync!(pose)
                println("Synching")
                ProtoSyn.request_c2i!(pose.state, all = true)
                sync!(pose)
            end

            ProtoSyn.replace_by_fragment!(pose, random_atom, fcn, 
                remove_downstream_graph = false,
                spread_excess_charge    = true)

            return pose
        end
    end

end
