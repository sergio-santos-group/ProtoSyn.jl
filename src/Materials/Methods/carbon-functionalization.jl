using Random
using StatsBase
using ProtoSyn.Units

base           = as"C"
edge           = base & BondCountSelection(3, <)
edge_or_center = base & BondCountSelection(4, <)
center         = base & BondCountSelection(3)

known_functional_groups = Dict{String, AbstractSelection}(
    "ine" => edge,
    "nyl" => edge,
    "eth" => edge,
    "grn" => center & !BondedToSelection(as"N")
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
    while length(functional_group_list) > 0
        fcn_id = pop!(functional_group_list)
        fcn = functional_group_id[fcn_id]
        ProtoSyn.write(pose, "current.pdb")
        if !(fcn.graph.name in keys(known_functional_groups))
            @warn "Tried to add unknown function group: $fcn"
        else
            fcn_sele    = known_functional_groups[fcn.graph.name]

            mask = (fcn_sele)(pose)
            count(mask.content) === 0 && begin
                @warn "Tried to add $(fnc.graph.name) functional group to the given Pose, but no matching anchoring Atom was available."
                continue
            end

            random_atom = StatsBase.sample(ProtoSyn.gather(mask, pose.graph))
            
            @info "Adding $(fcn.graph.name) to $random_atom"
            add_functionalization(pose, fcn, random_atom)
        end
    end
end


"""
Assumes z = 0
"""
function add_functionalization(pose::Pose, fcn::Fragment, atom::Atom)

    fcn = copy(fcn) # So that any changes don't apply to the template group

    # Measure the dihedral between the selected atom parents and a random bond
    # (In this case, the first bond found is picked). If this dihedral is π, the
    # functional group's third atom (if existent) is set to rotate by 180°.
    t  = [atom.parent.parent, atom.parent, atom, [a for a in atom.bonds if a !== atom.parent][1]]
    d = ProtoSyn.dihedral(map(a -> pose.state[a], t)...)
    if d ≈ π && ProtoSyn.count_atoms(fcn.graph) > 2
        fcn.state[fcn.graph[1, 3]].ϕ += 180°
    end

    ProtoSyn.replace_by_fragment!(pose, atom, fcn, 
        remove_downstream_graph = false,
        spread_excess_charge    = true)
end