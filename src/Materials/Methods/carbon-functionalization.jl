using Random
using StatsBase
using LinearAlgebra
using ProtoSyn.Units
using Printf

base           = as"C" & ChargeSelection(0.0) & !DownstreamTerminalSelection{Atom}() & !UpstreamTerminalSelection{Atom}()
edge           = base & BondCountSelection(3, <)
edge_or_center = base & BondCountSelection(4, <)
center         = base & BondCountSelection(3)

"""
# TODO: Documentation
"""
function functionalize!(pose::Pose, functional_groups::Dict{Fragment, T}; normalize_frequencies::Bool = false) where {T <: AbstractFloat}

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

    return functionalize!(pose, functional_groups_n)
end


"""
# TODO: Documentation
Warning: Based on charge
Warning: Expects correctly ordered segments (1 is bottom, end is top)
"""
function functionalize!(pose::Pose, functional_groups::Dict{Fragment, Int})

    L_layer = SerialSelection{Segment}(ProtoSyn.count_segments(pose.graph), :id)

    known_functional_groups = Dict{String, AbstractSelection}(
        "ine" => edge,                                             # Pyridine
        "nyl" => edge,                                             # Carbonyl
        "eth" => edge,                                             # Ether
        "grn" => center & !BondedToSelection(as"N"),               # Graphitic-N
        "xyl" => ((center & (sid"1" | L_layer)) | edge),           # Carboxyl
        "amn" => (center & (sid"1" | L_layer)) | edge,             # Amine-N
        "hyl" => (center & (sid"1" | L_layer)) | edge,             # Hydroxyl
        "oxn" => edge,                                             # Oxidized-N
    )

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
    _count = 0
    while length(functional_group_list) > 0

        for atom in eachatom(pose.graph)
            ϕ = round(rad2deg(pose.state[atom].ϕ), digits = 4)
        end

        fcn_id = pop!(functional_group_list)
        fcn = functional_group_id[fcn_id]
        ProtoSyn.write(pose, "current.pdb")
        if !(fcn.graph.name in keys(known_functional_groups))
            @warn "Tried to add unknown function group: $fcn"
        else
            fcn_sele = known_functional_groups[fcn.graph.name]

            mask = (fcn_sele)(pose)
            count(mask.content) === 0 && begin
                @warn "Tried to add $(fnc.graph.name) functional group to the given Pose, but no matching anchoring Atom was available."
                continue
            end

            random_atom = StatsBase.sample(ProtoSyn.gather(mask, pose.graph))
            
            @info "Adding $(fcn.graph.name) to $random_atom"
            add_functionalization!(pose, fcn, random_atom)
        end
    end
end


"""
Assumes z = 0
"""
function add_functionalization!(pose::Pose, fcn::Fragment, atom::Atom)

    function clockwise(atom::Atom)
        a = collect(pose.state[atom].t[1:2] .- pose.state[atom.parent].t[1:2])
        if length(atom.children) > 0
            b = collect(pose.state[atom.children[1]].t[1:2] .- pose.state[atom].t[1:2])
        else
            b = collect(pose.state[[a for a in atom.bonds if a !== atom.parent][1]].t[1:2] .- pose.state[atom].t[1:2])
        end
        return det(hcat(a, b)') > 0
    end

    function define_functionalization_up_down(up::Bool)

        fcn.state[fcn.graph[1, 3]].ϕ =  90°
        if clockwise(atom)
            fcn.state[fcn.graph[1, 3]].θ = up ? -90° :  90°
        else
            fcn.state[fcn.graph[1, 3]].θ = up ?  90° : -90°
        end

        if length(atom.children) > 0
            t  = [atom.parent.parent, atom.parent, atom, atom.children[1]]
        else
            t = [atom.parent.parent, atom.parent, atom, [a for a in atom.bonds if a !== atom.parent][1]]
        end
        d = ProtoSyn.dihedral(map(a -> pose.state[a], t)...)
        if ((d ≈ π) | (d ≈ -π)) && ProtoSyn.count_atoms(fcn.graph) > 2
            @info " | 180° parenthood"
            fcn.state[fcn.graph[1, 3]].ϕ += 180°
        end
    end

    function move_root_to_layer(layer_ID::Int; r::T = 1.4, d::T = 3.4) where {T <: AbstractFloat}

        root = ProtoSyn.root(pose.graph)
        pose.state[root].t[3:3]               .+= (layer_ID-1) * d
        pose.state[root.parent].t[3:3]        .+= (layer_ID-1) * d
        pose.state[root.parent.parent].t[3:3] .+= (layer_ID-1) * d

        if layer_ID % 2 === 0 # set multiple layer offset
            v2 = [r/2, cos(60°)*r, 0.0]
            pose.state[root].t[1:3]               .+= v2
            pose.state[root.parent].t[1:3]        .+= v2
            pose.state[root.parent.parent].t[1:3] .+= v2
        end

        ProtoSyn.request_c2i!(pose.state, all = true)
        sync!(pose)
    end

    fcn = copy(fcn) # So that any changes don't apply to the template group
    atom_ID = atom.id

    @info "\n---\nAdding functional group \"$(fcn.graph.name)\" to $atom\nPlacement:"

    # Perform adjustments to the template, based on the assigned atom to replace 
    # The functionalization can be assigned to a "center" atom (with 3 bonds) or
    # an "edge" atom (with 2 bonds)
    center = length(atom.bonds) === 3

    _root = ProtoSyn.root(pose.graph)
    root  = [_root, _root.parent, _root.parent.parent]
    if ProtoSyn.count_atoms(fcn.graph) > 2
        if center
            @info " | Center position"
            # In a "center" position, the functional group can be appended either
            # "upwards" or "downwards". In multi-layer carbons, the first layer
            # (bottom one) can only receive "downwards" oriented functional groups,
            # and the last layer (top one) can only receive "upwards" oriented
            # functional groups. Any "middle" layer cannot receive center function
            # groups. If only 1 layer is present, both "upwards" and "downwards"
            # functional groups are accepted. The orientation is picked randomly.
            # The top and bottom layers are inferred from the segment ID (this is
            # expected to be consistent with the output of `generate_carbon`)
            N_layers = ProtoSyn.count_segments(pose.graph)
            layer_ID = atom.container.container.id
            if N_layers === 1
                @info " | Single layer"
                if rand([1, 2]) === 1
                    define_functionalization_up_down(true)
                else
                    define_functionalization_up_down(false)
                end
            elseif layer_ID === 1 # Top layer | "downwards" orientation
                @info " | Top layer"
                define_functionalization_up_down(true)
            elseif layer_ID === N_layers # Bottom layer | "upwards" orientation
                @info " | Bottom layer"
                define_functionalization_up_down(false)
            else
                @error "No suitable placement for functional group was found."
            end
        else
            @info " Edge position"
            # Case on an "edge":
            # Measure the dihedral between the selected atom parents and a random bond
            # (In this case, the first bond found is picked). If this dihedral is π, the
            # functional group's third atom (if existent) is set to rotate by 180°.
            if length(atom.children) > 0
                t  = [atom.parent.parent, atom.parent, atom, atom.children[1]]
                @info " Atom has children. Measuring \"t\" dihedral between atoms $([a.id for a in t])"
            else
                t = [atom.parent.parent, atom.parent, atom, [a for a in atom.bonds if a !== atom.parent][1]]
                @info " Atom has NO children. Measuring \"t\" dihedral between atoms $([a.id for a in t])"
            end
            if any([a in root for a in t])
                @info " Moving root to meet atom on layer $(atom.container.container.index)."
                move_root_to_layer(atom.container.container.index)
            end

            d = ProtoSyn.dihedral(map(a -> pose.state[a], t)...)
            if ((d ≈ π) | (d ≈ -π))
                @info " | 180° parenthood"
                fcn.state[fcn.graph[1, 3]].ϕ += 180°
            end
        end
    end

    ProtoSyn.set_parenthood_as_forces!(pose)
    ProtoSyn.write_forces(pose, "parenthood.txt")
    ProtoSyn.write(pose, "carbon.pdb")

    ProtoSyn.replace_by_fragment!(pose, atom, fcn, 
        remove_downstream_graph = false,
        spread_excess_charge    = true)

    # Block "close-by" points # TODO
    # println(collect(eachatom(pose.graph))[atom.id].children)

    reindex(pose)
    ProtoSyn.request_i2c!(pose.state, all = false)
    sync!(pose)

    ProtoSyn.set_parenthood_as_forces!(pose)
    ProtoSyn.write_forces(pose, "parenthood.txt")
    ProtoSyn.write(pose, "carbon.pdb")
end