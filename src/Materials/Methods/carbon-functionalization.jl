using Random
using Printf
using StatsBase
using ProgressMeter
using LinearAlgebra
using ProtoSyn.Units

base           = as"C" & ChargeSelection(0.0)
edge           = base & BondCountSelection(3, <)
edge_or_center = base & BondCountSelection(4, <)
center         = base & BondCountSelection(3)

"""
    functionalize!(pose::Pose, functional_groups::Dict{Fragment, T}; [normalize_frequencies::Bool = false], [attempt_minimization::Bool = true]) where {T <: AbstractFloat}

Add N functional groups to the given [`Pose`](@ref) `pose`. The number of added
functional groups is given by the `functional groups` dictionary, mapping
[`Fragment`](@ref) instances (for example, from the default
`ProtoSyn.modification_grammar`) to a percentage of achoring points available.
By default, anchoring points are non-charged carbon [`Atom`](@ref) instances.
If the `normalize_frequencies` flag is set to `true` (`false`, by default),
ProtoSyn normalizes the input `functional_groups` dictionary so that all
possible non-charged carbon [`Atom`](@ref) instances are functionalized (in
proportional percentages to the input dictionary). Note that this method expects
correctly ordered [`Segment`](@ref) instances (1 is bottom, :end is top). If
`attempt_minimization` is set to `true` (is, by default), a quick Monte Carlo
simulation attempts to minimize the prevalence of inter-atomic clashes.

!!! ukw "Note:"
    The reason only non-charged carbon [`Atom`](@ref) instances are considered for functionalization is because ProtoSyn automatically assigns a charge when adding a functional group. Therefore, only non-charged [`Atom`](@ref) instances are left "open" for functionalization.


---


    functionalize!(pose::Pose, functional_groups::Dict{Fragment, Int}; [attempt_minimization::Bool = true])

In an alternative syntax, the `functional_groups` fictionary directly maps
[`Fragment`](@ref) instances to the actual number of desired functional groups
to add. If `attempt_minimization` is set to `true` (is, by default), a quick
Monte Carlo simulation attempts to minimize the prevalence of inter-atomic
clashes.

# See also
[`add_functionalization!`](@ref)

# Examples
```
julia> ProtoSyn.Materials.functionalize!(pose, Dict(ProtoSyn.modification_grammar.variables["XYL"] => 10))

julia> pose
Pose{Topology}(Topology{/CRV:42474}, State{Float64}:
 Size: 380
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function functionalize!(pose::Pose, functional_groups::Dict{Fragment, T}; normalize_frequencies::Bool = false, attempt_minimization::Bool = true) where {T <: AbstractFloat}

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

    return functionalize!(pose, functional_groups_n,
        attempt_minimization = attempt_minimization)
end

function functionalize!(pose::Pose, functional_groups::Dict{Fragment, Int}; attempt_minimization::Bool = true)

    L_layer = SerialSelection{Segment}(ProtoSyn.count_segments(pose.graph), :id)

    # Start fcns added dict
    fcns_added = Dict{String, Int}()
    for key in keys(functional_groups)
        fcns_added[key.graph.name] = 0
    end

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
    p = Progress(length(functional_group_list), dt = 0.2, desc = "Functionalizing carbon ...", color = :gray)
    while length(functional_group_list) > 0

        next!(p)

        fcn_id = pop!(functional_group_list)
        fcn = functional_group_id[fcn_id]
        if !(fcn.graph.name in keys(known_functional_groups))
            @warn "Tried to add unknown function group: $fcn"
        else
            fcn_sele = known_functional_groups[fcn.graph.name]

            mask = (fcn_sele)(pose)
            count(mask.content) === 0 && begin
                @warn "Tried to add $(fcn.graph.name) functional group to the given Pose, but no matching anchoring Atom was available."
                continue
            end

            random_atom = StatsBase.sample(ProtoSyn.gather(mask, pose.graph))
            
            @info "Adding $(fcn.graph.name) to $random_atom"
            fcn_added = add_functionalization!(pose, fcn, random_atom,
                attempt_minimization = attempt_minimization)

            if fcn_added !== nothing
                fcns_added[fcn.graph.name] += 1
            end
        end
    end

    # Adjust residual partial charge (should be close to 0.0)
    spread_charges!(pose)

    return pose, fcns_added
end

function spread_charges!(pose::Pose)
    c = sum([a.δ for a in pose.state.items[4:end]])
    Δc = c / ProtoSyn.count_atoms(pose.graph)
    for atom in pose.state.items[4:end]
        atom.δ -= Δc
    end
end


"""
    add_functionalization!(pose::Pose, fcn::Fragment, atom::Atom)

Add a single functional group `fcn` (a [`Fragment`](@ref) instance) to the given
[`Atom`](@ref) instance `atom` in the [`Pose`](@ref) `pose`. The [`Atom`](@ref)
instance `atom` is replaced (using the
[`replace_by_fragment!`](@ref ProtoSyn.replace_by_fragment!) method).

# See also
[`functionalize!`](@ref)

# Examples
```
julia> ProtoSyn.Materials.add_functionalization!(pose, ProtoSyn.modification_grammar.variables["XYL"], pose.graph[1, 1, 16])

julia> pose
Pose{Topology}(Topology{/CRV:14425}, State{Float64}:
 Size: 353
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function add_functionalization!(pose::Pose, fcn::Fragment, atom::Atom; attempt_minimization::Bool = true)

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
        pose.state[root].t[3:3]               .= (layer_ID-1) * d
        pose.state[root.parent].t[3:3]        .= (layer_ID-1) * d
        pose.state[root.parent.parent].t[3:3] .= (layer_ID-1) * d

        if layer_ID % 2 === 0 # set multiple layer offset
            v2 = [r/2, cos(60°)*r, 0.0]
            pose.state[root].t[1:3]               .-= v2
            pose.state[root.parent].t[1:3]        .-= v2
            pose.state[root.parent.parent].t[1:3] .-= v2
        end

        # pose.state.i2c && sync!(pose)
        ProtoSyn.request_c2i!(pose.state, all = true)
        sync!(pose)
    end

    fcn = copy(fcn) # So that any changes don't apply to the template group
    child_of_root = atom in ProtoSyn.root(pose.graph).children

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
                @info "No suitable placement for functional group was found."
                return nothing
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

            if child_of_root
                @info "Child of root!"
                # The next values can't be 0.0 because of trignometry errors
                fcn.state[fcn.graph[1, 3]].θ = 0.001°
                for child in fcn.graph[1, 3].children
                    fcn.state[child].ϕ += 180.001°
                end
            end

            d = ProtoSyn.dihedral(map(a -> pose.state[a], t)...)
            if ((d ≈ π) | (d ≈ -π))
                @info " | 180° parenthood"
                fcn.state[fcn.graph[1, 3]].ϕ += 180°
            end
        end
    end

    _id = atom.id
    ProtoSyn.replace_by_fragment!(pose, atom, fcn, 
        remove_downstream_graph = false,
        spread_excess_charge    = true)

    # Neat indexation
    reindex(pose)
    ProtoSyn.request_i2c!(pose.state, all = false)
    sync!(pose)

    # Block "close-by" points (marked with +0.001 charge)
    if !center && ProtoSyn.count_atoms(fcn.graph) > 2
        @info "Blocking atoms ..."
        child_atom = [c for c in atom.children if c.name !== "C"][1]
        c_sele     = SerialSelection{Atom}(child_atom.id, :id)
        n_sele     = (1.6:c_sele) & !BondedToSelection(c_sele) & !c_sele
        for atom in n_sele(pose, gather = true)
            @info " Blocked atom $atom"
            pose.state[atom].δ += 0.001
        end
    end

    # (Optional) Minimize structure
    if attempt_minimization && ProtoSyn.count_atoms(fcn.graph) > 3

        # 2. Minimization
        energy_fcn = ProtoSyn.Calculators.EnergyFunction([
            ProtoSyn.Calculators.Restraints.get_default_all_atom_clash_restraint()
        ])
        energy_fcn.selection = 4.0:SerialSelection{Atom}(_id, :id)
        fcn_names  = [a.name for a in eachatom(fcn.graph)]
        child_atom = [c for c in atom.children if (c.name in fcn_names) && (length(c.children) > 0)][1]
        @info "Atom: $atom\n Atom children: $(atom.children)\n Child atom: $child_atom\n Child atom children: $(child_atom.children)"
        dh_sele    = SerialSelection{Atom}(child_atom.children[1].id, :id)
        sampler    = ProtoSyn.Mutators.DihedralMutator(randn, 1.0, 0.5, dh_sele)
        thermostat = ProtoSyn.Drivers.get_constant_temperature(0.0)
        minimizer  = ProtoSyn.Drivers.MonteCarlo(energy_fcn, sampler, nothing, 10, thermostat)
        minimizer(pose)
    end

    sync!(pose)

    # Recover original root
    move_root_to_layer(1)

    return pose
end


"""
"""
function add_hydrogens!(pose::Pose, res_lib::LGrammar, selection::Opt{AbstractSelection} = nothing)

    hydro = res_lib.variables["HID"]

    atoms = (selection & an"C|C0"r & BondCountSelection(3, <))(pose, gather = true)

    for atom in atoms
        add_functionalization!(pose, hydro, atom; attempt_minimization = false)
        sync!(pose)
    end

    spread_charges!(pose)

    return pose
end

add_hydrogens!(pose::Pose) = add_hydrogens!(pose, ProtoSyn.modification_grammar, nothing)