module HydrogenBonds

    using ProtoSyn
    using LinearAlgebra
    using ProtoSyn.Calculators: EnergyFunctionComponent

    """
        HydrogenBondPair(charged::Atom, base::Atom)

    Define a new [`HydrogenBondPair`](@ref). An [`HydrogenBondPair`](@ref) is a
    set of two bonded [`Atom`](@ref) instances, the `charged` and `base` atom.
    As an example, in a carbonyl group (-C=O), the carbon is the `base` and the
    oxygen is the `charged` atom. It's expected that this
    [`HydrogenBondPair`](@ref) is involved in hydrogen bonding interactions.

    # See also
    [`HydrogenBondNetwork`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.HydrogenBonds.HydrogenBondPair(pose.graph[1, 1, "C"], pose.graph[1, 1, "O"])
    MET:1:O - MET:1:C
    ```
    """
    struct HydrogenBondPair
        charged::Atom
        base::Atom
    end

    function Base.show(io::IO, hbp::HydrogenBondPair, level_code::Opt{ProtoSyn.LevelCode} = nothing)
        level_code = level_code === nothing ? ProtoSyn.LevelCode() : level_code
        lead       = ProtoSyn.get_lead(level_code)
    
        b = "$(hbp.base.container.name):$(hbp.base.container.index):$(hbp.base.name)"
        c = "$(hbp.charged.container.name):$(hbp.charged.container.index):$(hbp.charged.name)"
        println(io, lead*"$b - $c")
    end
    

    """
        HydrogenBondNetwork(donors::Vector{HydrogenBondPair}, acceptors::Vector{HydrogenBondPair})

    Define a new [`HydrogenBondNetwork`](@ref), a set of `donors`
    [`HydrogenBondPair`](@ref) instances (strongly electronegative
    [`Atom`](@ref) instances, such as N, O or F) and `acceptors`
    [`HydrogenBondPair`](@ref) instances (electronegative [`Atom`](@ref)
    instances with a lone electron pair).

        HydrogenBondNetwork()

    Define an empty [`HydrogenBondNetwork`](@ref).

    # See also
    [`calc_hydrogen_bond_network`](@ref) [`generate_hydrogen_bond_network`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.HydrogenBonds.generate_hydrogen_bond_network(pose)
    Donors: 131 | Acceptors: 104

    julia> ProtoSyn.Calculators.HydrogenBonds.HydrogenBondNetwork()
    Donors: 0 | Acceptors: 0
    ```
    """
    struct HydrogenBondNetwork
        donors::Vector{HydrogenBondPair}
        acceptors::Vector{HydrogenBondPair}
    end

    function Base.show(io::IO, hbn::HydrogenBondNetwork, level_code::Opt{LevelCode} = nothing)
        println(io, "Donors: $(length(hbn.donors)) | Acceptors: $(length(hbn.acceptors))")
    end
    
    HydrogenBondNetwork() = begin
        HydrogenBondNetwork(Vector{HydrogenBondPair}(), Vector{HydrogenBondPair}())
    end


    """
        generate_hydrogen_bond_network(pose::Pose, selection::Opt{AbstractSelection} = nothing)

    Attempts to predict the [`HydrogenBondNetwork`](@ref) of a given
    [`Pose`](@ref) `pose` (restricted to the selected region provided by the
    `AbstractSelection` `selection`). The following simple criteria are used:
    * Donors are hydrogen (H) [`Atom`](@ref) instances connected to a single nitrogen (N) [`Atom`](@ref) instance;
    * Acceptors are oxygen (O) [`Atom`](@ref) instances connected to a single carbon (C) [`Atom`](@ref) instance.

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.HydrogenBonds.generate_hydrogen_bond_network(pose)
    Donors: 131 | Acceptors: 104
    ```
    """
    function generate_hydrogen_bond_network(pose::Pose, selection::Opt{AbstractSelection} = nothing)

        if selection === nothing
            sele = TrueSelection{Atom}()
        else
            sele = ProtoSyn.promote(selection, Atom)
        end

        hbn  = HydrogenBondNetwork()

        for atom in sele(pose, gather = true)
            
            if (atom.symbol === "H") && (length(atom.bonds) === 1) && (atom.bonds[1].symbol === "N")
                push!(hbn.donors, HydrogenBondPair(atom, atom.bonds[1]))
            end

            if (atom.symbol === "O") && (length(atom.bonds) === 1) && (atom.bonds[1].symbol === "C")
                push!(hbn.acceptors, HydrogenBondPair(atom, atom.bonds[1]))
            end
        end

        return hbn
    end


    """
        calc_hydrogen_bond_network([::Type{A}], pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; [hydrogen_bond_network::Union{HydrogenBondNetwork, Function} = HydrogenBondNetwork()], [potential::Function = (x; qi = 0.0, qj = 0.0) -> 0.0]) where {A <: ProtoSyn.AbstractAccelerationType}
        
    Calculate the [`Pose`](@ref) `pose` hydrogen bond energy according to the
    given `potential` function (make sure the [`Pose`](@ref) `pose` is synched,
    see [`sync!`](@ref)). ProtoSyn loops through the provided
    [`HydrogenBondNetwork`](@ref) `hydrogen_bond_network` and applied the
    `potential` to each pair of charged [`Atom`](@ref) instances (based on the
    inter-atomic distance). Besides this component, `calc_hydrogen_bond_network`
    adds a geometric potential based on the angle between the charged
    [`Atom`](@ref) instances and each base (rewards 180º conformations).
    Optionally, `hydrogen_bond_network` can be a `Function`, in which case a new
    [`HydrogenBondNetwork`](@ref) is calculated. Such a custom function should
    have the following signature:

    # ```my_hydrogen_bond_predictor(pose::Pose, selection::Opt{AbstractSelection} = nothing)```

    If provided, an `AbstractSelection` `selection` limits the selected
    [`Atom`](@ref) instances considered for [`HydrogenBondNetwork`](@ref)
    prediction (if `hydrogen_bond_network` is a `Function`, otherwise the
    provided [`HydrogenBondNetwork`](@ref) is static and `selection` has no
    effect). An optional parameter `A` (`Type{<: AbstractAccelerationType}`) can
    be provided, stating the acceleration type used to calculate this energetic
    contribution (See [ProtoSyn acceleration types](@ref)). Uses
    `ProtoSyn.acceleration.active` by default. Note that, depending on the
    `potential` employed, atomic charges may be required (See
    [`assign_default_charges`](@ref ProtoSyn.Calculators.Electrostatics.assign_default_charges!),
    for example).

    # See also
    [`get_default_hydrogen_bond_network`](@ref)

    # Examples
    ```
    julia> e, f, hb_pairs = ProtoSyn.Calculators.HydrogenBonds.calc_hydrogen_bond_network(pose, nothing, false, hydrogen_bond_network = hbn, potential = potential)
    (-9.329167508668307, nothing, (...))

    julia> hb_pairs
    122-element Vector{Tuple{Atom, Atom}}:
    (Atom{/2a3d:41940/A:1/MET:1/O:19}, Atom{/2a3d:41940/A:1/TRP:4/H:39})
    (Atom{/2a3d:41940/A:1/MET:1/O:19}, Atom{/2a3d:41940/A:1/ALA:5/H:63})
    (Atom{/2a3d:41940/A:1/GLY:2/O:26}, Atom{/2a3d:41940/A:1/ALA:5/H:63})
    (...)
    ```
    """
    function calc_hydrogen_bond_network(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; hydrogen_bond_network::Union{HydrogenBondNetwork, Function} = HydrogenBondNetwork(), potential::Function = (x; qi = 0.0, qj = 0.0) -> 0.0) where {A <: ProtoSyn.AbstractAccelerationType}
        
        if isa(hydrogen_bond_network, Function)
            hydrogen_bond_network = hydrogen_bond_network(pose, selection)
        end

        T = eltype(pose.state)
        ehb = T(0.0)
        hb_list = Vector{Tuple{Atom, Atom}}()
        for acceptor in hydrogen_bond_network.acceptors
    
            c = pose.state[acceptor.base].t
            o = pose.state[acceptor.charged].t
            qi = pose.state[acceptor.charged].δ
            voc = c .- o
            doc = norm(voc)
    
            for donor in hydrogen_bond_network.donors
    
                n = pose.state[donor.base].t
                h = pose.state[donor.charged].t
                qj = pose.state[donor.charged].δ
                vhn = n .- h
                dhn = norm(vhn)
    
                vho = o .- h
                voh = .- vho
                dho = norm(vho)
    
                cos1 = dot(vho, vhn) / (dho * dhn)
                cos2 = dot(voc, voh) / (doc * dho)
    
                d = potential(dho, qi = qi, qj = qj)
                c1 = cos1 <= 0.0 ? cos1 : 0.0
                c2 = cos2 <= 0.0 ? cos2 : 0.0
    
                ehb_i = d * c1 * c2
                ehb_i = ehb_i < 0.0 ? ehb_i : 0.0
                ehb += ehb_i
                
                if ehb_i < 0.0
                    push!(hb_list, (acceptor.charged, donor.charged))
                end
            end
        end
        
        return ehb, nothing, hb_list
    end

    calc_hydrogen_bond_network(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; hydrogen_bond_network::Union{HydrogenBondNetwork, Function} = HydrogenBondNetwork(), potential::Function = (x; qi = 0.0, qj = 0.0) -> 0.0) = begin
        calc_hydrogen_bond_network(ProtoSyn.acceleration.active, pose, selection, update_forces; hydrogen_bond_network = hydrogen_bond_network, potential = potential)
    end


    """
        get_default_hydrogen_bond_network(;[α::T = 1.0]) where {T <: AbstractFloat}

    Return the default hydrogen bond [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent). `α` sets
    the component weight (on an
    [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance). This
    component employs the [`calc_hydrogen_bond_network`](@ref) method, therefore
    defining a [`Pose`](@ref) energy based on a given potential function
    multiplied by a geometric angle component for each pair defined in an
    [`HydrogenBondNetwork`](@ref). By default, this
    [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) uses the
    [`generate_hydrogen_bond_network`](@ref) function to predict a new
    [`HydrogenBondNetwork`](@ref) each call. Consider setting
    `efc.settings[:hydrogen_bond_network]` as an [`HydrogenBondNetwork`](@ref)
    to employ a static list of interacting [`Atom`](@ref) instances (for
    improved performance).

    # See also
    [`fixate_hydrogen_bond_network!`](@ref)

    # Settings
    * `hydrogen_bond_network::Union{HydrogenBondNetwork, Function}` - Defines either the [`HydrogenBondNetwork`](@ref) predictor function or static [`HydrogenBondNetwork`](@ref);
    * `potential::Function` - Define the potential to apply (calculates energy and force value from inter-atomic distance - and optionally atomic charges `qi` & `qj`);

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.HydrogenBonds.get_default_hydrogen_bond_network()
    🞧  Energy Function Component:
    +---------------------------------------------------+
    | Name           | Hydrogen_Bonds                   |
    | Alpha (α)      | 1.0                              |
    | Update forces  | false                            |
    | Calculator     | calc_hydrogen_bond_network       |
    +---------------------------------------------------+
     |    +----------------------------------------------------------------------------------+
     ├──  ● Settings                      | Value                                            |
     |    +----------------------------------------------------------------------------------+
     |    | potential                     | bump_potential_charges                           |
     |    | hydrogen_bond_network         | generate_hydrogen_bond_network                   |
     |    +----------------------------------------------------------------------------------+
     |    
     └──  ○  Selection: nothing
    ```
    """
    function get_default_hydrogen_bond_network(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Hydrogen_Bonds",
            calc_hydrogen_bond_network,
            nothing,
            Dict{Symbol, Any}(
                :hydrogen_bond_network => generate_hydrogen_bond_network,
                :potential => ProtoSyn.Calculators.get_bump_potential_charges(c = 3.0, r = 1.5)),
            α,
            false)
    end


    """
        fixate_hydrogen_bond_network!(efc::EnergyFunctionComponent, pose::Pose)

    If the given [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) `efc` is an Hydrogen Bonds
    [`EnergyFunctionComponent`](@ref ProtoSyn.Calculators.EnergyFunctionComponent) (and the `:hydrogen_bond_network` setting
    is a `Function`), calculate a new [`HydrogenBondNetwork`](@ref) and apply it
    as a static list of interaction [`Atom`](@ref) instances (improved
    performance).

    # See also
    [`get_default_hydrogen_bond_network`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Calculators.HydrogenBonds.fixate_hydrogen_bond_network!(efc, pose)
    Donors: 131 | Acceptors: 104
    ```
    """
    function fixate_hydrogen_bond_network!(efc::EnergyFunctionComponent, pose::Pose)
        @assert :hydrogen_bond_network in keys(efc.settings) ":hydrogen_bond_network setting not found in the provided EnergyFunctionComponent. Are you sure this is an Hydrogen Bonds EnergyFunctionComponent?"
        @assert isa(efc.settings[:hydrogen_bond_network], Function) ":hydrogen_bond_network doesn't seem to be a function. Perhaps it is already set as an HydrogenBondNetwork?"
        efc.settings[:hydrogen_bond_network] = efc.settings[:hydrogen_bond_network](pose)
    end
end