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
    # TODO
    H connected to only 1 N
    O connected to only 1 C
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

    generate_hydrogen_bond_network(pose::Pose) = generate_hydrogen_bond_network(pose, nothing)

    
    """
    # TODO

    If a hydrogen bond network is provided (instead of a function), this
    function assumes all provided pairs are within the selection (no check is
    made).

    """
    function calc_hydrogen_bond_network(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; hydrogen_bond_network::Union{HydrogenBondNetwork, Function} = HydrogenBondNetwork(), potential::Function = (x) -> 0.0) where {A <: ProtoSyn.AbstractAccelerationType}
        
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

    calc_hydrogen_bond_network(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; hydrogen_bond_network::Union{HydrogenBondNetwork, Function} = HydrogenBondNetwork(), potential::Function = (x) -> 0.0) = begin
        calc_hydrogen_bond_network(ProtoSyn.acceleration.active, pose, selection, update_forces; hydrogen_bond_network = hydrogen_bond_network, potential = potential)
    end

    """
    # TODO
    """
    function fixate_hydrogen_bond_network!(efc::EnergyFunctionComponent, pose::Pose)
        efc.settings[:hydrogen_bond_network] = efc.settings[:hydrogen_bond_network](pose)
    end


    """
    # TODO
    """
    function get_default_hydrogen_bond_network_comp(;α::T = 1.0) where {T <: AbstractFloat}
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
end