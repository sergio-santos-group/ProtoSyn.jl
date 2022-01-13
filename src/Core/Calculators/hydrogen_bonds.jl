module HydrogenBonds

    using ProtoSyn
    using LinearAlgebra
    using ProtoSyn.Calculators: EnergyFunctionComponent

    """
    # TODO
    """
    struct HydrogenBondPair
        charged::Atom
        base::Atom
    end

    function Base.show(io::IO, hbp::HydrogenBondPair, level_code::Opt{ProtoSyn.LevelCode} = nothing)
        level_code = level_code === nothing ? ProtoSyn.LevelCode() : level_code
        lead       = ProtoSyn.get_lead(level_code)
    
        b = "$(hbp.base.parent.name):$(hbp.base.parent.index):$(hbp.base.name)"
        c = "$(hbp.charged.parent.name):$(hbp.charged.parent.index):$(hbp.charged.name)"
        println(io, lead*"$b - $c")
    end
    
    """
    # TODO
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

    If a hydrogen bond network is provided (instead of a function), this
    function assumes all provided pairs are within the selection (no check is
    made).

    """
    function calc_hydrogen_bond_network(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; hydrogen_bond_network::Union{HydrogenBondNetwork, Function} = HydrogenBondNetwork(), d0::T = 3.0) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        
        if isa(hydrogen_bond_network, Function)
            hydrogen_bond_network = hydrogen_bond_network(pose, selection)
        end

        ehb = T(0.0)
        for acceptor in hydrogen_bond_network.acceptors
    
            c = pose.state[acceptor.base].t
            o = pose.state[acceptor.charged].t
            voc = c .- o
            doc = norm(voc)
    
            for donor in hydrogen_bond_network.donors
    
                n = pose.state[donor.base].t
                h = pose.state[donor.charged].t
                vhn = n .- h
                dhn = norm(vhn)
    
                vho = o .- h
                voh = .- vho
                dho = norm(vho)
    
                cos1 = dot(vho, vhn) / (dho * dhn)
                cos2 = dot(voc, voh) / (doc * dho)
    
                d = (dho - d0)^2
                c1 = cos1 <= 0.0 ? cos1 : 0.0
                c2 = cos2 <= 0.0 ? cos2 : 0.0
    
                ehb_i = d - c1 * c2
                ehb_i = ehb_i < 0.0 ? ehb_i : 0.0
                ehb += ehb_i
            end
        end
        
        return ehb, nothing
    end

    calc_hydrogen_bond_network(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; hydrogen_bond_network::Union{HydrogenBondNetwork, Function} = HydrogenBondNetwork(), d0::T = 3.0) where {T <: AbstractFloat} = begin
        calc_hydrogen_bond_network(ProtoSyn.acceleration.active, pose, selection, update_forces; hydrogen_bond_network = hydrogen_bond_network, d0 = d0)
    end

    """
    # TODO
    H connected to only 1 N
    O connected to only 1 C
    """
    function generate_hydrogen_bond_network(pose::Pose, selection::Opt{AbstractSelection})

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
    """
    function fixate_hydrogen_bond_network!(efc::EnergyFunctionComponent, pose::Pose)
        efc.settings[:hydrogen_bond_network] = efc.settings[:hydrogen_bond_network](pose)
    end


    """
    # TODO
    """
    function get_default_hydrogen_bond_network_comp(;α::T = 1.0) where {T <: AbstractFloat}
        return EnergyFunctionComponent(
            "Hydrogen_Bond_Network",
            calc_hydrogen_bond_network,
            nothing,
            Dict{Symbol, Any}(
                :hydrogen_bond_network => generate_hydrogen_bond_network,
                :d0 => 3.0),
            α,
            false)
    end
end