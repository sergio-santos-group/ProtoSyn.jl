module HydrogenBonds

    using ProtoSyn
    using LinearAlgebra

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
    """
    function calc_hydrogen_bond_network(::Type{A}, pose::Pose, update_forces::Bool; hydrogen_bond_network::Union{HydrogenBondNetwork, Function} = HydrogenBondNetwork(), d0::T = 3.0) where {A <: ProtoSyn.AbstractAccelerationType, T <: AbstractFloat}
        
        if isa(hydrogen_bond_network, Function)
            hydrogen_bond_network = hydrogen_bond_network(pose)
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
    
                # d = ((0.5 * ((d0 / dho)^12)) - (1.5 * ((d0 / dho)^10)))
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

    calc_hydrogen_bond_network(pose::Pose, update_forces::Bool; hydrogen_bond_network::Union{HydrogenBondNetwork, Function} = HydrogenBondNetwork(), d0::T = 3.0) where {T <: AbstractFloat} = begin
        calc_hydrogen_bond_network(ProtoSyn.acceleration.active, pose, update_forces; hydrogen_bond_network = hydrogen_bond_network, d0 = d0)
    end
end