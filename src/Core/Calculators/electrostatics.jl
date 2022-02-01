module Electrostatics

    using ProtoSyn
    using ProtoSyn.Calculators: EnergyFunctionComponent, VerletList, MaskMap
    
    # TODO: Documentation

    function assign_acc2_eem_charges!(pose::Pose, filename::String)
        data = split(open(f->read(f, String), filename), "\n")[2]
        T = eltype(pose.state)
        charges = map((value) -> parse(T, value), split(data))
        for (i, atom) in enumerate(eachatom(pose.graph))
            pose.state[atom].δ = charges[i]
        end

        return charges
    end    

    # ---

    function calc_coulomb(::Type{A}, pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing) where {A <: ProtoSyn.AbstractAccelerationType}
        cp = ProtoSyn.Calculators.coulomb_potential
        e, f = ProtoSyn.Calculators.apply_potential(A, pose, cp, update_forces, vlist, selection, mask)
        return e, f
    end # function

    calc_coulomb(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool; mask::MaskMap = nothing, vlist::Opt{VerletList} = nothing) = begin
        calc_coulomb(ProtoSyn.acceleration.active, pose, selection, update_forces; mask = mask, vlist = vlist)
    end

    function get_default_coulomb(;α::T = 1.0) where {T <: AbstractFloat}
        _mask = ProtoSyn.Calculators.get_bonded_mask()
        return EnergyFunctionComponent(
            "Coulomb",
            calc_coulomb,
            nothing,
            Dict{Symbol, Any}(:mask => _mask, :vlist => nothing),
            α,
            true)
    end

end