# TODO Documentation
function apply!(state::State, rotamer::Rotamer, residue::Residue)
    @assert rotamer.name == residue.name "Tried to apply a $(rotamer.name) rotamer to residue $(residue.name)."
    T = eltype(state)

    chis = Dict{DihedralType, Tuple{T, T}}()
    for (chi, (value, sd)) in rotamer.chis
        _value = (randn() * sd) + value
        setdihedral!(state, chi(residue), _value)
        chis[chi] = (_value, T(0))
    end

    return Rotamer(rotamer.name, chis)
end

function sample(rs::RotamerStack{T}, n::Int = -1) where {T <: AbstractFloat}
    if n <= 0
        return StatsBase.sample(rs.rotamers, rs.weights)
    else
        return StatsBase.sample(rs.rotamers[1:n], rs.weights[1:n])
    end
end

Base.push!(rs::RotamerStack{T}, rotamer::Rotamer{T}, weight::T) where {T <: AbstractFloat} = begin

    push!(rs.rotamers, rotamer)
    rs.weights = Weights(vcat(rs.weights, Weights([weight])))
end

function findnearest(a, x)
    length(a) > 0 || return nothing
    r = searchsorted(a,x)
    length(r) > 0 && return first(r)
    last(r) < 1 && return first(searchsorted(a,a[first(r)]))
    first(r) > length(a) && return first(searchsorted(a,a[last(r)]))
    x-a[last(r)] <= a[first(r)]-x && return first(searchsorted(a,a[last(r)]))
    x-a[last(r)] > a[first(r)]-x && return first(searchsorted(a,a[first(r)]))
end


Base.getindex(rl::R, phi::T, psi::T) where {R <: RotamerLibrary, T <: AbstractFloat} = begin
    phi_index = findnearest(rl.phis, phi)
    psi_index = findnearest(rl.psis, psi)
    return rl.rotamer_stacks[phi_index, psi_index]
end

function get_rotamer(pose::Pose, residue::Residue)
    T = eltype(pose.state)
    chis = Dict{DihedralType, Tuple{T, T}}()
    for chi_index in 1:(length(Dihedral.chi_dict[residue.name.content]) - 1)
        chi = getfield(ProtoSyn.Peptides.Dihedral, Symbol("chi$chi_index"))
        #            Value                                  Standard Deviation
        chis[chi] = (getdihedral(pose.state, chi(residue)), T(0))
    end
    return Rotamer{T}(residue.name, chis)
end