"""
    apply!(state::State, rotamer::Rotamer, residue::Residue)

Apply the given [`Rotamer`](@ref) `rotamer` to the sidechain of
[`Residue`](@ref) `residue`, on the provided [`State`](@ref) `state`. Return the
applied [`Rotamer`](@ref) `rotamer`.

# Examples
```
julia> ProtoSyn.Peptides.apply!(pose.state, rot_lib["GLU"][35°, -35°][1], pose.graph[1][2])
Rotamer{Float64}: GLU | Chi1:   -55.6° ±  0.0 | Chi2:   187.7° ±  0.0 | Chi3:    -2.6° ±  0.0 | Chi4: --   
```
"""
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


"""
    sample(rs::RotamerStack{T}, [n::Int = -1]) where {T <: AbstractFloat}

Sample a [`Rotamer`](@ref) instance from the given [`RotamerStack`](@ref) `rs`,
taking the natural probability of occurrence into consideraction. If a `n` value
is given, sample only from the `n` most likely [`Rotamer`](@ref) instances. If
`n` is 0 or lower (-1, by default), sample from all [`Rotamer`](@ref) instances.
Return the sampled [`Rotamer`](@ref) instance.

# Examples
```
julia> ProtoSyn.Peptides.sample(rot_lib["GLU"][35°, -35°])
Rotamer{Float64}: GLU | Chi1:   -66.8° ±  8.2 | Chi2:   179.1° ± 11.8 | Chi3:   -35.0° ±  8.5 | Chi4: --         
```
"""
function sample(rs::RotamerStack{T}, n::Int = -1) where {T <: AbstractFloat}
    if n <= 0
        return StatsBase.sample(rs.rotamers, rs.weights)
    else
        n = min(n, length(rs.rotamers))
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

Base.getindex(rs::RotamerStack, i::Int) = begin
    return rs.rotamers[i]
end

Base.getindex(rl::R, phi::T, psi::T) where {R <: RotamerLibrary, T <: AbstractFloat} = begin
    phi_index = findnearest(rl.phis, phi)
    psi_index = findnearest(rl.psis, psi)
    return rl.rotamer_stacks[phi_index, psi_index]
end

"""
    get_rotamer(pose::Pose, residue::Residue)

Measure the existent chi dihedral angles in [`Residue`](@ref) `residue` on the
given [`Pose`](@ref) `pose`, saving the information in a new [`Rotamer`](@ref)
instance.

# Examples
```jldoctest
julia> ProtoSyn.Peptides.get_rotamer(pose, pose.graph[1][2])
Rotamer{Float64}: GLU | Chi1:   180.0° ±  0.0 | Chi2:  -180.0° ±  0.0 | Chi3:    90.0° ±  0.0 | Chi4: --     
```
"""
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