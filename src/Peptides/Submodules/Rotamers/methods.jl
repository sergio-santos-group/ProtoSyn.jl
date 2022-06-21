"""
    apply!(state::State, rotamer::Rotamer, residue::Residue)

Apply the given [`Rotamer`](@ref) `rotamer` to the sidechain of
[`Residue`](@ref) `residue`, on the provided [`State`](@ref) `state`. Return the
applied [`Rotamer`](@ref) `rotamer`, without applying the
[`sync!`](@ref ProtoSyn.sync!) method to the [`State`](@ref) `state`.

# Examples
```
julia> ProtoSyn.Peptides.apply!(pose.state, rot_lib["GLU"][35°, -35°][1], pose.graph[1][2])
Rotamer{Float64}: GLU | Chi1:   -55.6° ±  0.0 | Chi2:   187.7° ±  0.0 | Chi3:    -2.6° ±  0.0 | Chi4: --   
```
"""
function apply!(state::State, rotamer::Rotamer, residue::Residue)
    @assert rotamer.name == residue.name "Tried to apply a $(rotamer.name) rotamer to residue $(residue.name)."
    T = eltype(state)

    chis = Dict{AbstractSelection, Tuple{Opt{T}, T}}()
    for (chi, (value, sd)) in rotamer.chis
        value === nothing && begin
            chis[chi] = (value, T(0))
            continue
        end
        _value = (randn() * sd) + value
        if chi(residue, gather = true)[1] === nothing
            @warn "Tried to apply $chi on residue $residue, but the requested atom was not found. Not performing any action."
            chis[chi] = (NaN, T(0))
        else
            setdihedral!(state, chi(residue, gather = true)[1], _value)
            chis[chi] = (_value, T(0))
        end
    end

    return Rotamer(rotamer.name, chis)
end

apply!(state::State, rotamer::Nothing, residue::Residue) = begin
    @warn "Tried to apply a rotamer to $(residue.name), but no Rotamer was provided (rotamer = nothing)"
    return nothing
end


"""
    sample_rotamer(rs::BBI_RotamerLibrary{T}, [n::Int = -1]) where {T <: AbstractFloat}

Sample a [`Rotamer`](@ref) instance from the given [`BBI_RotamerLibrary`](@ref)
`rs`, taking the natural probability of occurrence into consideraction. If a `n`
value is given, sample only from the `n` most likely [`Rotamer`](@ref)
instances. If `n` is 0 or lower (-1, by default), sample from all
[`Rotamer`](@ref) instances. Return the sampled [`Rotamer`](@ref) instance.

# Examples
```
julia> ProtoSyn.Peptides.sample_rotamer(rot_lib["GLU"][35°, -35°])
Rotamer{Float64}: GLU | Chi1:   -66.8° ±  8.2 | Chi2:   179.1° ± 11.8 | Chi3:   -35.0° ±  8.5 | Chi4: --         
```
"""
function sample_rotamer(rs::BBI_RotamerLibrary{T}, n::Int = -1) where {T <: AbstractFloat}
    if n <= 0
        return StatsBase.sample(rs.rotamers, rs.weights)
    else
        n = min(n, length(rs.rotamers))
        return StatsBase.sample(rs.rotamers[1:n], rs.weights[1:n])
    end
end

Base.push!(rs::BBI_RotamerLibrary{T}, rotamer::Rotamer{T}, weight::T) where {T <: AbstractFloat} = begin

    push!(rs.rotamers, rotamer)
    rs.weights = Weights(vcat(rs.weights, Weights([weight])))
end


"""
    sample_rotamer(rl::Dict{String, BBI_RotamerLibrary}, residue::Residue, n::Int = -1)

Sample a [`Rotamer`](@ref) instance from the corresponding
[`BBI_RotamerLibrary`](@ref) for the given [`Residue`](@ref) `residue` name.
This method takes the natural probability of occurrence into consideraction. If
a `n` value is given, sample only from the `n` most likely [`Rotamer`](@ref)
instances. If `n` is 0 or lower (-1, by default), sample from all
[`Rotamer`](@ref) instances. Return the sampled [`Rotamer`](@ref) instance.

# Examples
```
julia> ProtoSyn.Peptides.sample_rotamer(rot_lib, pose.graph[1][24])
Rotamer{Float64}: SER | Chi1:   179.7° ±  8.4 | Chi2: --              | Chi3: --              | Chi4: --      
```
"""
function sample_rotamer(rl::Dict{String, BBI_RotamerLibrary}, residue::Residue, n::Int = -1)
    
    !(residue.name in keys(rl)) && return nothing
    rotamer_stack = rl[residue.name]

    return Peptides.sample_rotamer(rotamer_stack, n)
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

Base.getindex(rs::BBI_RotamerLibrary, i::Int) = begin
    return rs.rotamers[i]
end

Base.getindex(rl::BBD_RotamerLibrary, phi::T, psi::T) where {T <: AbstractFloat} = begin
    phi_index = findnearest(rl.phis, phi)
    psi_index = findnearest(rl.psis, psi)
    return rl.rotamer_stacks[phi_index, psi_index]
end


"""
    sample_rotamer(pose::Pose, rl::Dict{String, BBD_RotamerLibrary}, residue::Residue, [n::Int = -1], [random_inexistent_phi_psi::Bool = false])

Sample a [`Rotamer`](@ref) instance from the corresponding
[`BBD_RotamerLibrary`](@ref) for the given [`Residue`](@ref) `residue` name.
Since the [`BBD_RotamerLibrary`](@ref) is backbone-dependent, the current
[State](@ref state-types) of the provided [`Pose`](@ref) `pose` is read to sample a
[`Rotamer`](@ref). This method takes the natural probability of occurrence into
consideraction. If a `n` value is given, sample only from the `n` most likely
[`Rotamer`](@ref) instances. If `n` is 0 or lower (-1, by default), sample from
all [`Rotamer`](@ref) instances. If `random_inexistent_phi_psi` is set to `true`
(`false`, by default), randomly assigns a phi or psi angle value when none
exists (for example, in terminal [`Residue`](@ref) instances), allowing for a
[`Rotamer`](@ref) to still be sampled. Return the sampled [`Rotamer`](@ref)
instance.

# Examples
```
julia> ProtoSyn.Peptides.sample_rotamer(pose, rot_lib, pose.graph[1][6])
Rotamer{Float64}: GLU | Chi1:   -65.5° ±  9.5 | Chi2:    81.7° ±  9.6 | Chi3:     7.3° ±  8.5 | Chi4: --       
```
"""
function sample_rotamer(pose::Pose, rl::Dict{String, BBD_RotamerLibrary}, residue::Residue, n::Int = -1, random_inexistent_phi_psi::Bool = false)

    phi_atom = Peptides.phi(residue)
    if phi_atom === nothing
        if !random_inexistent_phi_psi
            @warn "Tried to sample a rotamer for residue $residue but the required atom for phi determination was not found. Consider setting the `random_inexistent_phi_psi` flag to true."
            return nothing
        else
            phi = rand(-π:π)
        end
    else
        phi = getdihedral(pose.state, phi_atom)
    end

    psi_atom = Peptides.psi(residue)
    if psi_atom === nothing
        if !random_inexistent_phi_psi
            @warn "Tried to sample a rotamer for residue $residue but the required atom for psi determination was not found. Consider setting the `random_inexistent_phi_psi` flag to true."
            return nothing
        else
            psi = rand(-π:π)
        end
    else
        psi = getdihedral(pose.state, psi_atom)
    end

    if !(residue.name in keys(rl))
        @warn "Tried to sample a rotamer for residue $residue but no entries were found in the provided rotamer library for this residue type."
        return nothing
    end
    rotamer_stack = rl[residue.name][phi, psi]

    return Peptides.sample_rotamer(rotamer_stack, n)
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
    chis = Dict{AbstractSelection, Tuple{Union{T, Nothing}, T}}()
    for chi_index in 1:(length(chi_dict[residue.name.content]) - 1)
        chi = ChiSelection(chi_index)
        c = chi(residue, gather = true)[1]
        if c === nothing
            chis[chi] = (c, T(0))
        else
            chis[chi] = (getdihedral(pose.state, c), T(0))
        end
    end
    return Rotamer{T}(residue.name, chis)
end