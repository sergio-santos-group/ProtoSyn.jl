using YAML

using ..ProtoSyn
using ..ProtoSyn.Units: tonumber

"""
    build(grammar::LGrammar{T}, derivation)

Build a new [`Pose`](@ref) instance using the given `derivation` sequence on the
provided [`LGrammar`](@ref) `grammar` instructions. Return the generated
[`Pose`](@ref) after synching (using the [`sync!`](@ref) method).

# See Also
[`fragment`](@ref)

# Examples
```
julia> res_lib = ProtoSyn.Peptides.grammar;

julia> pose = ProtoSyn.build(res_lib, seq"GME")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 39
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function build(grammar::LGrammar{T}, derivation) where {T<:AbstractFloat}
    top = Topology("UNK", 1)
    state = State(T)
    state.id = top.id
    pose = Pose(top, state)

    if !isempty(derivation)
        frag = fragment(grammar, derivation)
        # Appending the fragment (which is a segment) to the Topology
        append_fragment_as_new_segment!(pose, frag)
        
        ProtoSyn.request_i2c!(state; all=true)
    end

    sync!(pose)
    pose
end


export @seq_str

"""
    @seq_str(s::String)

Construct a vector of strings from the provided string. Helpful when providing a
_derivation_ to any building method (such as [`build`](@ref)).

# Short syntax
* seq"..."

# Examples
```jldoctest
julia> seq"ABC"
3-element Vector{String}:
 "A"
 "B"
 "C"
```
"""
macro seq_str(s::String); [string(c) for c in s]; end


export append_fragment!
"""
    append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op::Any = "α")

Based on the provided `grammar`, add the residue sequence from `derivation` to
the given [`Pose`](@ref) `pose`, appending it after the given [`Residue`](@ref)
`residue`. This residue and the new [`Fragment`](@ref) will be connected using
operation `op` ("α" by default). Request internal to cartesian coordinate
conversion and return the altered [`Pose`](@ref) `pose`.

# See also
[`insert_fragment!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.append_fragment!(pose, pose.graph[1][36], res_lib, seq"MMM")
Pose{Topology}(Topology{/2a3d:532}, State{Float64}:
 Size: 628
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)```
"""
function append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op::Any = "α")

    # Build the fragment to append
    frag = ProtoSyn.fragment(grammar, derivation)

    return ProtoSyn.append_fragment!(pose, residue, grammar, frag, op = op)
end


export insert_fragment!
"""
    insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op::Any = "α", connect_upstream::Bool = true)

Based on the provided `grammar`, add the residue sequence from `derivation` to
the given `pose`, inserting it on the position of the given [`Residue`](@ref)
instance `residue` (the `residue` gets shifted downstream). The first downstream
[`Residue`](@ref) and the new [`Fragment`](@ref) will be connected using
operation `op` ("α" by default). If `connect_upstream` is set to true (is, by
default), also connect to the upstream [`Residue`](@ref) instances using the
same operation `op`. Request internal to cartesian coordinate conversion and
return the altered [`Pose`](@ref) `pose`.

# See also
[`append_fragment!`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.unbond!(pose, pose.graph[1][1]["C"], pose.graph[1, 2, "N"])
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 343
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)

julia> ProtoSyn.insert_fragment!(pose, pose.graph[1][2], res_lib, seq"A")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 353
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; op::Any = "α", connect_upstream::Bool = true)

    # Build the fragment to insert
    frag = ProtoSyn.fragment(grammar, derivation)

    return insert_fragment!(pose, residue, grammar, frag, op = op)
end