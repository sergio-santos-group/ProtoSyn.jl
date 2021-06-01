"""
    build(grammar::LGrammar{T}, derivation, ss::NTuple{3,Number} = SecondaryStructure[:linear]) where {T <: AbstractFloat}

Build a [`Pose`](@ref) using the given `derivation` sequence on the provided
`grammar` instructions. If an `ss` is provided, automatically apply it to the
built pose (linear secondary structure, by default).

!!! ukw "Note:"
    This function is an extension of [`ProtoSyn.build`](@ref).

# See also
[`setss!`](@ref)

# Examples
```
julia> pose = ProtoSyn.Peptides.build(res_lib, seq"QQQ")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 51
 i2c: false | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function build(grammar::LGrammar{T}, derivation, ss::SecondaryStructureTemplate = SecondaryStructure[:linear]) where {T <: AbstractFloat}

    pose = ProtoSyn.build(grammar, derivation)
    ProtoSyn.Peptides.setss!(pose, ss)
    sync!(pose)
    pose
end


"""
    append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; [ss::Opt{SecondaryStructureTemplate} = nothing], [op = "α"])

Add the a [`Fragment`](@ref) built from the given `derivation` to the provided
[`Pose`](@ref) `pose`, appending it after the given [`Residue`](@ref) `residue`.
This residue and the new [`Fragment`](@ref) `frag` will be connected using
operation `op` ("α" by default) of the given [`LGrammar`](@ref) `grammar`. If
given, a [`SecondaryStructureTemplate`](@ref) `ss` can be applied to the new
appendage (using the [`setss!`](@ref) method). In either case, the `C=O` bond
position is re-calculated and set (in the anchor for the first residue of the
appendage). Request internal to cartesian coordinate conversion and return the
altered [`Pose`](@ref) `pose`.

!!! ukw "Note:"
    This function is an extension of [`ProtoSyn.append_fragment!`](@ref).

# See also
[`insert_fragment!`](@ref ProtoSyn.Peptides.insert_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Opt{SecondaryStructureTemplate}, ::Any))

# Examples
```jldoctest
julia> ProtoSyn.Peptides.append_fragment!(pose, pose.graph[1][end], res_lib, seq"AAA")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 373
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function append_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")
    frag = ProtoSyn.fragment(grammar, derivation)
    return ProtoSyn.Peptides.append_fragment!(pose, residue, grammar, frag, ss = ss,
        op = op)
end


"""
    insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")

Insert the [`Fragment`](@ref) built from the given `derivation` in the provided
`pose`, on the position of the provided [`Residue`](@ref) instance `residue`
(the `residue` gets shifted downstream). This first downstream [`Residue`](@ref)
and the new [`Fragment`](@ref) will be connected using operation `op` ("α" by
default) from [`LGrammar`] `grammar`. Also connects to the upstream
[`Residue`](@ref) instance, using the same operation. If given, a
[`SecondaryStructureTemplate`](@ref) `ss` can be applied to the new appendage
(using the [`setss!`](@ref) method). If the appendage is not being inserted at
the [`root`](@ref), the `C=O` bond position is re-calculated and set (in the
anchor for the first residue of the appendage). If the appendage is being
inserted at the [`root`](@ref), perform a soft uncap of the terminal hydrogen
atoms (removes "H2" and "H3", leaves "H1", renames it to "H") and recalculate
the N-H bond position (at the first downstream [`Residue`](@ref)). Request
internal to cartesian coordinate conversion and return the altered
[`Pose`](@ref) `pose`.

!!! ukw "Note:"
    This function is an extension of [`ProtoSyn.insert_fragment!`](@ref).

# See also
[`append_fragment!`](@ref ProtoSyn.append_fragment!(::Pose{Topology}, ::Residue, ::LGrammar, ::Pose{Segment}; ::Opt{SecondaryStructureTemplate}, ::Any))

# Examples
```jldoctest
julia> ProtoSyn.Peptides.insert_fragment!(pose, pose.graph[1][1], res_lib, seq"AAA")
Pose{Topology}(Topology{/UNK:1}, State{Float64}:
 Size: 373
 i2c: true | c2i: false
 Energy: Dict(:Total => Inf)
)
```
"""
function insert_fragment!(pose::Pose{Topology}, residue::Residue, grammar::LGrammar, derivation; ss::Opt{SecondaryStructureTemplate} = nothing, op = "α")
    frag = ProtoSyn.fragment(grammar, derivation)
    return ProtoSyn.Peptides.insert_fragment!(pose, residue, grammar, frag, ss = ss,
        op = op)
end