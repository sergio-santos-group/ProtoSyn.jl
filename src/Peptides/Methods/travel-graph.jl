struct IUPACSA <: ProtoSyn.SearchAlgorithm end

function (sa::IUPACSA)(atom::Atom, stack::Vector{Atom})
    bonds = copy(ProtoSyn.sort_children(atom))
    if atom.name == "CA"
        c_index = findfirst((atom) -> atom.name == "C", bonds)
        if c_index === nothing
            @warn "Tried to sort $atom children following IUPAC convention, but no \"C\" atom was found."
        else
            push!(bonds, bonds[c_index])
            deleteat!(bonds, c_index)
        end
    end
    return vcat(bonds, stack)
end

"""
    (ProtoSyn.IUPAC)(atom::Atom, stack::Vector{Atom})

IUPAC-like search algorithm for [`travel_graph`](@ref). Correctly sorts the
given [`Atom`](@ref) `atom` children instances and concatenates with the current
`stack`. This method attempts to identify each [`Residue`](@ref) sidechain (from
alpha-carbon CA children names) and sorts sidechain before the remaining
backbone (C=O & N, etc). Consider using
[`assign_default_atom_names!`](@ref ProtoSyn.Peptides.assign_default_atom_names!)
to make sure all [`Atom`](@ref) names are as expected.

# Examples
```
julia> ProtoSyn.Peptides.IUPAC(pose.graph[1][1]["CA"], Vector{Atom}())
3-element Vector{Atom}:
 Atom{/test:36441/A:1/MET:1/HA:6}
 Atom{/test:36441/A:1/MET:1/CB:8}
 Atom{/test:36441/A:1/MET:1/C:7}
```
"""
IUPAC = IUPACSA()


"""
    travel_graph(start::Atom; [stop::Opt{Atom} = nothing], [search_algorithm::F = Peptides.IUPAC]) where {F <: SearchAlgorithm})

Overload of the [`travel_graph`](@ref ProtoSyn.travel_graph) method, where the
default `search_algorithm` is [`IUPAC`](@ref ProtoSyn.Peptides.IUPAC).

# Examples
```
julia> Peptides.travel_graph(pose.graph[1][1]["CA"])
15-element Vector{Atom}:
 Atom{/test:37905/A:1/MET:1/CA:5}
 Atom{/test:37905/A:1/MET:1/HA:6}
 Atom{/test:37905/A:1/MET:1/CB:8}
 Atom{/test:37905/A:1/MET:1/HB3:10}
 Atom{/test:37905/A:1/MET:1/HB2:11}
 Atom{/test:37905/A:1/MET:1/CG:12}
 Atom{/test:37905/A:1/MET:1/HG3:13}
 Atom{/test:37905/A:1/MET:1/HG2:14}
 Atom{/test:37905/A:1/MET:1/SD:15}
 Atom{/test:37905/A:1/MET:1/CE:16}
 Atom{/test:37905/A:1/MET:1/HE3:17}
 Atom{/test:37905/A:1/MET:1/HE2:18}
 Atom{/test:37905/A:1/MET:1/HE1:19}
 Atom{/test:37905/A:1/MET:1/C:7}
 Atom{/test:37905/A:1/MET:1/O:9}
```
"""
travel_graph(start::Atom; stop::Opt{Atom} = nothing, search_algorithm::F = IUPAC) where {F <: ProtoSyn.SearchAlgorithm} = begin
    ProtoSyn.travel_graph(start, stop = stop, search_algorithm = search_algorithm)
end