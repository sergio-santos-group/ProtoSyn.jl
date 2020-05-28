# using ...ProtoSyn
Atom=Int
Segment=Int
include("types.jl")

#const BondedList{N} where N = Dict{Union{Int,Symbol}, Vector{NTuple{N, Int}}}
export genbonded, loadparm

const BondedList = Dict{Int, Vector{Tuple{Vararg{Atom}}}}

function dfs(f::Function, pivot::Atom, path::Vector{Atom}, maxdepth::Int, depth::Int=0)
    # path[depth+1] = pivot
    push!(path, pivot)
    pivot.visited = true
    f(depth, path)
    if depth < maxdepth
        for atom in pivot.bonds
            !atom.visited && dfs(f, atom, path, maxdepth, depth+1)
        end
    end
    pivot.visited = false
    pop!(path)
end


function genbonded(graph::Segment, depth::Int)
    
    blist = BondedList()
    byatom = eachatom(graph)
    path = Vector{ProtoSyn.Atom}()
    # path = Vector{ProtoSyn.Atom}(undef, depth+1)

    for atom in byatom
        atom.visited = false
    end

    for i = 1:depth
        blist[i] = Vector{NTuple{i,Atom}}()
    end
    for atom in byatom
        dfs(atom, path, depth) do d,pth
            if pth[1].index < pth[end].index
               push!(blist[d], Tuple(pth))
            end
        end
    end
    blist
end

# struct ForcefieldSpec
#     name::String
#     genpairs::Bool
#     fudgeLJ::Float64
#     fudgeQQ::Float64
#     exclusion_depth::Int
#     components::Dict{Symbol, Dict{String, AbstractFFComponentType}} = Dict()
# end

@inline function genkey(types::T...) where {T<:AbstractString}
    #join(types, ":")
    types
    # tuple(types...)
end

struct ForcefieldParameters
    exclusion_depth::Int
    components::Dict{Symbol, Dict{KeyType, AbstractComponentType}}
    ForcefieldParameters() = new(3, Dict())
end

using YAML
function loadparm(::Type{T}, filename::AbstractString) where {T<:AbstractFloat}
    @info "Loading forcefield parameters from file '$filename'"
    yml = YAML.load(open(filename))
    
    parm = ForcefieldParameters()

    for (name,components) in yml["components"]
        container = get!(parm.components, Symbol(name), Dict())
        println(name, "------------------")

        for (typename, items) in components
            println(">> ", typename)
            type = getfield(Forcefield, Symbol(typename)){String,T}
            println("type: ", type)
            println("fieldcount: ", fieldcount(type))
            println("keylen: ", keylen(type))
            # nfields = parmcount(type)
            nfields = fieldcount(type)-keylen(type)
            println("nfields: ", nfields)
            println("fieldtypes: ", fieldtypes(type))
            ftypes = fieldtypes(type)[2:end]
            println(ftypes)
            klen = keylen(type)
            for item in items
                # println(item[end-nfields:end])
                #key = 
                #container[key] = type(item...)
                args = item[klen+1:end]
                args = [
                    convert(t, v isa Vector ? Tuple(v) : v)
                    # convert(t, v)
                    for (t,v) in zip(ftypes,args)
                ]
                key = genkey(item[1:klen]...)
               
                # println(type("", args...))
                println(key, ", ", args)
                println(type(key, args...))
                # push!(container, T(key, args...))
                container[key] = type(key, args...)
            end
        end

    end
    parm
end
loadparm(filename::AbstractString) = loadparm(Float64, filename)

Base.getindex(ffparms::ForcefieldParameters, comp::Symbol, keys::String...) = begin
    
    components = ffparms.components[comp]
    
    haskey(components, keys) && return components[keys]

    if comp === :improper
        key = ("X", keys[2:end]...)
        haskey(components, key) && return components[key]
        
        key = ("X", "X", keys[3], keys[4])
        haskey(components, key) && return components[key]
    elseif comp === :proper
        key = reverse(keys)
        haskey(components, key) && return components[key]
    
        key = ("X", keys[2], keys[3], "X")
        haskey(components, key) && return components[key]

        key = ("X", keys[3], keys[2], "X")
        haskey(components, key) && return components[key]
    else
        key = reverse(keys)
        haskey(components, key) && return components[key]
    end
    nothing
end


function genff(graph, ffparams::ForcefieldParameters, resmaps::Dict)
    blist = genbonded(graph, ffparams.exclusion_depth)
    #excluded = genexcluded(graph, ffparams.exclusion_depth)

    for residue in eachresidue(graph)
        resmap = get(resmaps, residue.name, nothing)
        resmap===nothing && error("unable to find map for residue $(residue.name)")

        # iterate over additional ff component names/values
        for (compname, addcomps) in resmap
            blistitem = get!(blist, Symbol(compname)) do; []; end
            for addcomp in addcomps
                push!(blistitem, cprod(residue, addcomp...))
            end
        end

    end

end

"""
!!!
    NOTE: requires attention
"""
# function cprod(residue::Residue, names::String...)
#     stack = []
#     for name in names
#         if startswith(name, '-') && hasparent(residue) && haskey(residue.parent, name)
#             push!(stack, residue.parent[name])
#         elseif startswith(name, '+') && haschildren(residue) && haskey(residue.children[1], name)
#             push!(stack, residue.children[1][name])
#         elseif haskey(residue, name)
#             push!(stack, residue[name])
#         else
#             return
#         end
#     end
#     stack
# end


#const HarmonicBond{T}     where {T<:AbstractFloat} = HarmonicThing{2,Int,T}

# const HarmonicBond{T}  where {T<:AbstractFloat} =  Bond{HarmonicType{T}}

# for subtype in subtypes(XX)
#     @eval begin
#     end
# end
