
abstract type AbstractPotential end
const KeyType = Tuple{Vararg{T,N}} where {T,N}

keylen(::Type{T}) where {T<:AbstractPotential} = length(T.types[1].types)

genkey(t::T...) where {T<:AbstractString} = join(t,":")


#------------------------------------------------
using YAML

export ForcefieldParameters

struct ForcefieldParameters{T<:AbstractFloat}
    name::String
    fudgeLJ::T
    fudgeQQ::T
    genpairs::Bool
    exclusion_depth::Int
    components::Dict{Symbol, Dict{String, AbstractPotential}}
    # ForcefieldParameters{T}() where {T} = new{T}(3, Dict())
end

ForcefieldParameters(::Type{T}, filename::AbstractString) where {T<:AbstractFloat} = begin
    @info "Loading forcefield parameters from file '$filename'"
    yml = YAML.load_file(filename)
    
    parm = ForcefieldParameters{T}(
        yml["name"],
        yml["fudgeLJ"],
        yml["fudgeQQ"],
        yml["genpairs"],
        yml["exclusion_depth"],
        Dict()
    )

    for (name,components) in yml["components"]
        container = get!(parm.components, Symbol(name), Dict())
        
        for (typename, items) in components
            # get the adequate parametric type and make it a
            # concrete type with string keys and fields of type T
            if !isdefined(Forcefield, Symbol(typename))
                error("Unknown component type '$typename'")
            end
            
            type = getfield(Forcefield, Symbol(typename)){String,T}
            
            ftypes = fieldtypes(type)[2:end]
            klen = keylen(type)
            for item in items
                args0 = item[klen+1:end]
                args1 = [
                    convert(t, v isa Vector ? Tuple(v) : v)
                    for (t,v) in zip(ftypes,args0)
                ]
                key = genkey(item[1:klen]...)
                container[key] = type(Tuple(item[1:klen]), args1...)
                # println(key, type(key, args1...))
            end
        end
    end
    parm

end
ForcefieldParameters(filename::AbstractString) = ForcefieldParameters(Float64, filename)




Base.eltype(::ForcefieldParameters{T}) where {T} = T

Base.getindex(ffparms::ForcefieldParameters, comp::Symbol, keys::String...) = begin
    # THIS SHOULD GO TO THE AMBER FOLDER
    
    components = ffparms.components[comp]
    
    key = genkey(keys...)
    haskey(components, key) && return components[key]

    if comp === :improper
        key = genkey("X", keys[2:end]...)
        haskey(components, key) && return components[key]
        
        key = genkey("X", "X", keys[3], keys[4])
        haskey(components, key) && return components[key]
    elseif comp === :proper
        key = genkey(reverse(keys)...)
        haskey(components, key) && return components[key]
    
        key = genkey("X", keys[2], keys[3], "X")
        haskey(components, key) && return components[key]

        key = genkey("X", keys[3], keys[2], "X")
        haskey(components, key) && return components[key]
    else
        key = genkey(reverse(keys)...)
        haskey(components, key) && return components[key]
    end
    nothing
end







#------------------------------------------------

const ExclusionList = Dict{Int, Vector{Int}}

struct Forcefield2
    exclusions::ExclusionList
    components::Dict{DataType, Vector{<:AbstractPotential}}
end
Forcefield2() = Forcefield2(Dict(), Dict())

Base.push!(ff::Forcefield2, item::T) where {T<:AbstractPotential} = begin
    container = get!(ff.components, T) do; Vector{T}(); end
    push!(container, item)
    ff
end


const BondedList = Dict{Union{Int,Symbol}, Vector{Tuple{Vararg{Atom}}}}
