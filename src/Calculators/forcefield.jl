using YAML

export ForcefieldParameters

struct ForcefieldParameters{T<:AbstractFloat}
    name::String
    fudgeLJ::T
    fudgeQQ::T
    genpairs::Bool
    exclusion_depth::Int
    components::Dict{Symbol, Dict{String, AbstractPotential}}
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
            #if !isdefined(Forcefield, Symbol(typename))
            #    error("Unknown component type '$typename'")
            #end
            
            # type = getfield(Forcefield, Symbol(typename)){String,T}
            type = typefor(Symbol(typename)){String,T}
            
            klen = keylen(type)
            for item in items
                key = genkey(item[1:klen]...)
                args = map(tonumber, item[klen+1:end])
                if haskey(container, key)
                    container[key] += type(Tuple(item[1:klen]), args...)
                else
                    container[key] = type(Tuple(item[1:klen]), args...)
                end
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




struct Forcefield
    exclusions::ExclusionList
    components::Dict{Symbol, Vector{<:AbstractPotential}}
end
Forcefield() = Forcefield(Dict(), Dict())


Base.push!(ff::Forcefield, item::T) where {T<:AbstractPotential} = begin
    container = get!(ff.components, name(T)) do; Vector{T}(); end
    push!(container, item)
    ff
end


