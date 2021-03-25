"""
    ResidueName(s::String) <: AbstractString

A string overload to accomodate different known and expected denominations of
certain residues ("HIS" == "HIE" for example).

# Example
```jldoctest
julia> residue.name = ResidueName("HIS")
```
"""
mutable struct ResidueName <: AbstractString
    content::String
end

Base.show(io::IO, rn::ResidueName) = begin
    print(rn.content)
end

Base.:(==)(s1::String, s2::ResidueName)      = s2 == s1
Base.:(==)(s1::ResidueName, s2::ResidueName) = s1 == s2.content
Base.:(==)(s1::ResidueName, s2::String)      = begin
    if s1.content in ["HIE", "HIS"] && s2 in ["HIE", "HIS"]
        return true
    else
        return s1.content == s2
    end
end

Base.iterate(rn::ResidueName) = Base.iterate(rn.content)
Base.iterate(rn::ResidueName, n::Int64) = Base.iterate(rn.content, n)