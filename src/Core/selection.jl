
export @rname_str, @name_str, select


const OR = Val{1}
const AND = Val{2}

struct ByAtom end
struct ByResidue end

abstract type AbstractSelection end
struct Selection{B} <: AbstractSelection
    name::AbstractString
end

struct CompositeSelection{B} <: AbstractSelection
    left::Union{Nothing, AbstractSelection}
    right::Union{Nothing, AbstractSelection}
end

# [sr]*(name|id)

macro rname_str(s)
    return Selection{ByResidue}(s)
end
macro name_str(s)
    return Selection{ByAtom}(s)
end

select(t::AbstractContainer, s::Selection{ByResidue}) = begin
    vcat(collect(r.items for r in eachresidue(t) if r.name==s.name)...)
end

select(t::AbstractContainer, s::Selection{ByAtom}) = begin
    collect(r for r in eachatom(t) if r.name==s.name)
end

Base.:&(l::AbstractSelection, r::AbstractSelection) = CompositeSelection{AND}(l, r)
Base.:|(l::AbstractSelection, r::AbstractSelection) = CompositeSelection{OR}(l, r)

select(t::AbstractContainer, s::CompositeSelection{AND}) = begin
    lsele = select(t, s.left)
    rsele = select(t, s.right)
    return intersect(lsele, rsele)
end
select(t::AbstractContainer, s::CompositeSelection{OR}) = begin
    lsele = select(t, s.left)
    rsele = select(t, s.right)
    return union(lsele, rsele)
end

# (s::AbstractSelection)(t::AbstractContainer) = select(t, s)
#(s::T)(t::AbstractContainer) where {T <: AbstractSelection} = select(t, s)