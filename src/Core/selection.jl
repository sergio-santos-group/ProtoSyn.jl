
export @sym_str, @sn_str, @rn_str, @an_str, @x_str
export select


struct ById end
struct ByName end
struct BySymbol end

const Static=Val{0}
const Dynamic=Val{1}

abstract type AbstractSelection end

mutable struct TrueSelection <: AbstractSelection
    isentry::Bool
    TrueSelection() = new(true)
end

mutable struct Selection{B} <: AbstractSelection
    isentry::Bool
    pattern::Union{AbstractString,Regex}
    match::Function
    Selection{B}(n::AbstractString) where B = begin
        s = new{B}(true)
        if startswith(n, '@')
            s.pattern = Regex(n[2:end])
            s.match = occursin
        else
            s.pattern = n
            s.match = isequal
        end
        # new{B}(true, startswith(n, '@') ? Regex(n[2:end]) : Regex("^$n\$"))
        s
    end
end

mutable struct NegateSelection <: AbstractSelection
    isentry::Bool
    sele::AbstractSelection
    NegateSelection(sele::AbstractSelection) = begin
        sele.isentry = false
        new(true, sele)
    end
end

mutable struct BinarySelection <: AbstractSelection
    isentry::Bool
    op::Function
    left::AbstractSelection
    right::AbstractSelection
    BinarySelection(op::Function, l::AbstractSelection, r::AbstractSelection) = begin
        l.isentry = false
        r.isentry = false
        new(true, op, l, r)
    end
end

abstract type AbstractDynamicSelection <: AbstractSelection end
mutable struct DistanceSelection <: AbstractDynamicSelection
    isentry::Bool
    distance::Number
    sele::AbstractSelection
    DistanceSelection(distance::Number, sele::AbstractSelection) = begin
        sele.isentry = false
        new(true, distance, sele)
    end
end


name(io::IO, s::NegateSelection, prefix="", suffix="") where T = begin
    println(io, "$(prefix)!(")
    name(io, s.sele, prefix*"  ")
    println(io, "$prefix)$suffix")
end

name(io::IO, s::Selection{T}, prefix="", suffix="") where T = begin
    println(io, "$(prefix)Selection{$(nameof(T))}($(s.pattern))$suffix")
end

name(io::IO, s::TrueSelection, prefix="", suffix="") = begin
    println(io, "$(prefix)TrueSelection()$suffix")
end

name(io::IO, s::BinarySelection, prefix="",suffix="") = begin
    println(io, "$prefix$(s.op)(")
    name(io, s.left,  prefix*"  ", ",")
    name(io, s.right, prefix*"  ")
    println(io, "$prefix)$suffix")
end

Base.show(io::IO, s::AbstractSelection) = name(io, s)


# [sra]*(n|id)

macro sym_str(s); Selection{BySymbol}(s); end
macro  sn_str(s); Selection{Segment}(s); end
macro  rn_str(s); Selection{Residue}(s); end
macro  an_str(s); Selection{Atom}(s); end

macro x_str(s)
    fields = reverse(split(s, '/'))
    nfields = length(fields)
    s = TrueSelection()
    if nfields > 0
        s = isempty(fields[1]) ? s : s & Selection{Atom}(fields[1])
    end
    if nfields > 1
        s = isempty(fields[2]) ? s : s & Selection{Residue}(fields[2])
    end
    if nfields > 2
        s = isempty(fields[3]) ? s : s & Selection{Segment}(fields[3])
    end
    s
end


@inline _collect(ac::AbstractContainer, s::AbstractSelection, m::BitVector) =
    s.isentry ? select(ac, m) : m


@inline select(ac::AbstractContainer, mask::BitVector) =
    (at for (m,at) in zip(mask,eachatom(ac)) if m)

select(ac::AbstractContainer, s::TrueSelection) = _collect(ac, s, trues(count_atoms(ac)))

select(top::Topology, s::Selection{Segment}) = begin
    mask = falses(count_atoms(top))
    for seg in top.items
        if s.match(s.pattern, seg.name)
            for atom in eachatom(seg)
                mask[atom.index] = true
            end
        end
    end
    _collect(top,s,mask)
end

select(ac::AbstractContainer, s::Selection{BySymbol}) = begin
    mask = falses(count_atoms(ac))
    for (i,atom) in enumerate(eachatom(ac))
        mask[i] = s.match(s.pattern, atom.symbol)
    end
    _collect(ac,s,mask)
end

select(ac::AbstractContainer, s::Selection{Residue}) = begin
    mask = falses(count_atoms(ac))
    for (i,atom) in enumerate(eachatom(ac))
        mask[i] = s.match(s.pattern, atom.container.name)
    end
    _collect(ac,s,mask)
end


select(ac::AbstractContainer, s::Selection{Atom}) = begin
    mask = falses(count_atoms(ac))
    for (i,atom) in enumerate(eachatom(ac))
        mask[i] = s.match(s.pattern, atom.name)
    end
    _collect(ac,s,mask)
end

Base.:&(l::AbstractSelection, r::AbstractSelection) = BinarySelection(&, l, r)
Base.:&(l::AbstractSelection, r::TrueSelection) = l
Base.:&(l::TrueSelection, r::AbstractSelection) = r
Base.:&(l::TrueSelection, r::TrueSelection) = l

Base.:|(l::AbstractSelection, r::AbstractSelection) = BinarySelection(|, l, r)
Base.:|(l::AbstractSelection, r::TrueSelection) = r
Base.:|(l::TrueSelection, r::AbstractSelection) = l
Base.:|(l::TrueSelection, r::TrueSelection) = l

Base.:!(l::AbstractSelection) = NegateSelection(l)


select(ac::AbstractContainer, s::BinarySelection) = begin
    lmask = select(ac, s.left)
    rmask = select(ac, s.right)
    mask = (s.op).(lmask, rmask)
    _collect(ac, s, mask)
end

select(ac::AbstractContainer, s::NegateSelection) = begin
    mask = .!select(ac, s.sele)
    _collect(ac, s, mask)
end


# Base.:(:)(l::Number, r::AbstractSelection) = DistanceSelection(l, r)
# Base.:(:)(l::AbstractSelection, r::Number) = DistanceSelection(r, l)
# Base.:&(l::AbstractDynamicSelection, r::AbstractSelection) =
#    DynamicBinarySelection(&, l, r)
# Base.:&(l::AbstractSelection, r::AbstractDynamicSelection) =
#     DynamicBinarySelection(&, l, r)

select(ac::AbstractContainer, s::DistanceSelection) = begin
    mask = select(ac, s.sele)
    return function (state::State)
        println(mask)
        println(state.id)
        println("selecting atom within $(s.distance) of $(s.sele)")
        _collect(ac,s,mask)
    end
end


DistanceSelection{Static}
DistanceSelection{Dynamic}
select(ac::AbstractContainer, s::DistanceSelection{Static}) = begin
    mask = select(ac, s.sele)
    return function (state::State)
        # <calc distances> & mask
        _collect(ac, s, mask)
    end
end
select(ac::AbstractContainer, s::DistanceSelection{Dynamic}) = begin
    selector = select(ac, s.sele)
    return function (state::State)
        mask = selector(state)
        # <calc distances> & mask
        _collect(ac, s, mask)
    end
end
select(ac::AbstractContainer, s::BinarySelection{Static,Static}) = begin
    lmask = select(ac, s.left)
    rmask = select(ac, s.right)
    mask = (s.op).(lmask, rmask)
    _collect(ac, s, mask)
end

select(ac::AbstractContainer, s::BinarySelection{Static,Dynamic}) = begin
    lmask = select(ac, s.left)
    rselector = select(ac, s.left)
    return function (state::State)
        rmask = rselector(state)
        # <calc distances> & mask
        mask = (s.op).(lmask, rmask)
        _collect(ac, s, mask)
    end
end

select(ac::AbstractContainer, s::BinarySelection{Dynamic,Dynamic}) = begin
    lselector = select(ac, s.left)
    rselector = select(ac, s.left)
    return function (state::State)
        lmask = lselector(state)
        rmask = rselector(state)
        # <calc distances> & mask
        mask = (s.op).(lmask, rmask)
        _collect(ac, s, mask)
    end
end