
export @sym_str, @sn_str, @rn_str, @an_str, @x_str
export select


struct ById end
struct ByName end
struct BySymbol end

struct Stateless end
struct Statefull end

abstract type AbstractSelection end


state_rule(::Type{Stateless}, ::Type{Stateless}) = Stateless
state_rule(::Type{Statefull}, ::Type{Stateless}) = Statefull
state_rule(::Type{Stateless}, ::Type{Statefull}) = Statefull
state_rule(::Type{Statefull}, ::Type{Statefull}) = Statefull


clear(s::AbstractSelection...) = foreach(x->x.isentry=false, s)

matcher(n::AbstractString) = begin
    startswith(n,'@') ? (Regex(n[2:end]), occursin) : (n, isequal)
end


#region LEAF SELECTORS

# =========================================================
mutable struct TrueSelection <: AbstractSelection
    isentry::Bool
    TrueSelection() = new(true)
end

state_type(::Type{TrueSelection}) = Stateless

name(io::IO, s::TrueSelection, prefix="", suffix="") = begin
    print(io, prefix)
    print(io, "TrueSelection()")
    println(io, suffix)
end



mutable struct TripleSelection <: AbstractSelection
    isentry::Bool
    apattern::Union{AbstractString,Regex}
    amatch::Function
    rpattern::Union{AbstractString,Regex}
    rmatch::Function
    spattern::Union{AbstractString,Regex}
    smatch::Function
    TripleSelection(an::AbstractString, rn::AbstractString, sn::AbstractString) = begin
        ap,am = matcher(an)
        rp,rm = matcher(rn)
        sp,sm = matcher(sn)
        new(true, ap, am, rp, rm, sp, sm)
    end
end

state_type(::Type{TripleSelection}) = Stateless

name(io::IO, s::TripleSelection, prefix="", suffix="") = begin
    print(io, prefix)
    print(io, "TripleSelection($(s.apattern),$(s.rpattern),$(s.spattern))")
    println(io, suffix)
end

# =========================================================
mutable struct Selection{T} <: AbstractSelection
    isentry::Bool
    pattern::Union{AbstractString,Regex}
    match::Function
    Selection{T}(n::AbstractString) where T = begin
        # if startswith(n, '@')
        #     return new{T}(true, Regex(n[2:end]), occursin)
        # end
        # new{T}(true, n, isequal)
        p, m = matcher(n)
        new{T}(true, p, m)
    end
end

state_type(::Type{Selection{T}}) where {T} = Stateless

name(io::IO, s::Selection{T}, prefix="", suffix="") where T = begin
    print(io, prefix)
    print(io, "Selection{$(nameof(T))}($(s.pattern))")
    println(io, suffix)
end

#endregion

#region COMPOUND SELECTORS
# =========================================================
mutable struct BinarySelection{L,R} <: AbstractSelection
    isentry::Bool
    op::Function
    left::AbstractSelection
    right::AbstractSelection
    BinarySelection(op::Function, l::L, r::R) where {L,R} = begin
        clear(l, r)
        new{state_type(L),state_type(R)}(true, op, l, r)
    end
end

state_type(::Type{BinarySelection{L,R}}) where {L,R} = state_rule(L,R)

name(io::IO, s::BinarySelection{L,R}, prefix="",suffix="") where {L,R} = begin
    print(io, prefix)
    println(io, "BinarySelection{$(nameof(L)),$(nameof(R))}($(s.op),")
    name(io, s.left,  prefix*"  ", ",")
    name(io, s.right, prefix*"  ")
    println(io, prefix, ")", suffix)
end

# =========================================================
mutable struct UnarySelection{T} <: AbstractSelection
    isentry::Bool
    op::Function
    sele::AbstractSelection
    element_wise::Bool
    UnarySelection(op::Function, sele::T; element_wise::Bool) where T = begin
        clear(sele)
        new{state_type(T)}(true, op, sele, element_wise)
    end
end

state_type(::Type{UnarySelection{T}}) where {T} = T

name(io::IO, s::UnarySelection{T}, prefix="",suffix="") where T = begin
    print(io, prefix)
    println(io, "UnarySelection{$(nameof(T))}($(s.op), element_wise=$(s.element_wise),")
    name(io, s.sele,  prefix*"  ")
    println(io, prefix, ")", suffix)
end

# =========================================================
mutable struct DistanceSelection{T} <: AbstractSelection
    isentry::Bool
    distance::Number
    sele::AbstractSelection
    DistanceSelection(distance::Number, sele::T) where T = begin
        clear(sele)
        new{state_type(T)}(true, distance, sele)
    end
end

state_type(::Type{DistanceSelection{T}}) where T = Statefull

name(io::IO, s::DistanceSelection{T}, prefix="", suffix="") where T = begin
    print(io, prefix)
    println(io, "DistanceSelection{$(nameof(T))}($(s.distance),")
    name(io, s.sele, prefix*"  ")
    println(io, prefix, ")", suffix)
end

#endregion COMPOUND SELECTORS


Base.show(io::IO, s::AbstractSelection) = name(io, s)


macro sym_str(s); Selection{BySymbol}(s); end
macro  sn_str(s); Selection{Segment}(s); end
macro  rn_str(s); Selection{Residue}(s); end
macro  an_str(s); Selection{Atom}(s); end

# PROBLEM 1 - This doesn't allow selection by other fields, such as 'symbol' or
# others that might be added later.

macro x_str(s)
    fields = reverse(split(s, '/'))
    nfields = length(fields)
    nfields > 3 && error("Invalid selector")
    if all(isempty, fields)
        selector = TrueSelection()
    else
        an = nfields > 0 && !isempty(fields[1]) ? fields[1] : "@"
        rn = nfields > 1 && !isempty(fields[2]) ? fields[2] : "@"
        sn = nfields > 2 && !isempty(fields[3]) ? fields[3] : "@"
        selector = TripleSelection(an, rn, sn)
    end
    selector
end



@inline _collect(ac::AbstractContainer, s::AbstractSelection, m::BitVector) =
    s.isentry ? select(ac, m) : m

@inline select(ac::AbstractContainer, mask::BitVector) =
    collect(at for (m,at) in zip(mask,eachatom(ac)) if m)

select(ac::AbstractContainer, s::TrueSelection) = _collect(ac, s, trues(count_atoms(ac)))
select(ac::AbstractContainer, s::Selection{BySymbol}) = _select(ac, s, :symbol)
select(ac::AbstractContainer, s::Selection{Residue}) = _select(ac, s, :container, :name)
select(ac::AbstractContainer, s::Selection{Atom}) = _select(ac, s, :name)
select(top::Topology, s::Selection{Segment}) = _select(top, s, :container, :container, :name)



#getf(ac, fields::Symbol...) = isempty(fields) ? ac : getf(getproperty(ac, fields[1]), fields[2:end]...)
getf(ac, fields::Symbol...) = begin
    for field in fields;
        ac = getproperty(ac, field)
    end
    ac
end


# QUESTION: Why more than 1 fields... ? Doesn't it just save the last one?
_select(ac::AbstractContainer, s::Selection, fields::Symbol...) = begin
    println("evaluating static")
    mask = falses(count_atoms(ac))
    for (i,atom) in enumerate(eachatom(ac))
        mask[i] = s.match(s.pattern, getf(atom, fields...))
    end
    _collect(ac, s, mask)
end



Base.:&(l::AbstractSelection, r::AbstractSelection) = BinarySelection(&, l, r)
Base.:&(l::AbstractSelection, r::TrueSelection) = l
Base.:&(l::TrueSelection, r::AbstractSelection) = r
Base.:&(l::TrueSelection, r::TrueSelection) = l

Base.:|(l::AbstractSelection, r::AbstractSelection) = BinarySelection(|, l, r)
Base.:|(l::AbstractSelection, r::TrueSelection) = r
Base.:|(l::TrueSelection, r::AbstractSelection) = l
Base.:|(l::TrueSelection, r::TrueSelection) = l

# Base.:(:)(l::Number, r::AbstractSelection) = DistanceSelection(l, r)
# Base.:(:)(l::AbstractSelection, r::Number) = DistanceSelection(r, l)
Base.:(:)(s::AbstractSelection, w::Number, of::AbstractSelection) = BinarySelection(&, s, DistanceSelection(w, of))

Base.:!(l::AbstractSelection) = UnarySelection(!, l; element_wise=true)
Base.any(l::AbstractSelection) = UnarySelection(any, l; element_wise=false)
Base.all(l::AbstractSelection) = UnarySelection(all, l; element_wise=false)


# select(ac::AbstractContainer, s::BinarySelection) = begin
#     lmask = select(ac, s.left)
#     rmask = select(ac, s.right)
#     mask = (s.op).(lmask, rmask)
#     _collect(ac, s, mask)
# end

# # select(ac::AbstractContainer, s::NegateSelection) = begin
# #     mask = .!select(ac, s.sele)
# #     _collect(ac, s, mask)
# # end


# # Base.:&(l::AbstractStatefullSelection, r::AbstractSelection) =
# #    StatefullBinarySelection(&, l, r)
# # Base.:&(l::AbstractSelection, r::AbstractStatefullSelection) =
# #     StatefullBinarySelection(&, l, r)

# # select(ac::AbstractContainer, s::DistanceSelection) = begin
# #     mask = select(ac, s.sele)
# #     return function (state::State)
# #         println(mask)
# #         println(state.id)
# #         println("selecting atom within $(s.distance) of $(s.sele)")
# #         _collect(ac,s,mask)
# #     end
# # end


# -------------------------
select(ac::AbstractContainer, s::TripleSelection) = begin
    mask = falses(count_atoms(ac))
    for (i,atom) in enumerate(eachatom(ac))
        mask[i] = s.amatch(s.apattern, atom.name) && 
                  s.rmatch(s.rpattern, atom.container.name) &&
                  s.smatch(s.spattern, atom.container.container.name)
    end
    _collect(ac, s, mask)
end

# -------------------------

select(ac::AbstractContainer, s::UnarySelection{Stateless}) = begin
    mask = select(ac, s.sele)
    if s.element_wise
        mask = (s.op).(mask)
    else
        mask .= s.op(mask)
    end
    _collect(ac, s, mask)
end

# -------------------------
select(ac::AbstractContainer, s::DistanceSelection{Stateless}) = begin
    mask = select(ac, s.sele)
    return function (state::State)
        println("processing state\n$s")
        # <calc distances> & mask
        _collect(ac, s, mask)
    end
end

select(ac::AbstractContainer, s::DistanceSelection{Statefull})::Function = begin
    selector = select(ac, s.sele)
    println("evaluating dynamic\n$s")
    return function (state::State)
        mask = selector(state)
        println("processing state $s")
        # <calc distances> & mask
        _collect(ac, s, mask)
    end
end

# -------------------------

select(ac::AbstractContainer, s::BinarySelection{Stateless,Stateless}) = begin
    lmask = select(ac, s.left)
    rmask = select(ac, s.right)
    _collect(ac, s, (s.op).(lmask, rmask))
end

select(ac::AbstractContainer, s::BinarySelection{Stateless,Statefull}) = _select(ac, s, s.left, s.right)
select(ac::AbstractContainer, s::BinarySelection{Statefull, Stateless}) = _select(ac, s, s.right, s.left)

_select(ac::AbstractContainer, s::BinarySelection, sless::AbstractSelection, sfull::AbstractSelection) = begin
    lmask = select(ac, sless)
    rselector = select(ac, sfull)
    return function (state::State)
        rmask = rselector(state)
        _collect(ac, s, (s.op).(lmask, rmask))
    end
end

select(ac::AbstractContainer, s::BinarySelection{Statefull,Statefull}) = begin
    lselector = select(ac, s.left)
    rselector = select(ac, s.right)
    return function (state::State)
        lmask = lselector(state)
        rmask = rselector(state)
        
        _collect(ac, s, (s.op).(lmask, rmask))
    end
end

#---------------------
abstract type AbstractSerialSelector end
struct IndexSelector <: AbstractSerialSelector
    val::Int
end

export ix
const ix = IndexSelector(1)

struct IdSelector <: AbstractSerialSelector
    val::Int
end

export id
const id = IdSelector(1)

Base.:(*)(x::Number, y::IndexSelector) = IndexSelector(x)
Base.:(*)(x::Number, y::IdSelector) = IdSelector(x)

function select(ac::AbstractContainer, s::IdSelector)
    for item in ac.items
        if item.id == s.val
            return item
        end
    end
    nothing
end

