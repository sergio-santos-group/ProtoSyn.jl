abstract type AbstractSelection end
abstract type AbstractStateMode end

# StateModes are used as parameters in the Selection structs and not as a
# parametric type, as this second option would require the existence of a
# resolving function for each combination of StateMode and Selection type:
# Ex: f(sele::Selection{Stateless}), f(sele::Selection{Statefull}), 
# f(sele::BinarySelection{Stateless}), f(sele::BinarySelection{Statefull}), etc
# Having them as parameters in the structs means a single resolving function
# needs to exist, who then does the necessary partition between Statefull and
# Stateless calls. No difference exists between Selection and BinarySelection
# calls, as they both resolve to calling either the existing saved mask, or
# calculating a new mask with sele.body(container).
struct Stateless <: AbstractStateMode end
struct Statefull <: AbstractStateMode end

export @resname, @segname, @atomname, @atomsymb, @atomid, @atomix

# -- MASKS
# A Mask struct is necessary, instead of a simple BitVector, to be able to know
# the type being enumerated, for promotion/demotion tasks
struct Mask{T <: AbstractContainer}
    content::BitVector

end
Mask() = Mask{AbstractContainer}(BitVector())


mutable struct Selection <: AbstractSelection
    body::Function
    is_exit_node::Bool
    state_mode::Type{<:AbstractStateMode}
    mask::Mask
end
Selection(body::Function, is_exit_node::Bool, state_mode::Type{<:AbstractStateMode}) = Selection(body, is_exit_node, state_mode, Mask())


mutable struct BinarySelection <: AbstractSelection

    left::AbstractSelection
    right::AbstractSelection
    is_exit_node::Bool
    body::Function
    state_mode::Type{<:AbstractStateMode}
    mask::Mask
end
BinarySelection(left::AbstractSelection, right::AbstractSelection, is_exit_node::Bool, body::Function, state_mode::Type{<:AbstractStateMode}) =
    BinarySelection(left, right, is_exit_node, body, state_mode, Mask())


# function (sele::Selection)(pose::Pose)
#     return sele(pose.graph)
# end


# Should be in types.jl
# Returns an array with all elements of DataType loopup in the AbstractContainer
# even if the given abstract container contains another container who in turns
# contains the desired DataType
# Ex. retrieve(pose.graph, ProtoSyn.AbstractAtom)
#   > returns all Atom instances on the Pose.
# function retrieve(container::AbstractContainer, lookup::DataType)
#     results = [container]
#     child = typeof(container).super

#     while child.super.parameters[1] != Nothing
#         if child == lookup
#             return results
#         else
#             child = child.super.parameters[1]
#             new_results = []
#             for item in results
#                 for _item in item.items
#                     push!(new_results, _item)
#                 end
#             end
#             results = new_results
#         end
#     end
#     if child == lookup
#         return results
#     end
#     return Nothing
# end

function or(left::Mask{T}, right::Mask{T})::Mask{T} where {T <: AbstractContainer}
    return Mask{T}(left.content .| right.content)
end

function and(left::Mask{T}, right::Mask{T})::Mask{T} where {T <: AbstractContainer}
    return Mask{T}(left.content .& right.content)
end

function generate_body(left::AbstractSelection, right::AbstractSelection, f::Function)
    return function(container, force_update::Bool = false)
        return f(left(container, force_update), right(container, force_update))
    end
end


# For some reason this isn't working !!! FIX THIS
# Base.show(io::IO, ::MIME"text/plain", ::Type{T}) where {T <: AbstractSelection} = begin
#     println("Instance of AbstractSelection")
# end

# --- AND/OR MERGING FUNCTIONS FOR BINARY SELECTORS ----------------------------

# These two functions take two selections and merge them in a BinarySelection,
# which then can be saved, applied, etc. When merging two Selection objects,
# the state_mode is set to Statefull if at least one of the selections occurs
# with that state_mode. This means that, when resolving, if at least one of the
# Selections needs to be updated, all of them are calculated. However, when
# calculating a Stateless selection (belonging to this BinarySelection), if
# force_update is not set to true, the call will always return the saved mask,
# and not calculate it again.

state_rule(::Type{Stateless}, ::Type{Stateless}) = Stateless
state_rule(::Type{Statefull}, ::Type{Stateless}) = Statefull
state_rule(::Type{Stateless}, ::Type{Statefull}) = Statefull
state_rule(::Type{Statefull}, ::Type{Statefull}) = Statefull

Base.:&(a::AbstractSelection, b::AbstractSelection, is_exit_node::Bool) = begin
    state_mode = state_rule(a.state_mode, b.state_mode)
    body       = generate_body(a, b, and)
    return BinarySelection(a, b, is_exit_node, body, state_mode)
end

Base.:|(a::AbstractSelection, b::AbstractSelection, is_exit_node::Bool) = begin
    state_mode = state_rule(a.state_mode, b.state_mode)
    body       = generate_body(a, b, or)
    return BinarySelection(a, b, is_exit_node, body, state_mode)
end

# -- RESOLVE FUNCTION ----------------------------------------------------------

# This is the Selection resolve function. It returns either an actual selection
# of AbstractContainers or a mask, depending on whether it is an exit node or
# not, respectively. 
# Note: This (or other) function should deal with promotion/demotion !
function(sele::AbstractSelection)(container::AbstractContainer, force_update::Bool = false)
    
    # 1. If there is already a calculated mask and the force_update flag was not
    # set to true, return the saved mask.
    if sele.state_mode == Statefull || length(sele.mask.content) == 0 || force_update
        sele.mask = sele.body(container, force_update)
    end

    # 2. Return a mask or an actual selection based on sele.is_exit_node
    if sele.is_exit_node
        # Change this to return actual selection!
        return sele.mask
    else
        return sele.mask
    end
end

# --- PROPERTY SELECTOR BODY GENERATOR -----------------------------------------

# Note all generators must return a function who, in turn, returns a Mask{T}
# This generator should return one such function to get the given 'field' equal
# to 'name', for the defined t type
const AcceptableTypes = Union{String, Int, Vector{Any}, Vector{String}, Regex}
function generate_property_selector_body(name::AcceptableTypes, t::DataType, field::Symbol)::Function

    lkm = Dict(Segment => eachsegment, Residue => eachresidue, Atom => eachatom)
    com = Dict(
        String         => isequal,
        Int            => isequal,
        Vector{Any}    => in,
        Vector{String} => in,
        Regex          => (s::String, r::Regex) -> occursin(r, s))

    comparison_function = com[typeof(name)]
    lookup_function     = lkm[t]

    # This function doesn't require force_update for anything. However, the
    # selection resolving function doesn't know whether it's calling a Selection
    # or a BinarySelection (in which case the force_update forces the included
    # selections to update regardless of StateMode or Mask status). Therefore,
    # this force_update exists to prevent errors.
    return function(container::AbstractContainer, force_update::Bool = false)

        m = falses(lookup_function(container).size[end])
        for (index, item) in enumerate(lookup_function(container))
            if comparison_function(getproperty(item, field), name)
                m[index] = true
            end
        end
        return Mask{t}(m)
    end
end

# --- PROPERTY SELECTOR MACRO GENERATOR ----------------------------------------

# This macro generates other macros, given the 'fname'. Such macros will search
# the final container for a given 'field' (of type 't'), and return a Mask with
# the results. The actual function called to calculate this Mask is generated by
# 'generate_property_selector_body', above.
macro generate_property_macro(fname, t, field)
    quote
        macro $(esc(fname))(expr::Expr, level::Int...)
            L = length(unique([typeof(i) for i in expr.args]))

            level = length(level) == 0 ? (0,) : level # Set default level to 0

            if typeof(expr.args[end]) == Expr
                # Composite entry (Ex: @resname "ALA" | @resname "GLN")
                
                # Increase the level when calling the next selector, so that it
                # is not marked as an exit_node
                push!(expr.args[3].args, level[1] + 1)

                # Evaluate the next selector. Should return a Selector to merge 
                # with the current one, using the exprs.args[1] call.
                extra_selections = eval(expr.args[3])

                # We need to evaluate expr.args[2] because it could be a Regex
                # entry, in which case the generate_property_selector_body
                # requires a Regex object, not the current Expression.
                rs = generate_property_selector_body(
                    eval(expr.args[2]), $(esc(t)), $(esc(field)))

                if level[1] == 0
                    return eval(expr.args[1])(
                        Selection(rs, false, Stateless),
                        extra_selections,
                        true)
                else
                    return eval(expr.args[1])(
                        Selection(rs, false, Stateless),
                        extra_selections,
                        false)
                end
            else
                if L == 1
                    # Multiple parameter entry (Ex: @resnamein ["GLN", "ALA"])
                    rs = generate_property_selector_body(
                        expr.args, $(esc(t)), $(esc(field)))
                else
                    # Regex parameter entry (Ex: @resname r"AL.")
                    rs = generate_property_selector_body(
                        eval(expr), $(esc(t)), $(esc(field)))
                end

                if level[1] == 0
                    return Selection(rs, true, Stateless)
                else
                    return Selection(rs, false, Stateless)
                end
            end
        end

        macro $(esc(fname))(expr::Union{String, Int}, level::Int...)
            # Simple parameter entry (Ex: @resname "ALA")

            level = length(level) == 0 ? (0,) : level # Set default level to 0
            rs = generate_property_selector_body(expr, $(esc(t)), $(esc(field)))

            if level[1] == 0
                return Selection(rs, true, Stateless)
            else
                return Selection(rs, false, Stateless)
            end
        end
    end
end

@generate_property_macro( segname, Segment, :name)
@generate_property_macro( resname, Residue, :name)
@generate_property_macro(atomname,    Atom, :name)
@generate_property_macro(atomsymb,    Atom, :symbol)
@generate_property_macro(  atomid,    Atom, :id)
@generate_property_macro(  atomix,    Atom, :index)

# ---

export @within
macro within(distance::Int, expr::Expr, level::Int = 0)

    # Level control: The current level is carried to the next call of this same
    # function (or similar) by pushing its value into the next expression
    # arguments 
    push!(expr.args, level + 1)
    # println("EXPRESSION: $(expr) ($(typeof(expr)))")
    # println("EXPRESSION: $(expr.args[1])")
    extra_selections = eval(expr)
    
    # TODO This function must be defined elsewhere, somehow
    # TODO This function must return a mask
    _teste = generate_double_selector(distance, extra_selections)

    if level == 0
        return Selection(_teste, true)
    else
        return Selection(_teste, false)
    end
end

# Here X and Y should be of ANY type because we don't know what are we giving to
# the function
function generate_double_selector(x, y)
    return function(container::AbstractContainer)
        println("Looking for residues within $x of $(y(container))")
        return [0, 1, 1]
    end
end

# TODO:
# - @dihedral selector
# - @within selector
# - Mask{T} struct DONE.
# - Statefull and Stateless functions. Some of these functions must receive not
# only an AbstractContainer but also a State with coordinates
# - Some Statefull functions don't need to select everything all the time, and
# should be able to pre-select the stateless part of such function as a seperate
# selection. This only applies when not doing design. This means that after
# selecting something, a stateless Selector object should keep information until
# asked to re-evaluate itself. DONE
# - Select names of residues/atoms with Regular Expressions DONE
# - Selectors for all type properties (name, symbol, etc) This includes Int
# properties, such as :id or :index DONE
# - @resnamein ["GLN", "GLY"] -> Property selectors with vectors DONE
# - When mixing Residue and Atom selections, a promotion/demotion step should
# happend. Residues become Atoms, per general rule.
# - Deal with situations where Stateless and Statefull selectors mix
# - Give variables to selectors. Ex. a = [1, 2 ,3]; (@atomid a)(pose.graph)
# - Get demoted results on demand. Ex. @resname "ALA" as list of all atoms whose
# residue they belong to is an "ALA".

# pose = ProtoSyn.Builder.build(ProtoSyn.Peptides.grammar(), ProtoSyn.Builder.seq"AAQG")