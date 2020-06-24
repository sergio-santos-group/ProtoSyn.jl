abstract type AbstractSelection end
abstract type AbstractStateMode end

# StateModes are used as parameters in the Selection structs and not as a
# parametric type. Parametric functors seems not possible in Julia as of yet:
# https://stackoverflow.com/questions/46063872/parametric-functors-in-julia
# This can be overcome by NOT using a functor, and instead calling a 'select'
# function. Ex. select(::Selection, ::Pose). This second option would require
# the existence of a resolving function for each combination of StateMode and
# Selection type, such as, for example:
# select(sele::Selection{Stateless}, ...);
# select(sele::Selection{Statefull}, ...);
# select(sele::BinarySelection{Stateless}, ...);
# select(sele::BinarySelection{Statefull}); etc
# Having them as parameters in the structs means a single resolving function
# needs to exist, who then does the necessary partition between Statefull and
# Stateless calls. No difference exists between Selection and BinarySelection
# calls, as they both resolve to calling either the existing saved mask, or
# calculating a new mask with sele.body(container).
struct Stateless <: AbstractStateMode end
struct Statefull <: AbstractStateMode end

export @resname, @segname, @atomname, @atomsymb, @atomid, @atomix
export @res, @seg, @atom

# Look-up menu
# Returns the correct function to loop over the given AbstractContainer type
const lkm = Dict(Segment => eachsegment, Residue => eachresidue, Atom => eachatom)

# --- MASK ---------------------------------------------------------------------
# A Mask struct is necessary, instead of a simple BitVector, to be able to know
# the type being enumerated, for promotion/demotion tasks
struct Mask{T <: AbstractContainer}
    content::BitVector

end
Mask() = Mask{AbstractContainer}(BitVector())
Mask{T}() where {T <: AbstractContainer} = Mask{T}(BitVector())

@inline Base.length(m::Mask) = length(m.content)

# --- SELECTION ----------------------------------------------------------------

mutable struct Selection <: AbstractSelection
    body::Function
    is_exit_node::Bool
    state_mode::Type{<:AbstractStateMode}
    mask::Mask
end
Selection(body::Function, is_exit_node::Bool, state_mode::Type{<:AbstractStateMode}) = Selection(body, is_exit_node, state_mode, Mask())

# Title function parses the selection function description to return a readable
# string with the identification of that function
function title(sele::Selection)::String
    s = string(sele)
    ss = vcat([split(x, "(") for x in split(s, ")")]...)
    if length(split(ss[7], ",")) <= 1
        return split(ss[5], ",")[1]
    else
        return *(split(ss[7], ",")[2:-1:1]...)[2:end]
    end
end

# Name function return a decription of the whole Selector for printing
function name(sele::Selection, level::Int = 0)::String
    prefix = repeat("  ", level)
    str  = prefix*"- Selection\n"
    str *= prefix*"Exit node : $(sele.is_exit_node)\n"
    str *= prefix*"State mode: $(sele.state_mode)\n"
    str *= prefix*"Name      : $(title(sele))"
    return str
end

# --- BINARY SELECTION ---------------------------------------------------------

mutable struct BinarySelection <: AbstractSelection

    left::AbstractSelection
    right::AbstractSelection
    is_exit_node::Bool
    body::Function
    state_mode::Type{<:AbstractStateMode}
    mask::Mask
end
BinarySelection(left::AbstractSelection, right::AbstractSelection,
    is_exit_node::Bool, body::Function, state_mode::Type{<:AbstractStateMode}) =
    BinarySelection(left, right, is_exit_node, body, state_mode, Mask())


function title(sele::BinarySelection)::String
    s = string(sele)
    m = [m.captures[1] for m in eachmatch(r"ProtoSyn\.(\w{2,3})\)", s)][end]
    return uppercasefirst(m)
end


function name(sele::BinarySelection, level::Int = 0)::String
    prefix = repeat("  ", level)
    str  = prefix*"- BinarySelection\n"
    str *= prefix*"Exit node : $(sele.is_exit_node)\n"
    str *= prefix*"State mode: $(sele.state_mode)\n"
    str *= prefix*"Name      : $(title(sele))\n"
    str *= name(sele.left, level + 1)*"\n"
    str *= name(sele.right, level + 1)
    return str
end


Base.show(io::IO, ::MIME"text/plain", s::T) where {T <: AbstractSelection} =
    print(io, name(s))


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

# ADD TITLE

function or(left::Mask{T}, right::Mask{T})::Mask{T} where {T <: AbstractContainer}
    return Mask{T}(left.content .| right.content)
end

function or(left::Bool, right::Bool)::Bool
    return left || right
end

function and(left::Mask{T}, right::Mask{T})::Mask{T} where {T <: AbstractContainer}
    return Mask{T}(left.content .& right.content)
end

function and(left::Bool, right::Bool)::Bool
    return left && right
end

# --- PROMOTION / DEMOTION OF MASKS --------------------------------------------

# type_rule always returns two values, in the following order:
# 1) the lowest ranking type of the two types given
# 2) the highest ranking type of the two types given
function type_rule(m1::Mask{T1}, m2::Mask{T2})::Tuple{Mask, Mask} where {T1, T2 <: AbstractContainer}
    levels = Dict(Atom => 1, Residue => 2, Segment => 3, Topology => 4)
    
    m1_type = typeof(m1).parameters[1]
    m2_type = typeof(m2).parameters[1]
    return levels[m1_type] < levels[m2_type] ? (m1, m2) : (m2, m1)
end


function demote(left_mask::Mask, right_mask::Mask, container::AbstractContainer, f::Function)::Mask
    # Promotion / Demotion of Mask Type
    # 1. Find lowest and highest ranking types in the conflict
    low_mask, high_mask = type_rule(left_mask, right_mask)
    low_type            = typeof(low_mask).parameters[1]
    high_type           = typeof(high_mask).parameters[1]
    
    if low_type == high_type
        return f(high_mask, low_mask)
    end 

    # 2. Create empy mask.content with same size of lowest ranking type
    mask = falses(length(low_mask.content))
    index = 1

    # 3. Iterate over both masks and perform the associated function
    # (and / or) on the correct items, whose result is saved in the new
    # Mask, to return.
    for (high_index, high_item) in enumerate(lkm[high_type](container))
        for low_item in lkm[low_type](high_item)
            mask[index] = f(high_mask.content[high_index], low_mask.content[index])
            index += 1
        end
    end

    return Mask{low_type}(mask)
end


# This function is responsible for generating the body of BinarySelections. This
# function is, in turn, responsible for performing the promotion/demotion of
# Masks based on their types, if they are different from each other. Right now,
# this function only produces demotion results.
function generate_body(left::AbstractSelection, right::AbstractSelection, f::Function, state_mode::Type{T}) where {T <: AbstractStateMode}
    if state_mode == Stateless
        return function(container::AbstractContainer, force_update::Bool = false)
            left_mask       = left(container, force_update)
            right_mask      = right(container, force_update)
            
            if typeof(left_mask).parameters[1] == typeof(right_mask).parameters[1]
                return f(left_mask, right_mask)
            else
                return demote(left_mask, right_mask, container, f)
            end
        end
    else
        return function(container::AbstractContainer, state::State, force_update::Bool = false)
            if left.state_mode == Stateless
                left_mask   = left(container, force_update)
            else
                left_mask   = left(container, state, force_update)
            end

            if right.state_mode == Stateless
                right_mask   = right(container, force_update)
            else
                right_mask   = right(container, state, force_update)
            end
            
            if typeof(left_mask).parameters[1] == typeof(right_mask).parameters[1]
                return f(left_mask, right_mask)
            else
                return demote(left_mask, right_mask, container, f)
            end
        end
    end
end

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
    body       = generate_body(a, b, and, state_mode)
    return BinarySelection(a, b, is_exit_node, body, state_mode)
end

Base.:|(a::AbstractSelection, b::AbstractSelection, is_exit_node::Bool) = begin
    state_mode = state_rule(a.state_mode, b.state_mode)
    body       = generate_body(a, b, or, state_mode)
    return BinarySelection(a, b, is_exit_node, body, state_mode)
end

# This functions are employed in cases where the user tries to group together
# 2 or more selections using parenthesis.
# Ex. (@resname "ALA" | @resname "GLN") & @atomname "CA"
# The function joins both selections, setting the is_exit_node of the resulting
# BinarySelection to true while setting both "children" selections is_exit_node
# to false.
Base.:&(a::AbstractSelection, b::AbstractSelection) = begin
    state_mode     = state_rule(a.state_mode, b.state_mode)
    body           = generate_body(a, b, and, state_mode)
    a.is_exit_node = false
    b.is_exit_node = false
    return BinarySelection(a, b, true, body, state_mode)
end


Base.:|(a::AbstractSelection, b::AbstractSelection) = begin
    state_mode     = state_rule(a.state_mode, b.state_mode)
    body           = generate_body(a, b, and, state_mode)
    a.is_exit_node = false
    b.is_exit_node = false
    return BinarySelection(a, b, true, body, state_mode)
end

# -- RESOLVE FUNCTION ----------------------------------------------------------

# This is the Selection resolve function. It returns either an actual selection
# of AbstractContainers or a mask, depending on whether it is an exit node or
# not, respectively. 
function(sele::AbstractSelection)(container::AbstractContainer, force_update::Bool = false)

    @assert sele.state_mode == ProtoSyn.Stateless "Statefull Selections require a State. Please use Selection(::AbstractContainer, ::State)"
    
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

function(sele::AbstractSelection)(container::AbstractContainer, state::State, force_update::Bool = false)
   
    # 1. If there is already a calculated mask and the force_update flag was not
    # set to true, return the saved mask.
    if sele.state_mode == Statefull || length(sele.mask.content) == 0 || force_update
        # 1.1 The function used to actually calculate the selection mask depends
        # on the state_mode of the selection. Stateless selections do not
        # require a State.
        if sele.state_mode == Statefull
            sele.mask = sele.body(container, state, force_update)
        else
            sele.mask = sele.body(container, force_update)
        end
    end

    # 2. Return a mask or an actual selection based on sele.is_exit_node
    if sele.is_exit_node
        # Change this to return actual selection!
        return sele.mask
    else
        return sele.mask
    end
end

# --- NOT SELECTION ------------------------------------------------------------

# Necessary functions to be able to select the negative of current selection
# Ex: (!@atom)(pose.graph)
# Ex: (!@resname "ALA")(pose.graph)

function Base.:!(sele::Selection)
    old_body = deepcopy(sele.body)
    sele.body = function(container::AbstractContainer, force_update::Bool = false)
        return !((old_body)(container, force_update))
    end

    sele.mask = Mask{typeof(sele.mask).parameters[1]}(.!sele.mask.content)
    return sele
end

function Base.:!(m::Mask{T}) where {T <: AbstractContainer}
    return Mask{T}(.!(m.content))
end


# --- GENERATE DEFAULT SELECTORS -----------------------------------------------
# --- 1. PROPERTY SELECTOR GENERATOR -------------------------------------------

# Note all generators must return a function who, in turn, returns a Mask{T}
# This generator should return one such function to get the given 'field' equal
# to 'name', for the defined t type
const AcceptableTypes = Union{String, Int, Vector{Any}, Vector{String}, Regex, Symbol}
function generate_property_selector_body(name::AcceptableTypes, t::DataType, field::Symbol)::Function

    # lkm = Dict(Segment => eachsegment, Residue => eachresidue, Atom => eachatom)
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

@generate_property_macro( segname, Segment,   :name)
@generate_property_macro( resname, Residue,   :name)
@generate_property_macro(atomname,    Atom,   :name)
@generate_property_macro(atomsymb,    Atom, :symbol)
@generate_property_macro(  atomid,    Atom,     :id)
@generate_property_macro(  atomix,    Atom,  :index)

# --- 2. TRUE SELECTIONS -------------------------------------------------------

function generate_true_selector_body(t::DataType)::Function

    lookup_function     = lkm[t]

    return function(container::AbstractContainer, force_update::Bool = false)

        m = trues(lookup_function(container).size[end])
        return Mask{t}(m)
    end
end

macro generate_true_macro(fname, t)
    quote
        macro $(esc(fname))(expr::Expr, level::Int...)
            # Composite entry (Ex: @atom & @resname "GLN")

            level = length(level) == 0 ? (0,) : level # Set default level to 0

            # Increase the level when calling the next selector, so that it
            # is not marked as an exit_node
            push!(expr.args[1].args, level[1] + 1)

            # Evaluate the next selector. Should return a Selector to merge 
            # with the current one, using the exprs.args[1] call.
            extra_selections = eval(expr.args[1])

            rs = generate_true_selector_body(eval($(esc(t))))

            if level[1] == 0
                return eval(expr.head)(
                    Selection(rs, false, Stateless),
                    extra_selections,
                    true)
            else
                return eval(expr.head)(
                    Selection(rs, false, Stateless),
                    extra_selections,
                    false)
            end
        end

        macro $(esc(fname))(level::Int...)
            # Simple parameter entry (Ex: @atom)

            level = length(level) == 0 ? (0,) : level # Set default level to 0
            rs = generate_true_selector_body($(esc(t)))

            if level[1] == 0
                return Selection(rs, true, Stateless)
            else
                return Selection(rs, false, Stateless)
            end
        end
    end
end

@generate_true_macro( seg, Segment)
@generate_true_macro( res, Residue)
@generate_true_macro(atom,    Atom)

# --- WITHIN SELECTOR ----------------------------------------------------------

export @within
macro within(distance::Real, expr::Expr, level::Int...)

    # Level control: The current level is carried to the next call of this same
    # function (or similar) by pushing its value into the next expression
    # arguments 
    level = length(level) == 0 ? (0,) : level # Set default level to 0
    push!(expr.args, level[1] + 1)

    extra_selections = eval(expr)
    
    rs = generate_within_selector(distance, extra_selections, or)

    if level[1] == 0
        return Selection(rs, true, Statefull)
    else
        return Selection(rs, false, Statefull)
    end
end

# Here X and Y should be of ANY type because we don't know what are we giving to
# the function (?)
function generate_within_selector(distance::Real, selection::AbstractSelection, rule::Function)
    
    co_sq = distance*distance

    return function(container::AbstractContainer, state::State, force_update::Bool = false)
        mask = selection(container)
        # 1) Demote whatever type of Mask 'selection' returns to Mask{Atom} 
        atom_mask = demote(mask, (@atom)(container), container, and)

        
        # 2) Loop over the selected Atoms and return neighbours
        # Note: Here, the distances between atoms are calculated in a naive
        # approach. In future versions, however, a global distance matrix may be
        # employed.
        masks = Array{Mask{Atom}, 1}()
        for (atom_i_index, atom_selected) in enumerate(atom_mask.content)
            if atom_selected
                _mask = Mask{Atom}(zeros(eachatom(container).size[end]))
                atom_i = state[atom_i_index].t
                # State.items includes the pseudoatoms from 'origin', so the
                # iterationmust start on index 4
                j_atoms = state.items[(3 + atom_i_index):end]
                for (atom_j_index, atom_j) in enumerate(state.items[4:end])
                    d_sq = sum(@. (atom_i - atom_j.t)^2)
                    if d_sq <= co_sq
                        _mask.content[atom_j_index] = true
                    end
                end
                push!(masks, _mask)
            end
        end

        # 3) Reduce the array of masks to a single mask based on a rule
        return reduce(rule, masks)
    end
end

# --- PRINT SELECTION TO FILE (TEMPORARY) --------------------------------------

export print_selection

# This temporary function is used as a hot fix to show selections in a PDB file
# format, where selected atoms are in red and non-selected atoms are displayed
# in blue, when loaded into a visualizer such as PyMOL.
function print_selection(io::IOStream, pose::Pose{Topology}, mask::Mask{T}) where {T <: AbstractContainer}
    
    atom_mask = demote(mask, (@atom)(pose.graph), pose.graph, and)

    Base.write(io, "MODEL\n")
    for (atom_index, atom) in enumerate(eachatom(pose.graph))
        sti = pose.state[atom.index] # returns and AtomState instance

        # In this file, selected atoms will be displayed in red while
        # non-selected atoms will be displayed in blue
        atom_symbol = atom_mask.content[atom_index] ? "O" : "N"

        s = @sprintf("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%24s\n",
            atom.index, atom_symbol,
            atom.container.name, atom.container.container.code,
            atom.container.id,
            sti.t[1], sti.t[2], sti.t[3],
            atom_symbol)
            Base.write(io, s)
    end

    for atom in eachatom(pose.graph)
        Base.write(io, @sprintf("CONECT%5d", atom.index))
       foreach(n->Base.write(io, @sprintf("%5d",n.index)), atom.bonds)
       Base.write(io,"\n")
    end
    Base.write(io, "ENDMDL")
end


# --- TODO: --------------------------------------------------------------------
# - @dihedral selector
# - @within selector DONE
# - @withinall selector DONE
# - Mask{T} struct DONE.
# - Statefull and Stateless functions. Some of these functions must receive not
# only an AbstractContainer but also a State with coordinates DONE
# - Deal with situations where Stateless and Statefull selectors mix DONE
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
# happen. Residues become Atoms, per general rule. DONE
# - Give variables to selectors. Ex. a = [1, 2 ,3]; (@atomid a)(pose.graph)
# - True selections (Ex. @atom) DONE
# - Get demoted results on demand. Ex. @atom & @resname "ALA" as list of all
# atoms whose residue they belong to is an "ALA". DONE
# - Get promoted results on demand. Ex. @res & @atomname "CA" as list of all
# residues who contain at least one "CA" atom. Does this make sense?
# - Return actual AbstractContainer instances when resolving an exit_node
# selection

# pose = sync!(ProtoSyn.Builder.build(ProtoSyn.Peptides.grammar(), ProtoSyn.Builder.seq"AAQG"))
# io = open("../teste.pdb", "w"); print_selection(io, pose, (@resname "ALA")(pose.graph)); close(io)