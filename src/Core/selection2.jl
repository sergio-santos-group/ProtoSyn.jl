abstract type AbstractSelection end

struct Selection <: AbstractSelection
    body::Function
    is_exit_node::Bool
end

# function (sele::Selection)(container::AbstractContainer)
#     if sele.is_exit_node
#         mask = sele.body(container)
#         # Apply mask to container and return actual selection
#         return mask # CHANGE
#     else
#         # Return mask only
#         return sele.body(container)
#     end
# end

# function (sele::Selection)(pose::Pose)
#     return sele.body(pose.graph)
# end


# Should be in types.jl
# Returns an array with all elements of DataType loopup in the AbstractContainer
# even if the given abstract container contains another container who in turns
# contains the desired DataType
# Ex. retrieve(pose.graph, ProtoSyn.AbstractAtom)
#   > returns all Atom instances on the Pose.
function retrieve(container::AbstractContainer, lookup::DataType)
    results = [container]
    child = typeof(container).super

    while child.super.parameters[1] != Nothing
        if child == lookup
            return results
        else
            child = child.super.parameters[1]
            new_results = []
            for item in results
                for _item in item.items
                    push!(new_results, _item)
                end
            end
            results = new_results
        end
    end
    if child == lookup
        return results
    end
    return Nothing
end

# --- RES NAME SELECTOR

export @resname
macro resname(residue_name::String)

    function return_resname(container::AbstractContainer)

        # 1. Verify that the provided container can contain Residues
        residues = retrieve(container, ProtoSyn.AbstractResidue)
        if residues == Nothing
            return Nothing
        end

        # 2. Loop over the residues and gather the desired residue names
        results = Vector{Residue}()
        for residue in residues
            if residue.name == residue_name
                push!(results, residue)
            end
        end
        return results
    end

    return Selection(return_resname)
end


macro resname(expr::Expr)
    # Macro with multiple selections

    function return_resname(container::AbstractContainer)

        # 1. Verify that the provided container can contain Residues
        residues = retrieve(container, ProtoSyn.AbstractResidue)
        if residues == Nothing
            return Nothing
        end

        # 2. Loop over the residues and gather the desired residue names
        results = Vector{Residue}()
        for residue in residues
            if residue.name == expr.args[2]
                push!(results, residue)
            end
        end
        return results
    end

    extra_selections = eval(expr.args[3])
    return eval(expr.args[1])(Selection(return_resname), extra_selections)
end


# --- AND SELECTOR

# struct AndSelection <: AbstractSelection

#     left::AbstractSelection
#     right::AbstractSelection
#     is_exit_node::Bool
#     body::Function

#     AndSelection(left::AbstractSelection, right::AbstractSelection, is_exit_node::Bool) = new(left, right, is_exit_node, generate_body(left, right, intersect))

# end

function generate_body(left::AbstractSelection, right::AbstractSelection, f::Function)
    return function(container)
        return f(left(container), right(container))
    end
end

# --- OR SELECTOR
struct BinarySelection <: AbstractSelection

    left::AbstractSelection
    right::AbstractSelection
    is_exit_node::Bool
    body::Function

    # AndSelection(left::AbstractSelection, right::AbstractSelection, is_exit_node::Bool) = new(left, right, is_exit_node, generate_body(left, right, union))
end

# For some reason this isn't working !!!
# Base.show(io::IO, ::MIME"text/plain", ::Type{T}) where {T <: AbstractSelection} = begin
#     println("Instance of AbstractSelection")
# end

# function (sele::AndSelection)(container::AbstractContainer)
#     if sele.is_exit_node
#         # Return actual selection
#         println("RETURN ACTUAL SELECTION")
#     else
#         # Return mask
#         return intersect(sele.left(container), sele.right(container))
#     end
# end

Base.:&(a::AbstractSelection, b::AbstractSelection, is_exit_node::Bool) = begin
    return BinarySelection(a, b, is_exit_node, generate_body(a, b, intersect))
end

function(sele::AbstractSelection)(container::AbstractContainer)
    # Decides whether to return a mask or an actual selection based on sele.is_exit_node
    println("Calculating selection ...")
    if sele.is_exit_node
        # Retun an actual selection after calculating the mask
        sele.body(container)
        return "Returning actual selection HERE !"
    else
        # Return a mask
        return sele.body(container)
    end
end

# function (sele::OrSelection)(container::AbstractContainer)
#     # TODO
#     # This function must deal with problems when Mask{Residue} meets Mask{Atom}

#     if sele.is_exit_node
#         println("RETURN ACTUAL SELECTION")
#         # Return actual selection
#     else
#         # Return mask
#         return union(sele.left(container), sele.right(container))
#     end
# end

Base.:|(a::AbstractSelection, b::AbstractSelection, is_exit_node::Bool) = begin
    return BinarySelection(a, b, is_exit_node, generate_body(a, b, union))
end

# ---
# NOTE: If we highjack the intersect and union functions we can return AbstractContainers (?)

export @teste
macro teste(expr::Expr, level::Int = 0)
    
    # TODO This function must be defined elsewhere, somehow
    # TODO This function must return a mask
    _teste = generate_resname_selector(expr.args[2])
    println("EXPRESSION: $(expr.args) ($(typeof(expr)))")

    push!(expr.args[3].args, level + 1)
    extra_selections = eval(expr.args[3])

    if level == 0
        return eval(expr.args[1])(Selection(_teste, false), extra_selections, true)
    else
        return eval(expr.args[1])(Selection(_teste, false), extra_selections, false)
    end
end


function generate_resname_selector(residue_name::String)
    return function(container::AbstractContainer)
        println("Searching container for all instances of residue $residue_name")
        return [1, 0, 0]
    end
end


macro teste(str::String, level::Int = 0)

    # TODO This function must be defined elsewhere, somehow
    # TODO This function must return a mask
    _teste = generate_resname_selector(str)

    if level == 0
        return Selection(_teste, true)
    else
        return Selection(_teste, false)
    end
end

# Generate a function with a blank space
# export generate_selector

# function generate_selector(selector::Expr, x::String)
#     selector.args[2].args[2] = x
#     return selector
# end

# const selector1 = quote
#     residue_name = to_fill
#     println(residue_name)
# end

export @within
macro within(distance::Int, expr::Expr, level::Int = 0)

    # TODO This function must be defined elsewhere, somehow
    # TODO This function must return a mask
    
    # Level control: The current level is carried to the next call of this same
    # function (or similar) by pushing its value into the next expression
    # arguments 
    push!(expr.args, level + 1)
    # println("EXPRESSION: $(expr) ($(typeof(expr)))")
    # println("EXPRESSION: $(expr.args[1])")
    extra_selections = eval(expr)

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
    end
end