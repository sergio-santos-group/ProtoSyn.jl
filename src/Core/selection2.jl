abstract type AbstractSelection end

struct Selection <: AbstractSelection
    body::Function
end

function (sele::Selection)(container::AbstractContainer)
    return sele.body(container)
end

function (sele::Selection)(pose::Pose)
    return sele.body(pose.graph)
end


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

struct AndSelection <: AbstractSelection

    left::AbstractSelection
    right::AbstractSelection

end

function (sele::AndSelection)(container::AbstractContainer)
    return intersect(sele.left(container), sele.right(container))
end

Base.:&(a::AbstractSelection, b::AbstractSelection) = begin
    return AndSelection(a, b)
end

# --- OR SELECTOR

struct OrSelection <: AbstractSelection

    left::AbstractSelection
    right::AbstractSelection

end

function (sele::OrSelection)(container::AbstractContainer)
    return union(sele.left(container), sele.right(container))
end

Base.:|(a::AbstractSelection, b::AbstractSelection) = begin
    return OrSelection(a, b)
end

# NOTE: If we highjack the intersect and union functions we can return AbstractContainers.