function promote(m1::Mask{T1}, m2::Mask{T2}, container::AbstractContainer) where {T1, T2}
    # Always lower type

    type1 = selection_type(m1)
    type2 = selection_type(m2)
    if type1 > type2
        m1 = promote(m1, type2, container)
    elseif type2 > type1
        m2 = promote(m2, type1, container)
    end

    m1, m2
end

function promote(mask::Mask{T1}, ::Type{T2}, container::AbstractContainer, f::Function = any)::Mask where {T1 <: AbstractContainer, T2 <: AbstractContainer}

    new_mask = Mask{T2}(counter(T2)(container))

    if T1 > T2
        count = counter(T2)
        iter = iterator(T1)
        i = 1
        for (entry, item) in zip(mask, iter(container))
            n = count(item)
            new_mask[i:(i+n-1)] = entry
            i += n
        end
    else
        count = counter(T1)
        iter = iterator(T2)
        i = 1
        for (index, item) in enumerate(iter(container))
            n = count(item)
            new_mask[index] = f(view(mask.content, i:(i+n-1)))
            i += n
        end
    end

    return new_mask
end