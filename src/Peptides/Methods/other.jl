using Printf

"""
    measure_similarity(sequence1::String, sequence2::String)

Showcase the total similarity and average similarity of two peptidic sequences
according to the `ProtoSyn.Peptides.aminoacid_similarity` map (this is a
mutation tolerance measure map by Stephenson et al. (2013) (see
[https://link.springer.com/article/10.1007/s00239-013-9565-0](https://link.springer.com/article/10.1007/s00239-013-9565-0))).

# Examples
```
julia> ProtoSyn.Peptides.measure_similarity(ProtoSyn.sequence(pose1), ProtoSyn.sequence(pose2))
 (...)
```
"""
function measure_similarity(sequence1::String, sequence2::String)
    if length(sequence1) !== length(sequence2)
        @warn "The provided sequences have different lengths ($(length(sequence1)) ≠ $(length(sequence2))). ProtoSyn will only consider the range of the smallest sequence."
    end

    s_total = 0.0
    N = minimum((length(sequence1), length(sequence2)))
    for i in 1:N
        aa1 = string(sequence1[i])
        aa2 = string(sequence2[i])
        s = Peptides.aminoacid_similarity[aa1][aa2]
        s_total += s
        @printf "%s - %s : %5.3f\n" aa1 aa2 s
    end
    s = @sprintf "Total similarity value: %8.3f | " s_total
    s *= @sprintf "Average similarity value: %8.3f\n" s_total / N
    printstyled(s, color = :green)

    return s_total, s_total / N
end