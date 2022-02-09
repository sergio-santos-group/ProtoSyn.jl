using Printf

function measure_similarity(sequence1::String, sequence2::String; show_results::Bool = true)
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