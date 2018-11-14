#TODO: Document function
# function calc_eSol!(st::Common.State, residues::Vector{Common.Residue})
    
#     # Create solvation pairs
#     d = Dict("Q" => -3.5, "W" => -0.9, "T" => -0.7, "C" =>  2.5, "P" => -1.6, "V" =>  4.2, "L" =>  3.8, "M" =>  1.9, "N" => -3.5, "H" => -3.2,
#              "A" =>  1.8, "D" => -3.5, "G" => -0.4, "E" => -3.5, "Y" => -1.3, "I" =>  4.5, "S" => -0.8, "K" => -3.9, "R" => -4.5, "F" =>  2.8)
#     solv_pairs = map(x -> SolvPair(x.cα, d[Aux.conv321(string(x.name))]), residues)
    
#     # Calculate solvation energy
#     n_atoms = length(solv_pairs)
#     l = zeros(Float64, 3)
#     e_sol::Float64 = 0.0
#     sum_f::Float64 = 0.0
#     for i in 1:(n_atoms - 1)
#         Ai = @view st.xyz[solv_pairs[i].i,:]
#         sum_f = 0.0
#         for j in (i+1):n_atoms
#             Aj = @view st.xyz[solv_pairs[j].i,:]
#             @. l[:] = Aj - Ai
#             dIJ = norm(l)
#             f = 1.0 / (1.0 + exp(-(12.0-dIJ) / 0.4))
#             sum_f += f
#         end
#         if ((sum_f < 21.0) && (solv_pairs[i].coef > 0.0)) || ((sum_f > 21.0) && (solv_pairs[i].coef < 0.0))
#             e_sol += solv_pairs[i].coef * (21.0 - sum_f)
#         end
#     end

#     st.energy.comp["eSol"] = e_sol
#     st.energy.eTotal = e_sol
#     return e_sol
# end


# function evaluate!(contact_pairs::Vector{ContactPair}, st::Common.State; threshold::Float64 = 0.5, k::Float64 = 1.0)

#     eContact::Float64 = 0.0
#     v12 = zeros(Float64, 3)
#     d12Sq = 0.0
#     tSq = threshold^2
#     for pair in contact_pairs
#         @views @. v12 = st.xyz[pair.c1, :] - st.xyz[pair.c2, :]
#         d12Sq = dot(v12, v12)
#         if d12Sq > tSq
#             eContact += pair.prob * (sqrt(d12Sq) - threshold)^2
#         end
#     end

#     eContact *= k
#     st.energy.comp["eContact"] = eContact
#     st.energy.eTotal = eContact
#     return eContact
# end

function evaluate!(contact_topology::Vector{ContactPair}, st::Common.State; r1::Float64 = 0.5, r2::Float64 = 0.8, k::Float64 = 1.0, do_forces::Bool = false)
    # All distances are in nm

    eContact::Float64 = 0.0
    v12 = zeros(Float64, 3)
    e2 = (r2 - r1) * (r2 - r1)
    for pair in contact_topology
        @views @. v12 = st.xyz[pair.c1, :] - st.xyz[pair.c2, :]
        d12 = norm(v12)
        if r1 <= d12 < r2
            dr = (d12 - r1)
            eContact += pair.prob * dr * dr
            if do_forces
                @. v12 *= k * (pair.prob * dr / d12)
                @. st.forces[pair.r1, :] -= v12'
                @. st.forces[pair.r2, :] += v12'
            end
        elseif r2 <= d12
            eContact += pair.prob * (e2 + d12 - r1)
            if do_forces
                @. v12 *= k * pair.prob * 0.5
                @. st.forces[pair.r1, :] -= v12'
                @. st.forces[pair.r2, :] += v12'
            end
        end
    end

    eContact *= k * 0.5
    st.energy.comp["eContact"] = eContact
    st.energy.eTotal = eContact
    return eContact
end