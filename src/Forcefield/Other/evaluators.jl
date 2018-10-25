function calc_eSol!(st::Common.State, phi_dihedrals::Vector{Common.Dihedral})

    # Verify input
    for dihedral in phi_dihedrals
        dihedral.dtype == Common.phi ? nothing : error("Tried to evaluate eSol but dihedral $dihedral is not a PHI.")
    end
    
    # Create solvation pairs
    d = Dict("Q" => -3.5, "W" => -0.9, "T" => -0.7, "C" =>  2.5, "P" => -1.6, "V" =>  4.2, "L" =>  3.8, "M" =>  1.9, "N" => -3.5, "H" => -3.2,
             "A" =>  1.8, "D" => -3.5, "G" => -0.4, "E" => -3.5, "Y" => -1.3, "I" =>  4.5, "S" => -0.8, "K" => -3.9, "R" => -4.5, "F" =>  2.8)
    solv_pairs = map(x -> SolvPair(x.a3, d[Aux.conv_aa_321(string(x.residue.name))]), phi_dihedrals)
    
    # Calculate solvation energy
    n_atoms = length(solv_pairs)
    l = zeros(Float64, 3)
    e_sol::Float64 = 0.0
    sum_f::Float64 = 0.0
    for i in 1:(n_atoms - 1)
        Ai = @view st.xyz[solv_pairs[i].i,:]
        sum_f = 0.0
        for j in (i+1):n_atoms
            Aj = @view st.xyz[solv_pairs[j].i,:]
            @. l[:] = Aj - Ai
            dIJ = norm(l)
            f = 1.0 / (1.0 + exp(-(12.0-dIJ) / 0.4))
            sum_f += f
        end
        if ((sum_f < 21.0) && (solv_pairs[i].coef > 0.0)) || ((sum_f > 21.0) && (solv_pairs[i].coef < 0.0))
            e_sol += solv_pairs[i].coef * (21.0 - sum_f)
        end
    end

    st.energy.comp["eSol"] = e_sol
    st.energy.eTotal += e_sol
    return e_sol
end