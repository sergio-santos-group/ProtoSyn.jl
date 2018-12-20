@doc raw"""
    evaluate!(solv_pairs::Vector{SolvPair}, state::Common.State[, do_forces::Bool = false, cut_off::Float64 = 2.0])::Float64

Evaluate an array of [`SolvPair`](@ref)'s using the current [`Common.State`](@ref),
calculate and update state.energy according to the coarse-grain model function defined.

```math
eSol = \sum_{i = 1}^{n_{i}}
\left\{\begin{matrix}
c_{i}\times (t - \sum f_{i}) & ~if~ (\sum f_{i} < t ~and~ c_{i}>0) ~or~ (\sum f_{i} > t ~and~ c_{i}<0) \\ 
0 & ~if~(\sum f_{i} > t ~and~ c_{i}<0) ~or~ (\sum f_{i} < t ~and~ c_{i}>0)
\end{matrix}\right.
```
```math
where \sum f_{i} = \sum_{j = i}^{n_{i}}\frac{1}{1+e^{-\frac{d_{t}-d_{ij}}{0.4}}}
```

- `dt` : Distance threshold (1.2 by default);
- `t`  : Minimum neighbouring residues (21.0 by default);
- `ci` : Residue solvation coeficient.

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.evaluate!(solv_pairs, state)
0.500
```
"""
function evaluate!(solv_pairs::Vector{SolvPair}, st::Common.State; do_forces::Bool = false)::Float64
    n_atoms = length(solv_pairs)
    l = zeros(Float64, 3)
    e_sol::Float64 = 0.0
    sum_f::Float64 = 0.0
    t::Float64 = 21.0
    distance_threshold::Float64 = 1.2
    for i in 1:(n_atoms - 1)
        Ai = @view st.xyz[solv_pairs[i].i, :]
        sum_f = 0.0
        for j in (i+1):n_atoms
            Aj = @view st.xyz[solv_pairs[j].i, :]
            @. l[:] = Aj - Ai
            dIJ = norm(l)
            f = 1.0 / (1.0 + exp(-(distance_threshold - dIJ) / 0.4))
            sum_f += f
        end
        if ((sum_f < t) && (solv_pairs[i].coef > 0.0)) || ((sum_f > t) && (solv_pairs[i].coef < 0.0))
            e_sol += solv_pairs[i].coef * (t - sum_f)
        end
    end
    st.energy.comp["eSol"] = e_sol
    st.energy.eTotal = e_sol
    return e_sol
end