@doc raw"""
    evaluate!(solv_pairs::Vector{SolvPair}, state::Common.State[, do_forces::Bool = false])::Float64

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
    n_res                       = length(solv_pairs)
    l::Vector{Float64}          = zeros(Float64, 3)
    e_sol::Float64              = 0.0
    sum_f::Float64              = 0.0
    t::Float64                  = 21.0
    distance_threshold::Float64 = 1.2

    for i in 1:n_res
        Ai = @view st.xyz[solv_pairs[i].i, :]
        sum_f = 0.0
        for j in 1:n_res
            if i == j
                continue
            end
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


@doc raw"""
    evaluate!(hb_groups::Vector{HbGroup}, state::Common.State[, do_forces::Bool = false])::Float64

Evaluate an array of [`HbGroup`](@ref)'s using the current [`Common.State`](@ref),
calculate and update state.energy according to the coarse-grain model function defined.

```math
eH = - \sum_{i = 1}^{n_{i}} \sum_{j = 1}^{n_{i}}  -3.1[cos(\theta _{1})cos(\theta _{2})]^{2}\times \left
[ 5.0 \left (\frac{0.2}{d_{OH}} \right )^{12} - 6.0 \left (\frac{0.2}{d_{OH}} \right )^{10}\right ]
```

- `dOH` : Distance between the hydrogen of residue `i` and oxygen of residue `j`;
- `θ1`  : Angle between NĤO of residue `i`;
- `θ2`  : Angle between OĈCα of residue `j`;

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.evaluate!(hb_groups, state)
-0.500
```
"""
function evaluate!(hb_groups::Vector{HbGroup}, st::Common.State; do_forces::Bool = false)::Float64

    n_res::Int64         = length(hb_groups)
    vHN::Vector{Float64} = zeros(Float64, 3)
    vOC::Vector{Float64} = zeros(Float64, 3)
    vHO::Vector{Float64} = zeros(Float64, 3)
    c1::Float64          = 0.0
    c2::Float64          = 0.0
    e_H::Float64         = 0.0 

    for i in 1:(n_res)
        N = @view st.xyz[hb_groups[i].n, :]
        H = @view st.xyz[hb_groups[i].h, :]
        @. vHN[:] = N - H

        for j in 1:n_res
            if i==j
                continue
            end
            C = @view st.xyz[hb_groups[j].c, :]
            O = @view st.xyz[hb_groups[j].o, :]
            @. vOC[:] = C - O
            @. vHO[:] = O - H
            dHO = norm(vHO)

            c1 = dot(vHN, vHO) / dHO
            c2 = dot(vOC, vHO) / dHO
            e_H -= -3.1 * ((c1 * c2)^2) * (5.0 * ((0.2 / dHO)^12) - 6.0 * ((0.2 / dHO)^10)) * hb_groups[i].coef
        end
    end

    st.energy.comp["eH"] = e_H
    st.energy.eTotal = e_H
    return e_H
end