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
function evaluate!(st::Common.State, solv_pairs::Vector{SolvPair}, do_forces::Bool = false)::Float64
    n_res              = length(solv_pairs)
    l                  = [.0, .0, .0]
    e_sol              = .0
    sum_f              = .0
    t                  = 21.0
    distance_threshold = 1.2

    for i in 1:n_res
        ci = solv_pairs[i].coef
        sum_f = .0
        for j in 1:n_res
            if i == j
                continue
            end
            dIJ = .0
            @inbounds for k=1:3
                delta = st.xyz[solv_pairs[i].i, k] - st.xyz[solv_pairs[j].i, k]
                dIJ  += delta ^ 2
            end
            dIJ = sqrt(dIJ)
            sum_f += 1.0 / (1.0 + exp(-(distance_threshold - dIJ) / 0.4))
        end
        if ((sum_f < t) && (ci > .0)) || ((sum_f > t) && (ci < .0))
            e_sol += ci * (t - sum_f)
        end
    end

    Common.set_energy_component(st.energy, :sol, e_sol)
    return e_sol
end


@doc raw"""
    evaluate!(hb_networks::HbNetwork, state::Common.State[, do_forces::Bool = false])::Float64

Evaluate a [`HbNetwork`](@ref)'s using the current [`Common.State`](@ref),
calculate and update state.energy according to the coarse-grain model function defined.

```math
eH = - \sum_{i = 1}^{n_{i}} \sum_{j = 1}^{n_{i}}  -[cos(\theta _{1})cos(\theta _{2})]^{2}\times \left
[ 5.0 \left (\frac{0.2}{d_{OH}} \right )^{12} - 6.0 \left (\frac{0.2}{d_{OH}} \right )^{10}\right ]
```

- `dOH` : Distance between the hydrogen of residue `i` and oxygen of residue `j`;
- `θ1`  : Angle between NĤO of residue `i`;
- `θ2`  : Angle between OĈCα of residue `j`;

# Examples
```julia-repl
julia> Forcefield.CoarseGrain.evaluate!(hb_network, state)
-0.500
```
"""
function evaluate!(st::Common.State, hb_network::HbNetwork, do_forces::Bool = false)::Float64

    vHN    = [.0, .0, .0]
    vOC    = [.0, .0, .0]
    vHO    = [.0, .0, .0]
    c1     = .0
    c2     = .0
    e_H    = .0 
    coords = st.xyz

    for donor in hb_network.donors
        for acceptor in hb_network.acceptors
            if acceptor.base == donor.base
                continue
            end

            dhn_ho = .0
            doc_ho = .0
            dho    = .0
            @inbounds for i=1:3
                delta_hn = coords[donor.base, i]       - coords[donor.charged, i]
                delta_ho = coords[acceptor.charged, i] - coords[donor.charged, i]
                delta_oc = coords[acceptor.base, i]    - coords[acceptor.charged, i]
                dhn_ho  += delta_hn*delta_ho
                doc_ho  += delta_oc*delta_ho
                dho     += delta_ho ^ 2
            end
            dho = sqrt(dho)
            c1 = dhn_ho / dho
            c2 = doc_ho / dho
            e_H -= - ((c1 * c2)^2) * (5.0 * ((0.2 / dho)^12) - 6.0 * ((0.2 / dho)^10))
        end
    end

    e_H *= hb_network.coef
    Common.set_energy_component(st.energy, :hb, e_H)
    return e_H
end