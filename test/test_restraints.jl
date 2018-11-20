using ProtoSyn
using Printf


fbr_union = Union{Forcefield.Restraints.DistanceFBR, Forcefield.Restraints.DihedralFBR}
@doc raw"""
    calc_numeric_forces(state::Common.State, fbr::Union{Forcefield.Restraints.DistanceFBR, Forcefield.Restraints.DihedralFBR}[, e1::Float64 = state.energy.eTotal, ϵ::Float64 = 1e-8]::Array{Float64, 2}

Calculate forces applied in the system to accomodate the given flat-bottomed restrain, using numeric derivatives.

# Examples
```julia-repl
julia> calc_numeric_forces(state, fbr)
[1.0 1.0 1.0; -1.0 -1.0 -1.0]
```
See also: [`test_distanceFBR`](@ref) [`test_dihedralFBR`](@ref)
"""
function calc_numeric_forces(state::Common.State, fbr::fbr_union, e1::Float64, ϵ::Float64 = 1e-8)::Array{Float64, 2}

    f2::Array{Float64, 2} = zeros(size(state.forces, 1), size(state.forces, 2))
    for atom_index in 1:size(state.xyz, 1)
        # Forces need to be measured and compared by moving one component at a time
        for component in [1, 2, 3]
            backup_component = state.xyz[atom_index, component]
            state.xyz[atom_index, component] += ϵ
            e2::Float64 = Forcefield.Restraints.evaluate!([fbr], state, do_forces=false)
            f2[atom_index, component] = - (e2 - e1) / ϵ
            state.xyz[atom_index, component] = backup_component
        end
    end
    return f2
end


@doc raw"""
    test_distanceFBR([workload::Int64 = 100, init_distance::Float64 = 10.0])

Calculate energy and forces applied in the system given a default flat-bottomed restraint (both using analytical and numeric derivatives);
A total of `workload` distances will be tested, equally spaced between `init_distance` and 0.0;
Print the energy, max_force felt and if the max_force is bellow a significant value (Default: 1e-3) for every distance tested.

# Examples
```julia-repl
julia> test_distanceFBR(2)
----------------------------------
         Distance FBR Test
Distance    Energy Max Force  Pass
----------------------------------
   10.00  6.00e+03  1.31e-04  true
    5.00  0.00e+00 -0.00e+00  true
    0.00  6.00e+03       NaN false
----------------------------------
```
See also: [`calc_numeric_forces`](@ref) [`test_dihedralFBR`](@ref)
"""
function test_distanceFBR(workload::Int64 = 100, init_distance::Float64 = 10.0)

    state = Common.State(2)
    state.xyz = [0.0 0.0 0.0; init_distance 0.0 0.0]
    fbr = Forcefield.Restraints.DistanceFBR(1, 2, 2.0, 4.0, 6.0, 8.0, 1e3)
    
    @printf "%s Distance FBR Test\n%8s %9s %9s %5s\n%s\n" " "^8 "Distance" "Energy" "Max Force" "Pass" "-"^34
    for d in init_distance:-(init_distance/workload):0.0
        state.xyz[2, 1] = d
        fill!(state.forces, zero(Float64))
        e1::Float64 = Forcefield.Restraints.evaluate!([fbr], state, do_forces=true)
        f2::Array{Float64, 2} = calc_numeric_forces(state, fbr, e1)
        @printf "%8.2f %9.2e %9.2e %5s\n" d e1 maximum(f2-state.forces) string(maximum(f2 - state.forces) < 1e-3)
    end
    @printf "%s\n" "-"^34
end


@doc raw"""
    test_dihedralFBR([workload::Int64 = 100])

Calculate energy and forces applied in the system given a default flat-bottomed restraint (both using analytical and numeric derivatives);
A total of `workload` angles will be tested, equally spaced between -180.0° and 180.0°;
Print the energy, max_force felt and if the max_force is bellow a significant value (Default: 1e-3) for every angle tested.

# Examples
```julia-repl
julia> test_dihedralFBR()
----------------------------------
         Dihedral FBR Test
Angle    Energy Max Force  Pass
---------------------------------
-180.00  1.54e+03  1.57e+03 false
-170.00  1.41e+03  4.56e-05  true
-160.00  1.27e+03  1.19e-05  true
(...)
----------------------------------
```
See also: [`calc_numeric_forces`](@ref) [`test_distanceFBR`](@ref)
"""
function test_dihedralFBR(workload::Int64 = 360)

    state = Common.State(4)
    state.xyz = [0.0 -1.0 0.0; 0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0]
    fbr = Forcefield.Restraints.DihedralFBR(1, 2, 3, 4, deg2rad(-90), deg2rad(-45), deg2rad(45), deg2rad(90), 1e3)
    f2 = zeros(Float64, state.size, state.size)
    
    @printf "%s Dihedral FBR Test\n%7s %9s %9s %5s\n%s\n" " "^8 "Angle" "Energy" "Max Force" "Pass" "-"^33
    for θ in -180:(360/workload):170.0
        θ = deg2rad(θ)
        state.xyz[1, 2] = cos(θ)
        state.xyz[1, 3] = sin(θ)
        fill!(state.forces, zero(Float64))
        e1::Float64 = Forcefield.Restraints.evaluate!([fbr], state, do_forces=true)
        f2::Array{Float64, 2} = calc_numeric_forces(state, fbr, e1)

        @printf "%7.2f %9.2e %9.2e %5s\n" rad2deg(θ) e1 maximum(f2-state.forces) string(maximum(f2 - state.forces) < 1e-3)
    end
    @printf "%s\n" "-"^34
    println(f2, " ", state.forces)
end

# ---------------------------------------------------------------------------
# NOTE: Running this file automatically runs both DistanceFBR and DihedralFBR test;
# ---------------------------------------------------------------------------

# test_distanceFBR()
test_dihedralFBR()