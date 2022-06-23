using Random
using ProtoSyn
using ProtoSyn.Drivers
using ProtoSyn.Peptides
using ProtoSyn.Calculators

"""
    RotamerBlitzState{T <: AbstractFloat}(step::Int = 0, converged::Bool = false, completed::Bool = false, stalled::Bool = false, acceptance_count = 0)
    RotamerBlitzState{T <: AbstractFloat}()

A [`RotamerBlitzState`](@ref) instance is a `DriverState` that keeps track of a
[`RotamerBlitz`](@ref) simulation status. Besides the default `:step`,
`:converged`, `:completed` and `:stalled` fields, a [`RotamerBlitzState`](@ref)
does not add any additional tracking fields specific for [`RotamerBlitz`](@ref)
simulations.

# See also
[`RotamerBlitz`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Peptides.Drivers.RotamerBlitzState()
ProtoSyn.Peptides.Drivers.RotamerBlitzState(0, false, false, false)
```
"""
Base.@kwdef mutable struct RotamerBlitzState <: ProtoSyn.Drivers.DriverState
    step::Int        = 0
    converged::Bool  = false
    completed::Bool  = false
    stalled::Bool    = false
end


"""
    RotamerBlitz(eval!::Union{Function, ProtoSyn.Calculators.EnergyFunction}, rotamer_library::Dict{String, ProtoSyn.Peptides.BBD_RotamerLibrary}, n_first::Int, max_steps::Int, callback::Opt{ProtoSyn.Drivers.Callback}, selection::Opt{AbstractSelection})

A [`RotamerBlitz`](@ref) `Driver` instance. As such, this object is callable as
a _functor_ with the following signature:

```julia
(driver::RotamerBlitz)(pose::Pose)
```

This simulation performs a sidechain packaging algorithm, where each step (until
a maximum of `n_steps`), all [`Residue`](@ref) instances in the provided
[`Pose`](@ref) `pose` (in random order) may suffer a rotamer conformation
change. As such, for each [`Residue`](@ref) instance, all `n_first` most likely
[`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) instances in the `rotamer_library`
(according to the current backbone phi and psi dihedral angles) are applied
(using the [`apply!`](@ref ProtoSyn.Peptides.apply!) method) and evaluated by
the provided `eval!`
[`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) or custom function.
Once all `n_first` most likely [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer)
instances are looped over, the most favourable (least energetic) one is
re-applied. Each step, a call to an optional
[`Callback`](@ref ProtoSyn.Drivers.Callback) `callback` is performed. A
companion [`RotamerBlitzState`](@ref) `DriverState` instance is also updated
each step and provided to the [`Callback`](@ref ProtoSyn.Drivers.Callback)
`callback`.

# Fields
* `eval!::Union{Function, EnergyFunction}` - The evaluator [`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) or custom function, receives a [`Pose`](@ref) `pose` as the single argument;
* `rotamer_library::Dict{String, ProtoSyn.Peptides.BBD_RotamerLibrary}` - The rotamer library used to sample new [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) instances from;
* `n_first::Int` - Maximum number of [`Rotamer`](@ref ProtoSyn.Peptides.Rotamer) instances to try, ordered from most to less likely;
* `max_steps::Int` - The total number of simulation steps to be performed;
* `callback::Opt{Callback}` - An optional [`Callback`](@ref ProtoSyn.Drivers.Callback) instance, receives two input arguments: the current [`Pose`](@ref) `pose` and the current `DriverState` `driver_state`;
* `selection::Opt{AbstractSelection}` - If provided, only selected [`Residue`](@ref) instances will be considered as target for the [`RotamerBlitz`](@ref) simulation, otherwise, consider all [`Residue`](@ref) instances.

# See also
[`RotamerBlitzState`](@ref)

# Examples
```
julia> rb = ProtoSyn.Peptides.Drivers.RotamerBlitz(energy_function, rot_lib, 1, 2, cb, an"CA")
⚒  Rotamer Blitz Driver:
 ├──  ●  Evaluator:
 |    └── 🗲  Energy Function (4 components):
 |        +----------------------------------------------------------------------+
 |        | Index | Component name                                | Weight (α)   |
 |        +----------------------------------------------------------------------+
 |        | 1     | TorchANI_ML_Model                             |      1.000   |
 |        | 2     | Caterpillar_Solvation                         |      0.010   |
 |        | 3     | Bond_Distance_Restraint                       |      1.000   |
 |        | 4     | Cα-Cα_Clash_Restraint                         |    100.000   |
 |        +----------------------------------------------------------------------+
 |   
 ├──  ● Callback:
 |    └── ✉  Callback:
 |        +----------------------------------------------------------------------+
 |        | Index | Field                     | Value                            |
 |        +----------------------------------------------------------------------+
 |        | 1     | Event                     | energy_step_frame                |
 |        | 2     | Frequency                 | 10                               |
 |        +----------------------------------------------------------------------+
 |   
 ├──  ● Selection:
 |    └── FieldSelection › Atom.name = CA
 |   
 └──  ● Settings:
             Rotamer Library: Set ✓
                   Max steps: 2
     Search N first rotamers: 1
```
"""
mutable struct RotamerBlitz <: ProtoSyn.Drivers.Driver
    eval!::Union{Function, ProtoSyn.Calculators.EnergyFunction}
    rotamer_library::Dict{String, ProtoSyn.Peptides.BBD_RotamerLibrary}
    n_first::Int
    max_steps::Int
    callback::Opt{ProtoSyn.Drivers.Callback}
    selection::Opt{AbstractSelection}
end


function (driver::RotamerBlitz)(pose::Pose)
    
    T = eltype(pose.state)
    driver_state = RotamerBlitzState()
    if driver.selection === nothing
        residues = collect(eachresidue(pose.graph))
    else
        residues = ProtoSyn.promote(driver.selection, Residue)(pose, gather = true)
    end

    driver.callback !== nothing && driver.callback(pose, driver_state)

    for step in 1:driver.max_steps
        for residue in shuffle(residues)
            
            # Skip if residue does not have a rotamer library entry
            !(residue.name in keys(driver.rotamer_library)) && continue

            # Residue results hold the energy of each rotamer conformation
            # The best rotamer conformation will be applied to this residue
            energy          = driver.eval!(pose)
            rotamer         = ProtoSyn.Peptides.get_rotamer(pose, residue)
            residue_results = Dict{ProtoSyn.Peptides.Rotamer, T}(
                rotamer => energy)

            # Get the correct rotamer stack
            phi = getdihedral(pose.state, Peptides.phi(residue))
            psi = getdihedral(pose.state, Peptides.psi(residue))
            rotamer_stack = driver.rotamer_library[residue.name][phi, psi]

            # Iterate the 'n_first' rotamers, saving the evaluated energy
            _n_first = min(length(rotamer_stack.rotamers), driver.n_first)
            for rotamer in rotamer_stack.rotamers[1:_n_first]
                # * Note: Peptides.apply! doesn't sync! the pose
                exact_rotamer = Peptides.apply!(pose.state, rotamer, residue)
                sync!(pose)
                energy = driver.eval!(pose)
                residue_results[exact_rotamer] = energy
            end

            # Re-apply the best conformation (exactly)
            _, rotamer = findmin(residue_results)
            Peptides.apply!(pose.state, rotamer, residue)
        end

        driver_state.step = step
        driver.callback !== nothing && driver.callback(pose, driver_state) 
    end
    
    driver_state.completed = true
    driver_state
end

# * Show -----------------------------------------------------------------------
function Base.show(io::IO, rbd::RotamerBlitz, level_code::Opt{LevelCode} = nothing)
    init_level_code = level_code === nothing ? LevelCode() : level_code
    init_lead       = ProtoSyn.get_lead(level_code)
    init_inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, init_lead*"⚒  Rotamer Blitz Driver:")

    level_code = vcat(init_level_code, 3)
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    if typeof(rbd.eval!) === ProtoSyn.Calculators.EnergyFunction
        println(io, lead*" ●  Evaluator:")
        Base.show(io, rbd.eval!, vcat(level_code, 4))
    else
        println(io, lead*" ●  Evaluator: $(rbd.eval!)")
    end

    println(io, init_inner_lead*level_code.code_table[1])

    if rbd.callback !== nothing
        println(io, lead*" ● Callback:")
        Base.show(io, rbd.callback, vcat(level_code, 4))
    else
        println(io, lead*" ○  Callback: Not set")
    end

    println(io, init_inner_lead*level_code.code_table[1])

    if rbd.selection !== nothing
        println(io, lead*" ● Selection:")
        Base.show(io, rbd.selection, vcat(level_code, 4))
    else
        println(io, lead*" ○  Selection: Not set")
    end

    level_code = vcat(init_level_code, 4)
    lead       = ProtoSyn.get_lead(level_code)
    inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, init_inner_lead*level_code.code_table[1])

    println(io, lead*" ● Settings:")
    println(io, inner_lead*"        Rotamer Library: Set ✓")
    println(io, inner_lead*"              Max steps: $(rbd.max_steps)")
    println(io, inner_lead*"Search N first rotamers: $(rbd.n_first)")
end