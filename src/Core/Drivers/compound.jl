using Printf

"""
    CompoundDriver(drivers::Vector{Union{Function, AbstractMutator, Driver}})

A [`CompoundDriver`](@ref) `Driver` instance. As such, this object is callable
as a _functor_ with the following signature:

```julia
(driver::CompoundDriver)(pose::Pose)
```

A [`CompoundDriver`](@ref) `Driver` groups several `Driver`, `AbstractMutator`
and custom functions, calling each of the components sequentially. Does not
[`sync!`](@ref) changes, this task is left to each individual `Driver` or
`AbstractMutator` instance.

!!! ukw "Note:"
    In contrast with [`CompoundMutator`](@ref ProtoSyn.Mutators.CompoundMutator),
    this `Driver` does not accept an `AbstractSelection` modifier. This may
    change in future versions of ProtoSyn.

# Examples
```jldoctest
julia> ProtoSyn.Drivers.CompoundDriver([monte_carlo, monte_carlo])
⚒  CompoundDriver Driver (2 elements):
 ├── ⚒  Monte Carlo Driver:
 |    ├──  ●  Evaluator:
 |    |    └── 🗲  Energy Function (4 components):
 |    |        +----------------------------------------------------------------------+
 |    |        | Index | Component name                                | Weight (α)   |
 |    |        +----------------------------------------------------------------------+
 |    |        | 1     | TorchANI_ML_Model                             |      1.000   |
 |    |        | 2     | Caterpillar_Solvation                         |      0.010   |
 |    |        | 3     | Bond_Distance_Restraint                       |      1.000   |
 |    |        | 4     | Cα-Cα_Clash_Restraint                         |    100.000   |
 |    |        +----------------------------------------------------------------------+
 |    |   
 |    ├──  ● Sampler:
 |    |    └── ⚯  Dihedral Mutator:
 |    |        +----------------------------------------------------------------------+
 |    |        | Index | Field                       | Value                          |
 |    |        +----------------------------------------------------------------------+
 |    |        | 1     | angle_sampler               | Function randn                 |
 |    |        | 2     | p_mut                       | 0.0100                         |
 |    |        | 3     | step_size                   | 0.5000                         |
 |    |        +----------------------------------------------------------------------+
 |    |         ● Selection: Set
 |    |         └── FieldSelection › Atom.name = r"C|N"
 |    |   
 |    ├──  ● Callback:
 |    |    └── ✉  Callback:
 |    |        +----------------------------------------------------------------------+
 |    |        | Index | Field                     | Value                            |
 |    |        +----------------------------------------------------------------------+
 |    |        | 1     | Event                     | energy_step                      |
 |    |        | 2     | Frequency                 | 1                                |
 |    |        +----------------------------------------------------------------------+
 |    |   
 |    └──  ● Settings:
 |          Max steps: 10
 |        Temperature: linear_quench
 |   
 └── ⚒  Monte Carlo Driver:
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
      ├──  ● Sampler:
      |    └── ⚯  Dihedral Mutator:
      |        +----------------------------------------------------------------------+
      |        | Index | Field                       | Value                          |
      |        +----------------------------------------------------------------------+
      |        | 1     | angle_sampler               | Function randn                 |
      |        | 2     | p_mut                       | 0.0100                         |
      |        | 3     | step_size                   | 0.5000                         |
      |        +----------------------------------------------------------------------+
      |         ● Selection: Set
      |         └── FieldSelection › Atom.name = r"C|N"
      |   
      ├──  ● Callback:
      |    └── ✉  Callback:
      |        +----------------------------------------------------------------------+
      |        | Index | Field                     | Value                            |
      |        +----------------------------------------------------------------------+
      |        | 1     | Event                     | energy_step                      |
      |        | 2     | Frequency                 | 1                                |
      |        +----------------------------------------------------------------------+
      |   
      └──  ● Settings:
            Max steps: 10
          Temperature: linear_quench
```
"""
mutable struct CompoundDriver <: Driver
    drivers::Vector{Union{Function, AbstractMutator, Driver}}
end

function (compound_driver::CompoundDriver)(pose::Pose)
    # Note: as any other mutator/driver, a CoumpundDriver doesn't sync! the pose
    # before applying conformational changes, using the system "as is". This is
    # the default behaviour, by design, since sequential changes might not need
    # to sync! between themselves. For example, two or more DihedralMutators in
    # sequence simply change the internal coordinates, who don't need to be
    # applied and calculated to cartesian coordinates. However, given the fact
    # that this struct applies multiple changes of unknown nature, a sync! is
    # performed between drivers/mutators. For example: a DihedralMutator
    # followed by a CrankshaftMutator. The DihedralMutator issues an internal to
    # cartesian synchornization, who needs to be applied since the
    # CrankshaftMutator uses cartesian coordinates to calculate the rotation
    # axis.

    compound_driver.drivers[1](pose)
    if length(compound_driver.drivers) > 1
        for driver in compound_driver.drivers[2:end]
            driver(pose)
        end
    end
    
    return pose
end

Base.push!(compound_driver::CompoundDriver, driver::Union{Function, AbstractMutator, Driver}) = begin
    push!(compound_driver.drivers, driver)
end

# * Show -----------------------------------------------------------------------

function Base.show(io::IO, drv::CompoundDriver, level_code::Opt{LevelCode} = nothing)
    init_level_code = level_code === nothing ? LevelCode() : level_code
    init_lead       = ProtoSyn.get_lead(level_code)
    init_inner_lead = ProtoSyn.get_inner_lead(level_code)

    println(io, init_lead*"⚒  CompoundDriver Driver ($(length(drv.drivers)) elements):")

    for driver in drv.drivers[1:(end-1)]
        if isa(driver, AbstractMutator)
            Base.show(io, driver, vcat(init_level_code, 3))
            println(io, init_inner_lead*init_level_code.code_table[1])
        elseif isa(driver, Driver)
            Base.show(io, driver, vcat(init_level_code, 3))
            println(io, init_inner_lead*init_level_code.code_table[1])
        else
            f_lead = ProtoSyn.get_lead(vcat(init_level_code, 3))
            println(io, f_lead*" ●  Function: $driver")
            println(io, init_inner_lead*init_level_code.code_table[1])
        end
    end

    if isa(drv.drivers[end], AbstractMutator)
        Base.show(io, drv.drivers[end], vcat(init_level_code, 4))
    elseif isa(drv.drivers[end], Driver)
        Base.show(io, drv.drivers[end], vcat(init_level_code, 4))
    else
        f_lead = ProtoSyn.get_lead(vcat(init_level_code, 4))
        println(io, f_lead*" ●  Function: $(drv.drivers[end])")
    end
end