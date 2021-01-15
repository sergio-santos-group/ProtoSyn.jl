using Printf

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
            sync!(pose) # apply last driver changes
            driver(pose)
        end
    end
end

Base.push!(compound_driver::CompoundDriver, driver::Union{Function, AbstractMutator, Driver}) = begin
    push!(compound_driver.drivers, driver)
end

function Base.show(io::IO, drv::CompoundDriver)
    println(" Compound Driver ($(length(drv.drivers)) elements):")
    for (i, driver) in enumerate(drv.drivers)
        Base.show(driver)
    end
end