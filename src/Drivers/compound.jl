mutable struct CompoundDriver <: Driver
    drivers::Vector{Union{Function, AbstractMutator, Driver}}
end

function (compound_driver::CompoundDriver)(pose::Pose)
    for driver in compound_driver.drivers
        driver(pose)
    end
end

Base.push!(compound_driver::CompoundDriver, driver::Union{Function, AbstractMutator, Driver}) = begin
    push!(compound_driver.drivers, driver)
end