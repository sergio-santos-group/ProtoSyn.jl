@doc """
    The Common module stores functions with higher levels of abstraction and
    that are common practice in ProtoSyn algorithms.
"""
module Common
    using Printf
    using ProtoSyn
    using ProtoSyn.Drivers: DriverState, Callback
    using ProtoSyn.Units: defaultFloat
    using ProtoSyn.Calculators
    using ProtoSyn.Peptides

    include("energy_function.jl")
    include("callback.jl")
end