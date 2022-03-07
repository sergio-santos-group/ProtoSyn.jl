module Calculators

    using ProtoSyn
    using Base.Cartesian
    using Printf

    MaskMap = Opt{Union{ProtoSyn.Mask{<: ProtoSyn.AbstractContainer}, Matrix{<: AbstractFloat}, Function}}

    include("verlet_list.jl")
    include("distance_matrix.jl")

    # Load energy function components
    include("energy_function_component.jl")

    if "NO_TORCHANI" in keys(ENV) && ENV["NO_TORCHANI"] === "true"
        @warn "Environment variable NO_TORCHANI set to `true`. Not loading torchani."
    else
        @info " | Loading TorchANI"
        include("torchani.jl")
    end

    @info " | Loading Hydrogen Bonds"
    include("hydrogen_bonds.jl")

    @info " | Loading SASA"
    include("sasa.jl")

    @info " | Loading GB"
    include("gb.jl")

    @info " | Loading Electrostatics"
    include("electrostatics.jl")

    @info " | Loading Restraint Models"
    include("Potentials/potentials.jl")
    include("restraints.jl")

    @info " | Loading Energy Function"
    include("energy_function.jl")

    include("ref15.jl")

end