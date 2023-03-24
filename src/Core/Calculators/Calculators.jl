module Calculators

    using ProtoSyn
    using Base.Cartesian
    using Printf

    MaskMap = Opt{Union{ProtoSyn.Mask{<: ProtoSyn.AbstractContainer}, Matrix{<: AbstractFloat}, Function}}

    include("verlet_list.jl")
    include("distance_matrix.jl") # Requires verlet_list.jl

    # Load energy function components
    include("energy_function_component.jl")

    if "USE_TORCHANI" in keys(ENV) && ENV["USE_TORCHANI"] === "false"
        @warn "Environment variable USE_TORCHANI set to `false`. Not loading torchani."
    else
        @info " | Loading TorchANI"
        include("torchani.jl")
    end

    @info " | Loading Hydrogen Bonds"
    include("hydrogen_bonds.jl")

    @info " | Loading SASA"
    include("sasa.jl")
    include("radius_gyration.jl")

    if "USE_IGBR_NN" in keys(ENV) && ENV["USE_IGBR_NN"] === "false"
        @warn "Environment variable USE_IGBR_NN set to `false`. Not loading iGBR-NN ONNX models."
    else
        @info " | Loading GB"
        include("gb.jl")
    end

    @info " | Loading Electrostatics"
    include("electrostatics.jl")

    @info " | Loading Restraint Models"
    include("Potentials/potentials.jl")
    include("restraints.jl")
    include("custom_ref_energy.jl")

    @info " | Loading Energy Function"
    include("energy_function.jl")

    # include("ref15.jl")

    # --- Show available Energy Function Components

    """
        get_available_energy_function_components([m::Module = ProtoSyn.Calculators])

    Returns all available [`EnergyFunctionComponent`](@ref) instances in the
    provided `Module` `m` (defaults to `ProtoSyn.Calculators`).
    Recursivelly searches any inner `Module`.

    # See also
    [`show_available_energy_function_components`](@ref)

    # Examples
    ```
    julia> ProtoSyn.Calculators.get_available_energy_function_components(ProtoSyn.Calculators)
    14-element Vector{Function}:
     get_default_custom_ref_energy (generic function with 1 method)
     get_default_coulomb (generic function with 1 method)
     (...)
     ```
    """
    function get_available_energy_function_components(m::Module)
        energy_function_components = Vector{Function}()
        
        # 1. Get all functions inside the module
        all_functions = [x for x in names(m, all=true) if x ∉ (:eval, :include, :proxy) && getproperty(m, x) isa Function && !occursin("#", string(x))]

        # 2. Filter for functions that return an EnergyFunctionComponent
        for _function in all_functions
            return_types = Base.return_types(getfield(m, _function))
            is_energy_function_component = false
            for return_type in return_types
                if return_type <: ProtoSyn.Calculators.EnergyFunctionComponent && return_type !== Union{}
                    is_energy_function_component = true
                    break
                end
            end

            # 3. Save function that returns EnergyFunctionComponent
            if is_energy_function_component
                push!(energy_function_components, getfield(m, _function))
            end
        end

        # 4. List loaded modules
        loaded_modules = [x for x in names(m, all=true) if x ∉ (:eval, :include, :proxy) && getproperty(m,x) isa Module && string(x) != split(string(m), ".")[end]]

        # 5. Recursively search in the inner modules
        for loaded_module in loaded_modules
            _m = getfield(m, loaded_module)
            nefc = get_available_energy_function_components(_m)
            energy_function_components = vcat(energy_function_components, nefc)
        end

        return energy_function_components
    end

    
    """
        show_available_energy_function_components([io::IO = stdout], [m::Module = ProtoSyn.Calculators])

    Prints all available [`EnergyFunctionComponent`](@ref) instances in the
    provided `Module` `m` (defaults to `ProtoSyn.Calculators`) to the
    given `IO` `io` (defaults to `stdout`). Recursivelly searches any inner
    `Module`.

    # See also
    [`get_available_energy_function_components`](@ref)

    # Examples
    ```jldoctest
    julia> ProtoSyn.Calculators.show_available_energy_function_components(ProtoSyn.Calculators)
    +------------------------------------------------------------------------------------------------+
    | Index | Component name            | Function                                                   |
    +------------------------------------------------------------------------------------------------+
    | 1     | Custom_Ref_Energy         | ProtoSyn.Calculators.get_default_custom_ref_energy         |
    | 2     | Coulomb                   | ProtoSyn.Calculators.Electrostatics.get_default_coulomb    |
    | 3     | GB_Solvation              | ProtoSyn.Calculators.GB.get_default_gb                     |
    +------------------------------------------------------------------------------------------------+
    └── Consider using the `?` menu to learn more about each EnergyFunctionComponent.
    ```
    """
    function show_available_energy_function_components(io::IO, m::Module)
        energy_function_components = get_available_energy_function_components(m)

        println(io, "+"*repeat("-", 120)*"+")
        @printf(io, "| %-5s | %-25s | %-82s |\n", "Index", "Component name", "Function")
        println(io, "+"*repeat("-", 120)*"+")
        for (index, component) in enumerate(energy_function_components)
            c = component()
            @printf(io, "| %-5d | %-25s | %-82s |\n", index, c.name, string(methods(component)[1].module)*"."*string(component))
        end
        println(io, "+"*repeat("-", 120)*"+")
        println("  └── Consider using the `?` menu to learn more about each EnergyFunctionComponent.\n")
    end

    show_available_energy_function_components(m::Module) = show_available_energy_function_components(stdout, m)
    show_available_energy_function_components()          = show_available_energy_function_components(stdout, ProtoSyn.Calculators)
    get_available_energy_function_components()           = get_available_energy_function_components(ProtoSyn.Calculators)
end