"""
    calc_torchani_internal_energy(M::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::AbstractSelection, [update_forces::Bool = false]; [static_ref_energy::T = 0.0], [use_ensemble::Bool = false], [model::Int = 3]) where {T <: AbstractFloat}
    
Calculate and return the [`Pose`](@ref) `pose` internal energy according to a
single TorchANI model neural network. The internal energy is the TorchANI energy
value for intra-residue interactions only. This should be used in conjuntion
with the [`calc_torchani_model`](@ref) and [`calc_torchani_ensemble`](@ref) in
design efforts (thus providing a correct ΔΔG of mutation). The employed model
can be defined using the `model` parameter (from model 1 to 8, default is 3).
Optionally, if the `use_ensemble` flag is set to true, the whole TorchANI
ensemble will be used to calculate the internal energy. The optional `A`
parameter defines the acceleration type used. If left undefined the default
`ProtoSyn.acceleration.active` mode will be used. By setting the
`update_forces` flag to `true` (`false` by default), this function will also
calculate and return the forces acting on each atom based on a single
TorchANI model neural network. If provided, an `AbstractSelection` `selection`
defines a sub-set of [`Residue`](@ref) instances to calculate the internal
energy on (if the provided `AbstractSelection` `selection` is not of
[`Residue`](@ref) type, it will be promoted using the
[`promote`](@ref ProtoSyn.promote) method). The parameter `static_ref_energy`
defines an energy value that is summed to the calculated energy. In order to
improve performance: it is often the case where a single [`Residue`](@ref)
instance is being design at a time. In such cases, by applying an
`AbstractSelection` `selection` and defining the `static_ref_energy` as the
internal energy of all non-mutated [`Residue`](@ref) instances, a faster
internal energy calculation can be accurately reproduced. 

# See also:
[`calc_torchani_model`](@ref) [`calc_torchani_ensemble`](@ref) [`get_default_torchani_internal_energy`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.TorchANI.calc_torchani_internal_energy(pose, nothing)
(-1.9374065455049276, nothing)
```
"""
function calc_torchani_internal_energy(M::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::AbstractSelection, update_forces::Bool = false; static_ref_energy::T = 0.0, use_ensemble::Bool = false, model::Int = 3) where {T <: AbstractFloat}

    sele = ProtoSyn.promote(selection, Residue)

    e = 0.0
    for residue in sele(pose, gather = true)
        if use_ensemble
            _e, _ = ProtoSyn.Calculators.TorchANI.calc_torchani_ensemble(M, pose, SerialSelection{Residue}(residue.id, :id), false)
        else
            _e, _ = ProtoSyn.Calculators.TorchANI.calc_torchani_model(M, pose, SerialSelection{Residue}(residue.id, :id), false, model = model)
        end
        e -= _e
    end

    return e + static_ref_energy, nothing
end

function calc_torchani_internal_energy(M::Type{<: ProtoSyn.AbstractAccelerationType}, pose::Pose, selection::Nothing, update_forces::Bool = false; static_ref_energy::T = 0.0, use_ensemble::Bool = false, model::Int = 3) where {T <: AbstractFloat}
    calc_torchani_internal_energy(M, pose, TrueSelection{Residue}(), update_forces, static_ref_energy = static_ref_energy, use_ensemble = use_ensemble, model = model)
end

calc_torchani_internal_energy(pose::Pose, selection::Opt{AbstractSelection}, update_forces::Bool = false; static_ref_energy::T = 0.0, use_ensemble::Bool = false, model::Int = 3) where {T <: AbstractFloat} = begin
    calc_torchani_internal_energy(ProtoSyn.acceleration.active, pose, selection, update_forces, static_ref_energy = static_ref_energy, use_ensemble = use_ensemble, model = model)
end


"""
    get_default_torchani_internal_energy(;[α::T = 1.0]) where {T <: AbstractFloat}

Return the default TorchANI internal energy model
[`EnergyFunctionComponent`](@ref). `α` sets the component weight (on an
[`EnergyFunction`](@ref ProtoSyn.Calculators.EnergyFunction) instance). This
component employs the [`calc_torchani_internal_energy`](@ref) method, therefore
predicting a structure's TorchANI internal energy based on a single model or
ensemble.

# Settings
* `use_ensemble::Bool` - Defines whether to use a single model or the whole TorchANI ensemble to calculate the structure's internal energy;
* `model::Int` - If using a single model to calculate the structure's internal energy, define which model (from 1 to 8);
* `static_ref_energy::T` - Define a static energy value to add to the calculated internal energy (where T <: AbstractFloat);

# See also
[`calc_torchani_model`](@ref) [`get_default_torchani_ensemble`](@ref)

# Examples
```jldoctest
julia> ProtoSyn.Calculators.TorchANI.get_default_torchani_internal_energy()
🞧  Energy Function Component:
+---------------------------------------------------+
| Name           | TorchANI_Ref_Energy              |
| Alpha (α)      | 1.0                              |
| Update forces  | true                             |
| Calculator     | calc_torchani_internal_energy    |
+---------------------------------------------------+
 |    +----------------------------------------------------------------------------------+
 ├──  ● Settings                      | Value                                            |
 |    +----------------------------------------------------------------------------------+
 |    | static_ref_energy             | 0.0                                              |
 |    | use_ensemble                  | false                                            |
 |    | model                         | 3                                                |
 |    +----------------------------------------------------------------------------------+
 |    
 └──  ○  Selection: nothing
```
"""
function get_default_torchani_internal_energy(;α::T = 1.0) where {T <: AbstractFloat}
    return EnergyFunctionComponent(
        "TorchANI_Ref_Energy",
        calc_torchani_internal_energy,
        nothing,
        Dict{Symbol, Any}(:static_ref_energy => 0.0, :use_ensemble => false, :model => 3),
        α,
        true)
end


"""
    fixate_static_ref_energy!(efc::EnergyFunctionComponent, pose::Pose, [selection::Opt{AbstractSelection} = nothing])

If the given [`EnergyFunctionComponent`](@ref) `efc` is a TorchANI Reference
Energy [`EnergyFunctionComponent`](@ref), calculate a new `static_ref_energy` on
the given [`Pose`](@ref) `pose` and apply it as a static internal energy value
(improved performance). If an `AbstractSelection` `selection` is provided, the
internal energy is only calculated in the subset of selected [`Residue`](@ref)
instances.

# See also
[`get_default_torchani_internal_energy`](@ref)

# Examples
```
julia> ProtoSyn.Calculators.TorchANI.fixate_static_ref_energy!(efc, pose, rid"1:20")
-0.3799925707280636
```
"""
function fixate_static_ref_energy!(efc::EnergyFunctionComponent, pose::Pose, selection::Opt{AbstractSelection} = nothing)
    if :static_ref_energy in keys(efc.settings)
        efc.settings[:static_ref_energy] = calc_torchani_internal_energy(pose, selection, use_ensemble = efc.settings[:use_ensemble], model = efc.settings[:model])[1]
    end
end