module Units

# When necessary, defaults all calls to Float64
const defaultFloat = Float64
const defaultCleanCacheEvery = Int16(10000)
const max_gpu_allocation = defaultFloat(0.5)
const defaultTorchANIport = 50000

export °, Å, kJ, mol, rad, nm, m, J, tonumber

const   ° = deg2rad(1)
const   Å = 1
const  kJ = 1
const mol = 1
const rad = 1

const nm = 10Å
const m = 1e10Å
const J = kJ/1000

export tonumber

tonumber(T::DataType, v::Number) = T(v)
tonumber(v::Number)      = tonumber(defaultFloat, v)

tonumber(T::DataType, v::String) = T(eval(Meta.parse(v)))
tonumber(v::String)      = tonumber(defaultFloat, v)

# * CONSTANTS

const max_bond_lengths = Dict{String, defaultFloat}(
    "CC" => 1.6,
    "CH" => 1.2,
    "CN" => 2.2,
    "CO" => 2.2,
    "CS" => 2.6,
    "NH" => 1.2,
    "OH" => 1.2,
    "HC" => 1.2,
    "NC" => 2.2,
    "OC" => 2.2,
    "SC" => 2.6,
    "HN" => 1.2,
    "HO" => 1.2)

end