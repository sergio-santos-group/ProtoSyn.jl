module Units

# When necessary, defaults all calls to Float64
const defaultFloat = Float64
const defaultCleanCacheEvery = Int16(20)

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

end