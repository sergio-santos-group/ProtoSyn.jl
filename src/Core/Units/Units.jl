module Units

# When necessary, defaults all calls to Float64
const defaultFloat = Float64

export °, Å, kJ, mol, rad, nm, m, J, tonumber

const   ° = deg2rad(1)
const   Å = 1
const  kJ = 1
const mol = 1
const rad = 1

const nm = 10Å
const m = 1e10Å
const J = kJ/1000

tonumber(v::Number) = v
tonumber(v::String) = eval(Meta.parse(v))

end