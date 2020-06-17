module Units

export °, Å, kJ, mol, rad, nm, m, J, tonumber

const ° = deg2rad(1)
const   Å = 1
const  kJ = 1
const mol = 1
const rad = 1

const nm = 10Å
const m = 1e10Å
const J = kJ/1000

tonumber(v::Number) = v
# tonumber(v::String) = (x=eval(Meta.parse(v)); println(v,"=",x); x)
tonumber(v::String) = eval(Meta.parse(v))
kJ/mol/nm^2
end