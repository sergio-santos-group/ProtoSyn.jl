module Units

export °, Å, nm, tonumber

const ° = deg2rad(1)
const Å = 1
const nm = 10Å

tonumber(v::Number) = v
tonumber(v::String) = eval(Meta.parse(v))

end