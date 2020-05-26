module Units

export °, Å, nm, m, tonumber

const ° = deg2rad(1)
const Å = 1
const nm = 10Å
const m = 1e10Å

tonumber(v::Number) = v
# tonumber(v::String) = (x=eval(Meta.parse(v)); println(v,"=",x); x)
tonumber(v::String) = eval(Meta.parse(v))

end