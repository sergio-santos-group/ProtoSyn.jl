# Atom


```plantuml

class Atom {
    id : Int = -1
    name : String = "?"
    symbol : String = "?"
    x : Float64 = -1.0
    y : Float64 = -1.0
    z : Float64 = -1.0
    parent : Opt{T} = nothing
}

class Residue {
    name : String
    bonds : ConnectGraph? = nothing
    atoms : Opt{Vector{_Atom{Residue}}} = nothing
    atomsbyname : Opt{Dict{String, _Atom}} = nothing
}

Atom "0..*" o-- "1" Residue
@enduml
```
