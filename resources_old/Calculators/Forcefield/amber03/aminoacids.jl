struct NameToTypeMapItem
    name::String
    type::String
    charge::Float64
end


let

_map = Dict()

ala = get!(_map, "ALA", Dict())

ala[:atoms] = Dict(
    "C"   => NameToTypeMapItem( 0.570224,  "C"),
    "CA"  => NameToTypeMapItem(-0.027733, "CT"),
    "CB"  => NameToTypeMapItem(-0.229951, "CT"),
    "H"   => NameToTypeMapItem( 0.294276,  "H"),
    "HA"  => NameToTypeMapItem( 0.120802, "H1"),
    "HB1" => NameToTypeMapItem( 0.077428, "HC"),
    "HB2" => NameToTypeMapItem( 0.077428, "HC"),
    "HB3" => NameToTypeMapItem( 0.077428, "HC"),
    "N"   => NameToTypeMapItem(-0.404773,  "N"),
    "O"   => NameToTypeMapItem(-0.555129,  "O"),
)
ala[:improper]

@ffdef ala :atoms NameToTypeMapItem begin
    (  "C",  "C",  0.570224),
    ( "CA", "CT", -0.027733),
    ( "CB", "CT", -0.229951),
    (  "H",  "H",  0.294276),
    ( "HA", "H1",  0.120802),
    ("HB1", "HC",  0.077428),
    ("HB2", "HC",  0.077428),
    ("HB3", "HC",  0.077428),
    (  "N",  "N", -0.404773),
    (  "O",  "O", -0.555129)
end

@ffdef ala :improper

ala[:atoms] = Dict(
    "C"   => Dict(:charge =>  0.570224, :type=>  "C"),
    "CA"  => Dict(:charge => -0.027733, :type=> "CT"),
    "CB"  => Dict(:charge => -0.229951, :type=> "CT"),
    "H"   => Dict(:charge =>  0.294276, :type=>  "H"),
    "HA"  => Dict(:charge =>  0.120802, :type=> "H1"),
    "HB1" => Dict(:charge =>  0.077428, :type=> "HC"),
    "HB2" => Dict(:charge =>  0.077428, :type=> "HC"),
    "HB3" => Dict(:charge =>  0.077428, :type=> "HC"),
    "N"   => Dict(:charge => -0.404773, :type=>  "N"),
    "O"   => Dict(:charge => -0.555129, :type=>  "O"),
)
ala[:improper] = (
    ()
)
end