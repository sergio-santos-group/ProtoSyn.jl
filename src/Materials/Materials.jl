@doc """
    Materials

The Materials modules introduces several methods related to different types of
materials (such as atomic lattices, etc).
"""
module Materials

    include("Methods/lattices.jl")
    include("Methods/perlin.jl")
    include("Methods/carbon.jl") # Requieres perlin.jl
end