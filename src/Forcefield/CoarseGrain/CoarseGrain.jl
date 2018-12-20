module CoarseGrain

using ..Common
using JSON
using LinearAlgebra

const default_aa_coef = Dict("Q"=>-3.5, "W"=>-0.9, "T"=>-0.7, "C"=> 2.5, "P"=>-1.6, "V"=>4.2, "L"=> 3.8, "M"=> 1.9, "N"=>-3.5, "H"=>-3.2,
                             "A"=> 1.8, "D"=>-3.5, "G"=>-0.4, "E"=>-3.5, "Y"=>-1.3, "I"=>4.5, "S"=>-0.8, "K"=>-3.9, "R"=>-4.5, "F"=> 2.8)

include("components.jl")
include("evaluators.jl")
include("loader.jl")

end