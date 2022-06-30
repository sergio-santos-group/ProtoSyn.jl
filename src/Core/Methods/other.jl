using Statistics
using YAML

"""
    rand_vector_in_sphere([::Type{T}]) where {T <: AbstractFloat}

Return a random vector in a sphere, with norm = 1. The X, Y and Z cartesian
coordinates are of type `T`. If not provided, the `ProtoSyn.Units.defaultFloat`
type will be used.

# Examples
```
julia> ProtoSyn.rand_vector_in_sphere()
3-element Array{Float64,1}:
 0.4182585394795512
 0.4965570141814753
 0.7605859095589471
```
"""
function rand_vector_in_sphere(::Type{T}) where {T <: AbstractFloat}
    theta = 2 * π * rand(T)
    phi   = acos(1 - 2 * rand(T))
    x     = sin(phi) * cos(theta)
    y     = sin(phi) * sin(theta)
    z     = cos(phi)
    return Vector{T}([x, y, z])
end

rand_vector_in_sphere() = rand_vector_in_sphere(ProtoSyn.Units.defaultFloat)


"""
    fibonacci_sphere([::Type{T}], n_points::Int = 1000) where {T <: AbstractFloat}

Generate a set `n_points` cartesian coordinates following the fibonnacci lattice
(see http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/).
If provided, use the `AbstractFloat` type `T` (otherwise, will use
`ProtoSyn.Units.defaultFloat`).

# Examples
```jldoctest
julia> ProtoSyn.fibonacci_sphere(10)
10-element Vector{Vector{Float64}}:
 [0.0, 1.0, 0.0]
 [-0.4634653634889746, 0.7777777777777778, 0.42457223795379545]
 [0.07269269081806179, 0.5555555555555556, -0.8282957185649263]
 (...)
```
"""
function fibonacci_sphere(::Type{T}, n_points::Int = 1000) where {T <: AbstractFloat}
    points = Vector{Vector{T}}()
    ϕ = pi * (3.0 - sqrt(5.0))

    for i in 0:(n_points - 1)
        y = 1 - (i / (n_points - 1)) * 2
        radius = sqrt(1 - y * y)
        θ = ϕ * i
        x = cos(θ) * radius
        z = sin(θ) * radius

        push!(points, [x, y, z])
    end

    return points
end

fibonacci_sphere(n_points::Int = 1000) = fibonacci_sphere(ProtoSyn.Units.defaultFloat, n_points)


"""
    rotation_matrix_from_axis_angle(axis::Vector{T}, angle::T) where {T <: AbstractFloat}

Return a rotation matrix based on the provided axis and angle (in radians).

# Examples
```jldoctest
julia> ProtoSyn.rotation_matrix_from_axis_angle([1.1, 2.2, 3.3], π/2)
3×3 Matrix{Float64}:
  0.0714286  -0.658927  0.748808
  0.944641    0.285714  0.16131
 -0.320237    0.695833  0.642857
```
"""
function rotation_matrix_from_axis_angle(axis::Vector{T}, angle::T) where {T <: AbstractFloat}
    q0 = cos(0.5 * angle)
    q1, q2, q3 = sin(0.5 * angle) * axis ./ norm(axis)
    return [1-2*q2*q2-2*q3*q3   2*q1*q2-2*q0*q3   2*q1*q3+2*q0*q2;
            2*q2*q1+2*q0*q3   1-2*q3*q3-2*q1*q1   2*q2*q3-2*q0*q1;
            2*q3*q1-2*q0*q2     2*q3*q2+2*q0*q1 1-2*q1*q1-2*q2*q2]
end


"""
    unit_circle(value::T) where {T <: AbstractFloat}

Maps an angle value (in radians) from ]-∞, +∞[ to [0, 2π[ range. Returns value
in radians.

# Examples
```jldoctest
julia> ProtoSyn.unit_circle(-2pi)
0.0

julia> ProtoSyn.unit_circle(-(3/2)pi)
1.5707963267948966
```
"""
function unit_circle(value::T) where {T <: AbstractFloat}
    r = deg2rad((value > 0 ? value : (2*pi + value)) * 360 / (2*pi))
    return mod(r, 2*pi)
end


"""
    half_unit_circle(value::T) where {T <: AbstractFloat}

Maps an angle value (in radians) from ]-∞, +∞[ to [-π, π[ range. Returns value
in radians.

# Examples
```jldoctest
julia> ProtoSyn.half_unit_circle(-2pi)
0.0

julia> ProtoSyn.half_unit_circle(-(3/2)pi)
-1.5707963267948966
```
"""
function half_unit_circle(value::T) where {T <: AbstractFloat}
    return value % pi
end


"""
    gpu_allocation()

Return the current fraction of the GPU memory allocated (in range [0, 1]).

# Examples
```
julia> ProtoSyn.gpu_allocation()
0.07821749405343822
```
"""
function gpu_allocation()
    mem = (CUDA.total_memory() - CUDA.available_memory()) / CUDA.total_memory()
    return ProtoSyn.Units.defaultFloat(mem)
end


"""
    read_yml(_filename::String)

Open and read the contents of a .yml file.

# Examples
```
julia> ProtoSyn.read_yml("resources/Peptides/grammars.yml")
Dict{Any, Any} with 2 entries:
 (...)
```
"""
function read_yml(_filename::String)
    open(_filename, "r") do io
        return YAML.load(io)
    end
end

# --- Tree display -------------------------------------------------------------

export LevelCode
mutable struct LevelCode
    levels::Vector{Int}
    code_table::Dict{Int, String}
    conv_table::Dict{Int, Int}

    LevelCode(levels::Vector{Int}, code_table::Dict{Int, String}, conv_table::Dict{Int, Int}) = begin

        code_table_keys = keys(code_table)
        for value in levels
            if !(value in code_table_keys)
                error("Tried to create LevelCode instance with level type $value, which is not found in the current code table. Possible values: $code_table_keys")
            end
        end
        for (key, value) in conv_table
            if !(key in code_table_keys) | !(value in code_table_keys)
                error("Tried to create LevelCode instance with conv_table $key => $value, which is not found in the current code table. Possible values: $code_table_keys")
            end
        end
        new(levels, code_table, conv_table)
    end
end

const default_code_table = Dict{Int, String}(
    1 => " |   ",
    2 => "     ",
    3 => " ├── ",
    4 => " └── ")

const default_conv_table = Dict{Int, Int}(
    3 => 1,
    4 => 2,
    1 => 1,
    2 => 2
)

LevelCode() = LevelCode(Vector{Int}(), default_code_table, default_conv_table)
LevelCode(levels::Vector{Int}) = LevelCode(levels, default_code_table, default_conv_table)

Base.vcat(level_code::LevelCode, value::Int) = begin
    code_table_keys = sort(collect(keys(level_code.code_table)))
    if !(value in code_table_keys)
        error("Tried to append level type $value to a LevelCode instance, but it was not found in the current code table. Possible values: $code_table_keys")
    end
    return LevelCode(vcat(level_code.levels, value), level_code.code_table, level_code.conv_table)
end

Base.copy(level_code::LevelCode) = begin
    return LevelCode(
        copy(level_code.levels),
        copy(level_code.code_table),
        copy(level_code.conv_table))
end

Base.getindex(level_code::LevelCode, i::Int) = level_code.levels[i]
Base.setindex!(level_code::LevelCode, val::Int, i::Int) = begin
    code_table_keys = sort(collect(keys(level_code.code_table)))
    if !(val in code_table_keys)
        error("Tried to append level type $val to a LevelCode instance, but it was not found in the current code table. Possible values: $code_table_keys")
    end
    level_code.levels[i] = val
end
Base.lastindex(level_code::LevelCode) = length(level_code.levels)


"""
    get_lead(level_code::Opt{LevelCode} = nothing)

Auxiliary function, uses level codes to show inheritance in the terminal for
complex types.
"""
function get_lead(level_code::Opt{LevelCode} = nothing)
    lead = ""
    level_code === nothing && return lead
    length(level_code.levels) === 0 && return lead

    for level in level_code.levels[1:(end-1)]
        lead *= level_code.code_table[level_code.conv_table[level]]
    end

    return lead * level_code.code_table[level_code.levels[end]]
end

function get_inner_lead(level_code::Opt{LevelCode} = nothing)
    lead = ""
    level_code === nothing && return lead
    length(level_code.levels) === 0 && return lead

    inner_level_code = copy(level_code)
    inner_level_code[end] = level_code.conv_table[level_code.levels[end]]
    return ProtoSyn.get_lead(inner_level_code)
end