# Taken from https://github.com/ararslan/ProceduralNoise.jl
# All credit goes to ararslan
# Based on https://gist.github.com/Flafla2/f0260a861be0ebdeef76
# Hash lookup table as defined by Ken Perlin
# Minor changes and updates by José Pereira

const PERMS1 = [151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7,
    225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247,
    120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
    88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134,
    139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220,
    105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80,
    73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86,
    164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38,
    147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189,
    28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101,
    155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232,
    178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12,
    191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181,
    199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236,
    205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180]

const PERMS = [PERMS1; PERMS1]


function fade(x::T) where {T<:AbstractFloat}
    return 6x^5 - 15x^4 + 10x^3
end


function smoothstep(a0::T, a1::T, x::T) where {T<:AbstractFloat}
    a0 != a1 || throw(ArgumentError("Arguments a0 and a1 cannot be equal"))
    x = clamp((x - a0) / (a1 - a0), 0.0, 1.0)
    return 3x^2 - 2x^3
end


function smootherstep(a0::T, a1::T, x::T) where {T<:AbstractFloat}
    a0 != a1 || throw(ArgumentError("Arguments a0 and a1 cannot be equal"))
    x = clamp((x - a0) / (a1 - a0), 0.0, 1.0)
    return fade(x)
end


function interpolate(a0::T, a1::T, w::T, method::Symbol = :linear) where {T<:AbstractFloat}
    w = round(w, digits = 10)
    0 ≤ w ≤ 1 || throw(ArgumentError("Expected 0 ≤ w ≤ 1, got $w"))
    if method ∉ [:linear, :smoothstep, :smootherstep]
        throw(ArgumentError("Unrecognized interpolation method $method"))
    end
    x = method == :linear ? w : eval(method)(a0, a1, w)
    return a0 + x * (a1 - a0)
end


function gradient(hash::Int, x::T, y::T, z::T) where {T<:AbstractFloat}
    h = hash & 15
    u = h < 8 ? x : y
    v = h < 4 ? y : h == 12 || h == 14 ? x : z
    return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v)
end


"""
    perlin(x::T, y::T, z::T)

Compute the perlin noise value at coordinates `x`, `y` and `z`. Taken from
[https://github.com/ararslan/ProceduralNoise.jl](https://github.com/ararslan/ProceduralNoise.jl).
All credit goes to ararslan. Based on [https://gist.github.com/Flafla2/f0260a861be0ebdeef76](https://gist.github.com/Flafla2/f0260a861be0ebdeef76).
Hash lookup table as defined by Ken Perlin.

# Examples
```
julia> ProtoSyn.Materials.perlin(0.5, 0.5, 0.5)
0.375
```
"""
function perlin(x::T, y::T, z::T) where {T<:AbstractFloat}
    xi = trunc(Int, x) & 255 + 1
    yi = trunc(Int, y) & 255 + 1
    zi = trunc(Int, z) & 255 + 1

    xf = abs(first(modf(x)))
    yf = abs(first(modf(y)))
    zf = abs(first(modf(z)))

    u = fade(xf)
    v = fade(yf)
    w = fade(zf)

    aaa = PERMS[PERMS[PERMS[xi] + yi] + zi]
    aba = PERMS[PERMS[PERMS[xi] + yi + 1] + zi]
    aab = PERMS[PERMS[PERMS[xi] + yi] + zi + 1]
    abb = PERMS[PERMS[PERMS[xi] + yi + 1] + zi + 1]
    baa = PERMS[PERMS[PERMS[xi + 1] + yi] + zi]
    bba = PERMS[PERMS[PERMS[xi + 1] + yi + 1] + zi]
    bab = PERMS[PERMS[PERMS[xi + 1] + yi] + zi + 1]
    bbb = PERMS[PERMS[PERMS[xi + 1] + yi + 1] + zi + 1]

    x1 = interpolate(gradient(aaa, xf, yf, zf), gradient(baa, xf - 1, yf, zf), u)
    x2 = interpolate(gradient(aba, xf, yf - 1, zf), gradient(bba, xf - 1, yf - 1, zf), u)
    y1 = interpolate(x1, x2, v)

    x1 = interpolate(gradient(aab, xf, yf, zf - 1), gradient(bab, xf - 1, yf, zf - 1), u)
    x2 = interpolate(gradient(abb, xf, yf - 1, zf - 1), gradient(bbb, xf - 1, yf - 1, zf - 1), u)
    y2 = interpolate(x1, x2, v)

    return (interpolate(y1, y2, w) + 1) / 2
end


"""
    perlinoctaves(x::T, y::T, z::T, octaves::Int = 1, persistence::T = 1.0)

Compute the perlin noise value at coordinates `x`, `y` and `z` with the given
octaves. In sum, generates a perlin noise with some randomization. Taken from
[https://github.com/ararslan/ProceduralNoise.jl](https://github.com/ararslan/ProceduralNoise.jl).
All credit goes to ararslan. Based on [https://gist.github.com/Flafla2/f0260a861be0ebdeef76](https://gist.github.com/Flafla2/f0260a861be0ebdeef76).
Hash lookup table as defined by Ken Perlin.

# Examples
```
julia> ProtoSyn.Materials.perlinoctaves(0.5, 0.5, 0.5, 8)
0.484375
```
"""
function perlinoctaves(x::T, y::T, z::T, octaves::Int = 1, persistence::T = 1.0) where {T <: AbstractFloat}
    total     = 0.0
    frequency = 1.0
    amplitude = 1.0
    maxval    = 0.0
    for i in 1:octaves
        total += perlin(x * frequency, y * frequency, z * frequency) * amplitude
        maxval += amplitude
        amplitude *= persistence
        frequency *= 2
    end
    return total / maxval
end