# module Common

# using ..Aux
# using LinearAlgebra

# include("energy.jl")
# include("residue.jl")
# include("metadata.jl")
# include("dihedral.jl")
# include("ss_blocks.jl")
# include("state.jl")
# include("callback.jl")
# include("loaders.jl")
# include("conf_gen.jl")
# include("macros.jl")

# const ss2bbd = Dict(
#         Common.SS.SHEET => Dict(Common.DIHEDRAL.phi => deg2rad(-115.0), Common.DIHEDRAL.psi => deg2rad(135.0)),
#         Common.SS.HELIX => Dict(Common.DIHEDRAL.phi => deg2rad(-60.0),  Common.DIHEDRAL.psi => deg2rad(-45.0)))
# end

module Common

macro customshow(type::Symbol, args::Symbol...)
quote
    function Base.show(io::IO, b::$(type))
        print(io, string($(type)), " {")
        for p in setdiff(fieldnames($(type)), [:children, :parent])
            v = getproperty(b,p)
            # if v isa AbstractContainerMetadata
            #     # s = "<container with $(length(v.children)) items>"
            #     s = "<$(length(v.children)) items of type $(eltype(v.children))>"
            # else
            if v == nothing
                s = "n.d."
            else
                s = string(v)
            end
            print(io, "\n   $(string(p)) = $s")
        end
        if b isa AbstractMetadata
            if :children in fieldnames($(type))
                v = b.children
                s = isa(v,Nothing) ? "n.d." : "<vector with $(length(v)) $(eltype(v))>"
                print(io, "\n   children = $s")
            end
            if :parent in fieldnames($(type))
                v = b.parent
                s = isa(v,Nothing) ? "n.d." : "@$(v.name).$(v.index)"
                print(io, "\n   parent = $s")
            end
        end
        print(io, "\n}")
    end
end
end


include("types.jl")

#include("methods.jl")
include("test.jl")

end


# @enum  DType begin
#  ϕ = 1
#  ψ = 2
# end

# using .Common
# d = Common.DihedralMetadata(1,2,3,4,ϕ)
# println(d)

