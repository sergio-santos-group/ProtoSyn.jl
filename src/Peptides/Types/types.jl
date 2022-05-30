using Printf
using ProtoSyn
using ProtoSyn.Units

"""
    Dihedral

An auxiliary struct for [`Dihedral`](@ref) selection
by name. Each of the available dihedrals (under the abstract type
`DihedralType`, see bellow) can be called using the following signature,
returning the respective representative atom:

```julia
(dihedral::DihedralType)(residue::Residue)
```

# Available dihedrals

`Dihedral.phi` `Dihedral.psi` `Dihedral.omega`
`Dihedral.chi1` `Dihedral.chi2` `Dihedral.chi3` `Dihedral.chi4`

# Examples
```
julia> ProtoSyn.Peptides.Dihedral.phi(pose.graph[1][1])
Atom{/UNK:1/UNK:1/SER:1/C:10}

julia> atom = ProtoSyn.Peptides.Dihedral.chi2(pose.graph[1][2])
Atom{/UNK:1/UNK:1/GLU:2/CD:22}

julia> ProtoSyn.getdihedral(pose.state, atom)
-3.141165
```

"""
module Dihedral

    using ProtoSyn

    abstract type DihedralType end

    # --- Backbone dihedral angles

    struct Phi <: DihedralType end; const phi = Phi()
    struct Psi <: DihedralType end; const psi = Psi()
    struct Omega <: DihedralType end; const omega = Omega()

    (phi::Phi)(residue::Residue)     = begin
        if residue.parent == ProtoSyn.root(residue)
            @warn "Residue $residue has no phi angle"
            return nothing
        end
        return residue["C"]
    end

    (psi::Psi)(residue::Residue)     = begin
        if length(residue.children) == 0 
            @warn "Residue $residue has no psi angle"
            return nothing
        end
        return residue.children[1]["N"]
    end

    (omega::Omega)(residue::Residue) = begin
        if residue.parent == ProtoSyn.root(residue)
            @warn "Residue $residue has no omega angle"
            return nothing
        end
        return residue["CA"]
    end

    # --- Chi dihedral angles

    struct Chi1 <: DihedralType end; const chi1 = Chi1()
    struct Chi2 <: DihedralType end; const chi2 = Chi2()
    struct Chi3 <: DihedralType end; const chi3 = Chi3()
    struct Chi4 <: DihedralType end; const chi4 = Chi4()

    function get_chi(residue::Residue, chi::Int)
        available_chis = chi_dict[residue.name.content]
        @assert chi < length(available_chis) "Tried to retrieve chi $chi on a residue with only $(length(available_chis) - 1) chi angles defined."
        chi_atom = available_chis[chi + 1] # Note the "+ 1"
        if residue[chi_atom] === nothing
            @warn "Chi $chi of residue $residue requires atom $(chi_atom), which was not found"
            return nothing
        end
        
        return residue[chi_atom]
    end

    (chi1::Chi1)(residue::Residue) = get_chi(residue, 1)
    (chi2::Chi2)(residue::Residue) = get_chi(residue, 2)
    (chi3::Chi3)(residue::Residue) = get_chi(residue, 3)
    (chi4::Chi4)(residue::Residue) = get_chi(residue, 4)

    # Dict built according to Dunbrack rotamer library.
    # Dihedral angles not in the dict can still be rotated, by manually defining
    # the controling atom. For example, ALA-CB dihedral can be rotated by
    # setting the dihedral value Δϕ on any of its children. By definition, 
    # dihedral angles which only rotate the position of hydrogens are not
    # considered in this dict.
    chi_dict = Dict{String, Vector{String}}(
        # "ALA" => ["CB"],
        # "CYS" => ["CB", "SG"],
        # "ASP" => ["CB", "CG", "OD1"],
        # "GLU" => ["CB", "CG", "CD", "OE1"],
        # "PHE" => ["CB", "CG", "CD1"],
        # "GLY" => [],
        # "HIS" => ["CB", "CG", "ND1"],
        # "HIE" => ["CB", "CG", "ND1"],
        # "ILE" => ["CB", "CG1", "CD1", "HD11"],
        # "LYS" => ["CB", "CG", "CD", "CE", "NZ"],
        # "LEU" => ["CB", "CG", "CD1"],
        # "MET" => ["CB", "CG", "SD", "CE"],
        # "ASN" => ["CB", "CG", "OD1", "HD21"],
        # "PRO" => ["CB", "CG", "CD", "HD2"], # Attention
        # "GLN" => ["CB", "CG", "CD", "NE2"],
        # "ARG" => ["CB", "CG", "CD", "NE", "CZ"],
        # "SER" => ["CB", "OG"],
        # "THR" => ["CB", "OG1"],
        # "VAL" => ["CB", "CG1"],
        # "TRP" => ["CB", "CG", "CD1"],
        # "TYR" => ["CB", "CG", "CD1"]
    )
end

include("secondary-structure.jl")