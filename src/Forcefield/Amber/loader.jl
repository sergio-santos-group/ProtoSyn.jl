@doc raw"""
    load_from_json(i_file::String)::Forcefield.Topology

Gather all topology components and return a [`Amber.Topology`](@ref Forcefield) object, parsing a JSON file.

# Examples
```julia-repl
julia> Forcefield.Amber.load_from_json(json_file)
Forcefield.Amber.Topology(
 atoms=ProtoSyn.Forcefield.Atom[Forcefield.Atom(name="N", σ=0.325, ϵ=0.711, q=0.0017, excls=[0, 1, 2, 3, 4, 5], pairs=[4, 5]), ...],
 bonds=ProtoSyn.Forcefield.HarmonicBond[Forcefield.HarmonicBond(a1=1, a2=2, k=2500.0, b0=0.19), ...],
 angles=ProtoSyn.Forcefield.HarmonicAngle[Forcefield.HarmonicAngle(a1=1, a2=2, a3=3, k=670.0, θ=1.92), ...],
 dihedralsCos=ProtoSyn.Forcefield.DihedralCos[Forcefield.DihedralCos(a1=1, a2=2, a3=3, a4=4, k=10.46, θ=180.0, mult=2.0), ...])
```
"""
function load_from_json(i_file::String)::Topology
    
    #Parse JSON file content to Dictionary
    ffield = Dict()
    open(i_file, "r") do f
        json_txt = read(f, String)
        ffield = JSON.parse(json_txt)
    end

    #Create atoms
    atoms = Vector{Atom}()
    for (index, atom) in enumerate(ffield["atoms"])
        index = index - 1
        name = ffield["atomtypes"][string(atom["type"])]["name"]
        q = atom["q"]
        σ = ffield["atomtypes"][string(atom["type"])]["sigma"]
        ε = ffield["atomtypes"][string(atom["type"])]["epsilon"]
        pre_excls = map(x -> x + 1, ffield["excls"][string(index)])
        excls = push!(pre_excls, -1)
        pairs = Vector{Int64}()
        for pair in ffield["pairs"]
            if pair["a1"] == index
                push!(pairs, (pair["a2"] + 1))
            end
        end

        #σ, ϵ and q are multiplied by the constants so that the geometric and arithmetic averages are
        #correct in the energy calculations (Inside the ATOM constructor)
        push!(atoms, Atom(name, σ, ε, q, excls, pairs))
    end

    # Fix bonds, angles and dihedrals numbering
    for component_type in ["bonds", "angles", "proper_dihedrals", "improper_dihedrals"]
        for component in ffield[component_type]
            for atom in ["a1", "a2", "a3", "a4"]
                if atom in keys(component)
                    component[atom] = component[atom] + 1
                end
            end
        end
    end

    #Create bonds
    bonds = Vector{HarmonicBond}()
    for bond in ffield["bonds"]
        cb = ffield["bondtypes"][string(bond["type"])]["cb"]
        b0 = ffield["bondtypes"][string(bond["type"])]["b0"]
        new_bond = HarmonicBond(bond["a1"], bond["a2"], cb, b0)
        push!(bonds, new_bond)
    end

    #Create angles
    angles = Vector{HarmonicAngle}()
    for angle in ffield["angles"]
        th = ffield["angletypes"][string(angle["type"])]["th"]
        ct = ffield["angletypes"][string(angle["type"])]["ct"]
        new_angle = HarmonicAngle(angle["a1"], angle["a2"], angle["a3"], ct, deg2rad(th))
        push!(angles, new_angle)
    end

    #Create dihedrals
    dihedralsCos = Vector{DihedralCos}()
    for dihedral in ffield["proper_dihedrals"]
        phi  = ffield["dihedraltypes"][string(dihedral["type"])]["phi"]
        cp   = ffield["dihedraltypes"][string(dihedral["type"])]["cp"]
        mult = ffield["dihedraltypes"][string(dihedral["type"])]["mult"]
        new_dihedralCos = DihedralCos(dihedral["a1"], dihedral["a2"], dihedral["a3"], dihedral["a4"], cp, deg2rad(phi), mult)
        push!(dihedralsCos, new_dihedralCos)
    end
    for dihedral in ffield["improper_dihedrals"]
        phi  = ffield["dihedraltypes"][string(dihedral["type"])]["phi"]
        cp   = ffield["dihedraltypes"][string(dihedral["type"])]["cp"]
        mult = ffield["dihedraltypes"][string(dihedral["type"])]["mult"]
        new_dihedralCos = DihedralCos(dihedral["a1"], dihedral["a2"], dihedral["a3"], dihedral["a4"], cp, deg2rad(phi), mult)
        push!(dihedralsCos, new_dihedralCos)
    end

    #Create topology
    return Topology(atoms, bonds, angles, dihedralsCos)
end