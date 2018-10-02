function load_from_json(json_file::String)
    
    #Parse JSON file content to Dictionary
    ffield = Dict()
    open(json_file, "r") do f
        json_txt = read(f, String)
        ffield = JSON.parse(json_txt)
    end

    #Create atoms
    atoms = Array{Forcefield.Atom, 1}()
    for (index, atom) in enumerate(ffield["atoms"])
        index = index - 1
        name = ffield["atomtypes"][string(atom["type"])]["name"]
        q = atom["q"]
        σ = ffield["atomtypes"][string(atom["type"])]["sigma"]
        ε = ffield["atomtypes"][string(atom["type"])]["epsilon"]
        pre_excls = map(x -> x + 1, ffield["excls"][string(index)])
        excls = push!(pre_excls, -1)
        pairs = Array{Int64, 1}()
        for pair in ffield["pairs"]
            if pair["a1"] == index
                push!(pairs, (pair["a2"] + 1))
            end
        end
        # fSI = 1/(4.pi.Eo) = 138.935485 kJ mol-1 nm e-2
        # sqrtfSI = sqrt(fSI);
        sqrtfSI=11.787089759563214

        #σ, ϵ and q are multiplied by the constants so that the geometric and arithmetic averages are
        #correct in the energy calculations
        push!(atoms, Forcefield.Atom(name, 0.5*σ, sqrt(ε), sqrtfSI*q, excls, pairs))
    end

    #Fix bonds, angles and dihedrals numbering
    for component_type in ["bonds", "angles", "dihedrals"]
        for component in ffield[component_type]
            for atom in ["a1", "a2", "a3", "a4"]
                if atom in keys(component)
                    component[atom] = component[atom] + 1
                end
            end
        end
    end

    #Create bonds
    bonds = Array{Forcefield.HarmonicBond, 1}()
    for bond in ffield["bonds"]
        cb = ffield["bondtypes"][string(bond["type"])]["cb"]
        b0 = ffield["bondtypes"][string(bond["type"])]["b0"]
        new_bond = Forcefield.HarmonicBond(bond["a1"], bond["a2"], cb, b0)
        push!(bonds, new_bond)
    end

    #Create angles
    angles = Array{Forcefield.Angle, 1}()
    for angle in ffield["angles"]
        th = ffield["angletypes"][string(angle["type"])]["th"]
        ct = ffield["angletypes"][string(angle["type"])]["ct"]
        new_angle = Forcefield.Angle(angle["a1"], angle["a2"], angle["a3"], ct, deg2rad(th))
        push!(angles, new_angle)
    end

    #Create dihedrals
    dihedralsCos = Array{Forcefield.DihedralCos, 1}()
    for dihedral in ffield["dihedrals"]
        phi  = ffield["dihedraltypes"][string(dihedral["type"])]["phi"]
        cp   = ffield["dihedraltypes"][string(dihedral["type"])]["cp"]
        mult = ffield["dihedraltypes"][string(dihedral["type"])]["mult"]
        new_dihedralCos = Forcefield.DihedralCos(dihedral["a1"], dihedral["a2"], dihedral["a3"], dihedral["a4"], cp, deg2rad(phi), mult)
        push!(dihedralsCos, new_dihedralCos)
    end

    #Create topology
    return Forcefield.Topology(atoms, bonds, angles, dihedralsCos)
end