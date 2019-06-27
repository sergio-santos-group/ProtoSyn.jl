include("structures.jl")

function get_atom_by_index(atoms::Vector{Atom}, index::Int64)::Opt{Atom}
    for atom in atoms
        if atom.index == index
            return atom
        end
    end
    return nothing
end


function read_pdb(input::String)::Tuple{Vector{Atom}, Dict{Int64, Vector{Int64}}} # Could be better -.-
    atoms = Vector{Atom}()
    conects = Dict{Int64, Vector{Int64}}()
    open(input, "r") do file
        pdb = readlines(file)
        map(pdb) do line
            m = match(r"ATOM\s*(\d*)\s*(\w*)\s*(\w*)", line)
            if m != nothing
                push!(atoms, Atom(m.captures))
            else
                m = match(r"CONECT\s*(\d*)\s*(\d*)\s*(\d*)\s*(\d*)", line)
                if m!= nothing
                    # atom = get_atom_by_index(atoms, parse(Int64, m.captures[1]))
                    local_conects = Vector{Int64}()
                    for index in m.captures[2:end]
                        # if length(index) > 0 && parse(Int64, index) > atom.index
                        if length(index) > 0
                            push!(local_conects, parse(Int64, index))
                        end
                    end
                    conects[parse(Int64, m.captures[1])] = local_conects
                end
            end
        end
    end
    return atoms, conects
end


function load_ff_rtp(input::String)::Dict{String, Dict{String, String}}
    name_2_atomtype = Dict{String, Dict{String, String}}()
    open(input*"/aminoacids.rtp", "r") do file
        current_aminoacid = ""
        pdb = readlines(file)
        map(pdb) do line
            m = match(r"\[\s(.{3})\s\]", line)
            if m != nothing
                current_aminoacid = String(m.captures[1])
                name_2_atomtype[current_aminoacid] = Dict{String, String}()
            else
                m = match(r"\s*(\w{1,3})\s*(\w*)\s*((?:-|\s)[0-9]\.[0-9]*)\s*\d*", line)
                if m != nothing
                    name_2_atomtype[current_aminoacid][String(m.captures[1])] = String(m.captures[2])
                end
            end
        end
    end
    return name_2_atomtype
end


function load_ff_non_bonded(input::String)::Forcefield

    parameter_name_2_type = Dict{Symbol, DataType}(:bondtypes => HarmonicBondType, :angletypes => HarmonicAngleType)

    forcefield = Forcefield()
    open(input*"/ffbonded.itp", "r") do file
        pdb = readlines(file)

        current_parameter = nothing
        current_collection = Vector{Any}()

        map(pdb) do line

            # Search for a new parameter on each line
            m = match(r"\[\s(\w*)\s\]", line)

            # If a new parameter has been found:
            if m != nothing

                # A parameter has already been filled and a new was found:
                # -> Save the current parameter
                if current_parameter != nothing && current_parameter in fieldnames(typeof(forcefield))
                    setproperty!(forcefield, current_parameter, current_collection)
                end

                # If a new parameter has been found:
                # -> Save this parameter as a Symbol
                # -> Start a new empty collection to save this parameter
                current_parameter = Symbol(m.captures[1])
                current_collection = Vector{Any}()

            # If a new parameter has not been found:
            else
                # Search for components belonging to the current parameter
                m = match(r"((?:\w+\**\s+)+)\d\s+((?:\d+.\d+\s+)+)(\d)*", line)

                # If a new component has been found:
                if m != nothing && current_parameter in keys(parameter_name_2_type)
                    results = map(split, filter(x -> x != nothing, m.captures))
                    results = collect(Iterators.flatten(results))
                    push!(current_collection, parameter_name_2_type[current_parameter](results))
                end
            end
        end
    end

    return forcefield
end


function group_iteration(atoms::Vector{Atom}, conects::Dict{Int64, Vector{Int64}}, group_every::Int64)::Vector{Vector{Atom}}
    
    function group(conects::Dict{Int64, Vector{Int64}}, group_every::Int64, start::Int64, original_result::Vector{Int64}, results::Vector{Vector{Int64}})::Vector{Vector{Int64}}

        for conect in conects[start]
            current_result = copy(original_result)
            if conect in current_result
                continue
            end
            push!(current_result, conect)
            l = length(current_result)
            if l == group_every && current_result[1] < current_result[l]
                push!(results, current_result)
            else
                results = group(conects, group_every, conect, current_result, results)
            end
        end

        return results
    end

    results = Vector{Vector{Atom}}()
    for atom_index in keys(conects)
        new = group(conects, group_every, atom_index, Int64[atom_index], Vector{Vector{Int64}}())
        results = vcat(results, map(x -> map(y -> get_atom_by_index(atoms, y), x), new))
    end
    return results
end


function compile_atomtypes!(atoms::Vector{Atom}, c::Dict{String, Dict{String, String}})
    for atom in atoms
        atom.atomtype = c[atom.residue][atom.name]
    end
end

function compile_topology(atoms::Vector{Atom}, conects::Dict{Int64, Vector{Int64}}, forcefield::Forcefield)::Topology

    parameter_name_2_type = Dict{Symbol, DataType}(:bondtypes => HarmonicBond, :angletypes => HarmonicAngle)

    function get_parameters_from_atomtypes(
        element::Vector{Atom},
        component_list::Vector{<:AbstractForcefieldComponent})

        element_atomtypes = map(x -> x.atomtype, element)

        for component in component_list
            if all(element_atomtypes == component.ordered_atomtypes) || all(element_atomtypes == component.reverse_ordered_atomtypes)
                return component
            end
        end
    end

    topology = Topology()

    for component in fieldnames(typeof(forcefield))

        current_collection = Vector{Any}()

        component_list = getproperty(forcefield, component)

        atom_count = component_list[1].atom_count
        for element in group_iteration(atoms, conects, atom_count)
            component_info = Any[x for x in element]
            parameters = get_parameters_from_atomtypes(element, component_list)
            for i in 4:length(fieldnames(typeof(parameters)))
                push!(component_info, getproperty(parameters, fieldnames(typeof(parameters))[i]))
            end
            push!(current_collection, parameter_name_2_type[component](component_info...))
        end

        topology_component = Symbol(replace(String(component), "type" => ""))
        setproperty!(topology, topology_component, current_collection)
    end

    return topology
end




# ----
function main()
    atoms, conects = read_pdb("sample.pdb")
    name_2_atomtype = load_ff_rtp("amber99sb-ildn.ff")
    ff = load_ff_non_bonded("amber99sb-ildn.ff")
    compile_atomtypes!(atoms, name_2_atomtype)
    topol = compile_topology(atoms, conects, ff)
    
    for entry in topol.angles
        println(entry)
    end
end

@time main()
@time main()