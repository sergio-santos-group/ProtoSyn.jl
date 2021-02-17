# using LightXML

res_lib  = ProtoSyn.Peptides.grammar()
sequence = seq"MGSWA"
pose     = ProtoSyn.Peptides.build(res_lib, sequence)
ProtoSyn.Peptides.remove_sidechains!(pose)


using ProtoSyn
pose = ProtoSyn.Peptides.load("../../2a3d.pdb")
@time e2, f2 = ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose, update_forces = false);

@time for i in 1:10000
    e2, f2 = ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose, update_forces = false);
    if i%50 == 0
        println("Step $i/10_000")
        GC.gc(false)
    end
end

using Profile

using PyCall
using CUDA
@pyimport torch
@pyimport torchani

device = torch.device("cuda")
model = torchani.models.ANI2x(periodic_table_index = true).to(device)
model(...)

function teste()
    coordinates = torch.tensor([pose.state.x.coords'], requires_grad = true, device = device).float()
    s           = ProtoSyn.Calculators.TorchANI.get_ani_species(pose)
    species     = torch.tensor([s], device = device)
    m1 = model.species_converter((species, coordinates))
    m2 = model.aev_computer(m1)
    torch.cuda.empty_cache()
end

function GPU_memory()
    return CUDA.total_memory() - CUDA.available_memory() / CUDA.total_memory()
end

Profile.clear_malloc_data()
teste()
teste()
teste()
teste()
exit(0)

c = collect(pose.state.x.coords')
s = ProtoSyn.Calculators.TorchANI.get_ani_species(pose)

proxy = ProtoSyn.XMLRPC.ServerProxy("http://localhost", 50000)

"""
    Recursively travel the XMLDocument or XMLElement and gather all entrys of
    label type `query`, pushing them to the `results` vector (a parse for the
    correct T type is attempted).
"""
function r_xml_travel!(xml::Union{XMLDocument, XMLElement}, query::String, results::Vector{T}) where {T <: AbstractFloat}

    if typeof(xml) == XMLDocument
        xml = LightXML.root(xml)
    end
    for element in child_elements(xml)
        if name(element) == query
            value = string(collect(child_nodes(element))[1])
            push!(results, parse(T, value))
        end
        is_elementnode(element) && r_xml_travel!(element, query, results)
    end
end

function calc(s::Vector{Int64}, c::Matrix{Float64}, update_forces::Bool = false)
    response = proxy.calc(s, c, update_forces)
    xml_string = replace(join(Char.(response.body)), "\n" => "")
    xml = LightXML.parse_string(xml_string)
    
    r = Vector{Float64}()
    r_xml_travel!(xml, "double", r)
    return splice!(r, 1), reshape(r, 3, :)
end

@time e, f = calc(s, c, true);

@time for i in 1:100
    e, f = calc(s, c, true);
end

@time e, f = calc(s, c, false);

@time for i in 1:100
    e, f = calc(s, c, false);
end

@time e2, f2 = ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose, update_forces = true);

@time for i in 1:100
    e2, f2 = ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose, update_forces = true);
end

@time e2, f2 = ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose, update_forces = false);

@time for i in 1:100
    e2, f2 = ProtoSyn.Calculators.TorchANI.calc_torchani_model(pose, update_forces = false);
end