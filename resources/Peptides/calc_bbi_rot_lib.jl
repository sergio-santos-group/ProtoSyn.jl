using ProtoSyn
using ProtoSyn.Peptides
using Base.Cartesian
using Images
using Printf

ncaa_residue_code = 'a'
N                 = 36
filename          = "ncaa_rotamers.lib"
energy_function = ProtoSyn.Calculators.EnergyFunction([
    ProtoSyn.Calculators.TorchANI.get_default_torchani_model()
])

res_lib    = ProtoSyn.Peptides.load_grammar_from_file("grammars.yml", "ncaa")
pose       = Peptides.build(res_lib, [string(ncaa_residue_code)])
step       = deg2rad(360)/N
n_chis     = length(Peptides.Dihedral.chi_dict[pose.graph[1, 1].name]) - 1
angles     = collect(Iterators.product(Tuple(repeat([0.0:step:2π], n_chis))...))
dimensions = length(size(angles)) === 1 ? (size(angles)[1], 1) : size(angles)
energies   = Matrix{Float64}(undef, dimensions)
rotamers   = Matrix{Peptides.Rotamer}(undef, dimensions)

_chis      = Vector{ProtoSyn.Peptides.Dihedral.DihedralType}([
    ProtoSyn.Peptides.Dihedral.chi1,
    ProtoSyn.Peptides.Dihedral.chi2,
    ProtoSyn.Peptides.Dihedral.chi3,
    ProtoSyn.Peptides.Dihedral.chi4
])
    
println("Apllying $N steps of $(rad2deg(step))° each.")
l = quote
    @nloops $n_chis i angles begin
        local a = @nref $n_chis angles i
        chi = Dict{Peptides.Dihedral.DihedralType, Tuple{Union{Nothing, Float64}, Float64}}()
        for n in 1:length(a)
            chi[_chis[n]] = (a[n], 0.0)
        end

        function f(pos...)
            rot = Peptides.Rotamer(pose.graph[1, 1].name.content, chi)
            Peptides.apply!(pose.state, rot, pose.graph[1, 1])
            energies[pos...] = energy_function(pose)
            rotamers[pos...] = rot
        end

        @ncall $n_chis f i
    end
end
eval(l)

local_minima         = Images.findlocalminima(energies)
selected_rotamers    = rotamers[local_minima]
selected_energies    = energies[local_minima]
_min                 = minimum(selected_energies)
_max                 = maximum(selected_energies)
minmax_energies      = map((x) -> (x - _min)/(_max - _min), selected_energies)
probabilities        = minmax_energies ./ sum(minmax_energies)
probabilities_sorted = sort(probabilities, rev = true)
perm                 = sortperm(probabilities, rev = true);
rotamer_sorted       = selected_rotamers[perm];

function get_chi(rot::Peptides.Rotamer, chi::Peptides.Dihedral.DihedralType)
    if chi in keys(rot.chis)
        return map(rad2deg, rot.chis[chi])
    else
        return (0.0, 0.0)
    end
end

open(filename, "w") do io
    write(io, @sprintf("%3s %8s %7s %7s %7s %7s %7s %7s %7s %7s\n", "# T",
        "Probabil", "chi1Val", "chi2Val", "chi3Val", "chi4Val", "chi1Sig",
        "chi2Sig", "chi3Sig", "chi4Sig"))

    for (probability, rotamer) in zip(probabilities_sorted, rotamer_sorted)
        probability === 0.0 && continue
        chi1 = get_chi(rotamer, Peptides.Dihedral.chi1)
        chi2 = get_chi(rotamer, Peptides.Dihedral.chi2)
        chi3 = get_chi(rotamer, Peptides.Dihedral.chi3)
        chi4 = get_chi(rotamer, Peptides.Dihedral.chi4)
        write(io, @sprintf("%3s %8.6f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n",
            rotamer.name, probability, chi1[1], chi2[1], chi3[1], chi4[1],
            chi1[2], chi2[2], chi3[2], chi4[2]))
    end
end