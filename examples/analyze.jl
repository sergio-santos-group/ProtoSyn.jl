using ProtoSyn
using ProtoSyn.Peptides
using Clustering
ENV["GKSwstype"]="100"
using Plots
using Plots.PlotMeasures
using PlotThemes
theme(:lime)

poses       = Vector{Pose}()
energies    = []
N           = 0

# Load poses based on table.csv
open("table.csv", "r") do file_in
    global N
    for line in eachline(file_in)
        push!(energies, parse(Float64, split(line)[6]))
        N += 1 
        push!(poses, ProtoSyn.Peptides.load("monte_carlo_$N.pdb"))
    end
end

# Calculate the RMSD matrix
rmsd_matrix = ProtoSyn.Clustering.rmsd_matrix(poses, an"CA")

function plot_results(assignment_function::Function, title::String)

    cdd_dbi, add_dbi = Vector{Float64}(), Vector{Float64}()
    cdd_di, add_di   = Vector{Float64}(), Vector{Float64}()
    lowest_count, highest_count = Vector{Int64}(), Vector{Int64}()

    for K in 2:N
        a = assignment_function(K)
        push!(lowest_count, minimum(counts(a))) 
        push!(highest_count, maximum(counts(a))) 

        # Get the complete- and average- diameter distances
        intra_cdd, inter_cdd = ProtoSyn.Clustering.complete_diameter_distances(a, rmsd_matrix)
        intra_add, inter_add = ProtoSyn.Clustering.average_diameter_distances(a, rmsd_matrix)

        # Calculate the complete- and average- davies-bouldin and dunn indexes
        push!(cdd_dbi, ProtoSyn.Clustering.davies_bouldin_index(intra_cdd, inter_cdd))
        push!(add_dbi, ProtoSyn.Clustering.davies_bouldin_index(intra_add, inter_add))

        push!(cdd_di, ProtoSyn.Clustering.dunn_index(intra_cdd, inter_cdd))
        push!(add_di, ProtoSyn.Clustering.dunn_index(intra_add, inter_add))
    end

    l = @layout [b;d]

    b = plot(2:N, cdd_dbi, label = "Complete Davies-Bouldin", xticks = collect(2:2:N), legend = :outertopright, fg_legend = :transparent, title = title, titlefontsize = 10)
    plot!(2:N, add_dbi, label = "Average Davies-Bouldin")
    plot!(2:N, cdd_di, label = "Complete Dunn")
    plot!(2:N, add_di, label = "Average Dunn")
    xlabel!("Nº clusters", xguidefontsize = 8)
    ylabel!("Index", yguidefontsize = 8)

    d = bar(2:N, lowest_count, widen = false, xticks = collect(2:2:N), linealpha = 0.0, c = :hotpink, bar_width = 0.2, label = "Lowest count cluster", bg_legend = :transparent, fg_legend = :transparent)
    bar!(1.8:1:N, highest_count, widen = false, linealpha = 0.0, bar_width = 0.2, label = "Highest count cluster")
    xlabel!("Nº clusters", xguidefontsize = 8)
    ylabel!("Nº of elements", yguidefontsize = 8)

    return plot(b, d, layout = l, link = :x)
end

function slh(K::Int)
    c = hclust(rmsd_matrix, linkage = :single, branchorder = :barjoseph)
    return cutree(c; k = K)
end

function alh(K::Int)
    c = hclust(rmsd_matrix, linkage = :average, branchorder = :barjoseph)
    return cutree(c; k = K)
end

function clh(K::Int)
    c = hclust(rmsd_matrix, linkage = :complete, branchorder = :barjoseph)
    return cutree(c; k = K)
end

function wlh(K::Int)
    c = hclust(rmsd_matrix, linkage = :ward, branchorder = :barjoseph)
    return cutree(c; k = K)
end

function km(K::Int)
    return kmedoids(rmsd_matrix, K, init = :kmcen).assignments
end

p1 = plot_results(slh, "Single linkage hierarchical")
p2 = plot_results(alh, "Average linkage hierarchical")
p3 = plot_results(clh, "Complete linkage hierarchical")
p4 = plot_results(wlh, "Ward linkage hierarchical")
p5 = plot_results(km, "K-medoids")
l = @layout [p1 p2 p3 p4 p5]
plot(p1, p2, p3, p4, p5, size = (3600, 720), margins = 10mm, layout = l)
savefig("fig.png")

plot(p1, size = (1080, 720), margins = 10mm)
savefig("p1.png")
plot(p2, size = (1080, 720), margins = 10mm)
savefig("p2.png")
plot(p3, size = (1080, 720), margins = 10mm)
savefig("p3.png")
plot(p4, size = (1080, 720), margins = 10mm)
savefig("p4.png")
plot(p5, size = (1080, 720), margins = 10mm)
savefig("p5.png")

clusters = kmedoids(rmsd_matrix, 20, init = :kmcen)

rg  = findmax(clusters.counts)[2]
hits = clusters.counts.==rg
gmin_e = Inf

if count(hits) != 1
    for (gi, hit) in enumerate(hits)
        if hit
            min_e = minimum(energies[clusters.assignments.==gi])
            if min_e < gmin_e
                rg = gi
                gmin_e = min_e
            end
        end
    end
end

gy = zeros(N)
ngy = zeros(N)
for (i, v) in enumerate(clusters.assignments)
    if v == rg
        gy[i] = energies[i]
    else
        ngy[i] = energies[i]
    end
end

begin
    l = @layout [a{0.1w} b{0.2h};c{0.1w} d]
    _hbar = bar(1:N, gy, legend = :none, widen = false, xaxis = nothing, linealpha = 0.0, c = :hotpink)
    bar!(1:N, ngy, legend = :none, widen = false, xaxis = nothing, linealpha = 0.0)

    _vbar = bar(1:N, gy, mirror = true, legend = :none, widen = false, yaxis = nothing, orientation = :horizontal, linealpha = 0.0, c = :hotpink)
    bar!(1:N, ngy, mirror = true, legend = :none, widen = false, yaxis = nothing, orientation = :horizontal, linealpha = 0.0)

    plot(
        plot(legend=false, grid=false, framestyle = :none), # Empty space
        _hbar,
        _vbar,
        heatmap(1:1:N, 1:1:N, rmsd_matrix, colorbar = :none, widen = false),
        layout = l, link = :x)
    plot!(size=(1080, 1080))
    savefig("fig2.png")
end


rg  = clusters.assignments[findmin(energies)[2]]

gy = zeros(N)
ngy = zeros(N)
for (i, v) in enumerate(clusters.assignments)
    if v == rg
        gy[i] = energies[i]
    else
        ngy[i] = energies[i]
    end
end

begin
    l = @layout [a{0.1w} b{0.2h};c{0.1w} d]
    _hbar = bar(1:N, gy, legend = :none, widen = false, xaxis = nothing, linealpha = 0.0, c = :hotpink)
    bar!(1:N, ngy, legend = :none, widen = false, xaxis = nothing, linealpha = 0.0)

    _vbar = bar(1:N, gy, mirror = true, legend = :none, widen = false, yaxis = nothing, orientation = :horizontal, linealpha = 0.0, c = :hotpink)
    bar!(1:N, ngy, mirror = true, legend = :none, widen = false, yaxis = nothing, orientation = :horizontal, linealpha = 0.0)

    plot(
        plot(legend=false, grid=false, framestyle = :none), # Empty space
        _hbar,
        _vbar,
        heatmap(1:1:N, 1:1:N, rmsd_matrix, colorbar = :none, widen = false),
        layout = l, link = :x)
    plot!(size=(1080, 1080))
    savefig("fig3.png")
end