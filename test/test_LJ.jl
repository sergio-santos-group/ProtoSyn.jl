using ProtoSyn

function test_LJ(;workload::Int64 = 100, init_distance::Float64 = 1.0)

    state = Common.State(2)
    state.xyz = [0.0 0.0 0.0; init_distance 0.0 0.0]
    d = init_distance
    data::Vector{Float64} = []
    distances::Vector{Float64} = []

    #                              SIGMA (σ)            EPSILON (ϵ)
    atoms = [
        Forcefield.Amber.Atom("N", 0.32499999149873800, 0.71127992436207890, 0.0, Int64[], Int64[]),
        Forcefield.Amber.Atom("H", 0.10690800172161044, 0.06568879294383578, 0.0, Int64[], Int64[])]
        
    for i in 1:workload
        println("\nDistance $i: $(ceil(d, digits = 1))")
        Forcefield.Amber.evaluate!(atoms, state, cut_off = 1.2, eCoulomb_λ = 0.0)
        println("Calculated: $(state.energy.comp["eLJ"])")

        ϵ = sqrt(0.71127992436207890 * 0.06568879294383578)
        σ = (0.32499999149873800 + 0.10690800172161044) / 2
        e = 4 * ϵ * ((σ/d)^12 - (σ/d)^6)
        println("    Teoric: $e")
        push!(data, state.energy.comp["eLJ"])
        push!(distances, d)
        d = -(init_distance/workload * i)+init_distance
        state.xyz[2, :] = [d 0.0 0.0]
    end
    # println(distances)
    # println(data)
end

test_LJ()