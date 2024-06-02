using Spins, Plots, Measures

probe = Dict(("z", [1])=>1.0); 
plots = []
for numSites in 2:2:6 
    basisStates = getBasisStates(numSites)
    leftDownState = [fill(0, Int(numSites/2)); fill(1, Int(numSites/2))]
    rightDownState = reverse(leftDownState)
    ham = Dict{Tuple{String, Vector{Int64}}, Float64}()
    for i in 1:numSites
        ham[("+-", [i, i==numSites ? 1 : i+1])] = 0.5
        ham[("+-", [i==numSites ? 1 : i+1, i])] = 0.5
        ham[("zz", [i, i==numSites ? 1 : i+1])] = 0.25
    end
    p = plot(thickness_scaling=1.1, linewidth=4, label="N=$numSites", size=(300, 150), xlabel="time", ylabel="local \$S_z\$", left_margin=-1mm, top_margin=-1mm, bottom_margin=0mm) 
    for state in [initState, reverse(initState)]
        mele = timeEvolvedMatrixElement(basisStates, Dict(BitVector(state) => 1.0), probe, ham, 30.0, 0.1)
        plot!(p, mele)
    end
    push!(plots, p)
    plot(plots...)
end
savefig("spinThermal.pdf")
