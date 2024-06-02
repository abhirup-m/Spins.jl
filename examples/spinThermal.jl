using Spins, Plots, Measures, ProgressMeter

timeMax = 20.0
deltaTime = 0.1
probe = Dict(("z", [1])=>1.0); 
plots = []
numSitesSet = 2:2:12
@showprogress for numSites in numSitesSet
    basisStates = Spins.getBasisStates(numSites)
    leftDownState = [fill(0, trunc(Int, numSites/2)); fill(1, numSites - trunc(Int, numSites/2))]
    rightDownState = reverse(leftDownState)
    hamiltonian = Dict{Tuple{String, Vector{Int64}}, Float64}()
    for i in 1:numSites
        hamiltonian[("+-", [i, i==numSites ? 1 : i+1])] = 0.5
        hamiltonian[("+-", [i==numSites ? 1 : i+1, i])] = 0.5
        hamiltonian[("zz", [i, i==numSites ? 1 : i+1])] = 0.25
    end
    p = plot(annotationfontsize=11, annotate=(timeMax * 0.9, 0.2, "N=$numSites")) 
    for initState in [leftDownState, rightDownState]
        mele = Spins.timeEvolvedMatrixElement(basisStates, Dict(BitVector(initState) => 1.0), probe, hamiltonian, timeMax, deltaTime)
        plot!(p, range(0, stop=timeMax, step=deltaTime), mele, linewidth=4, label="", thickness_scaling=1.3, xlabel="time \$\\rightarrow\$", ylabel="local \$S_z\$")
    end
    push!(plots, p)
end
p = plot(plots..., size=(1100, 1500), layout=(length(numSitesSet), 1))
savefig(p, "spinThermal.pdf")
