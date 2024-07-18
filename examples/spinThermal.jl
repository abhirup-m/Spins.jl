using Plots, Measures, ProgressMeter
include("../src/Base.jl")
include("../src/Correlations.jl")
theme(:dark)
default(linewidth=3, thickness_scaling=1.5, leftmargin=-5mm, bottommargin=-4mm)
timeMax = 30.0
deltaTime = 0.2
probe = Dict(("z", [1]) => 0.5);
numSitesSet = 4:6:16
saveNames = []
for numSites in numSitesSet
    basisStates = getBasisStates(numSites; magz=0)
    leftDownState = [fill(0, div(numSites, 2)); fill(1, numSites - div(numSites, 2))]
    rightDownState = reverse(leftDownState)
    hamiltonian = Dict{Tuple{String,Vector{Int64}},Float64}()
    for i in 1:numSites
        hamiltonian[("+-", [i, i == numSites ? 1 : i + 1])] = 0.5
        hamiltonian[("+-", [i == numSites ? 1 : i + 1, i])] = 0.5
        hamiltonian[("zz", [i, i == numSites ? 1 : i + 1])] = 0.25
    end
    p = plot(
    )
    mele = [[], []]
    @time for (i, initState) in enumerate([leftDownState, rightDownState])
        mele[i] = timeEvolvedMatrixElement(basisStates, Dict(BitVector(initState) => 1.0), probe, hamiltonian, timeMax, deltaTime)
    end
    plot!(p, range(0, stop=timeMax, step=deltaTime), mele;
        legend=:topright, title="\$N=$(numSites)\$",
        ylabel="\$\\langle\\mathcal{O}(t)\\rangle\$", xlabel="time \$(t)\$",
        labels=["\$\\uparrow\\uparrow\\downarrow\\downarrow\$" "\$\\downarrow\\downarrow\\uparrow\\uparrow\$"]
    )
    hline!(p, [0]; linestyle=:dot, label="\$\\mathcal{\\bar{O}}(E)\$")
    savefig(p, "spinThermal_$(numSites).pdf")
    push!(saveNames, "spinThermal_$(numSites).pdf")
end
run(`pdfunite $saveNames spinThermal.pdf`)
run(`rm $saveNames`)
