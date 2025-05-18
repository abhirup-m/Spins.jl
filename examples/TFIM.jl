include("../src/Dmrg.jl")
using Plots

####  H = -\sum_i S_i^z S^z_{i+1} - g\sum_i S_i^z  ####

maxSize = 25
steps = 20
yvals = []
initSites = 2
for g in 0.5:0.2:1.5
    initHam = [("zz", [i, i+1], -1.) for i in 1:(initSites-1)]
    append!(initHam, [("+", [i], -g) for i in 1:initSites])
    append!(initHam, [("-", [i], -g) for i in 1:initSites])
    bond = i -> [("zz", [i, i+1], -1.)]
    loc = i -> [("+", [i+1], -g/2), ("-", [i+1], -g/2), ("z", [i+1], 0.0),]
    results = InfiniteDMRG(initHam, bond, loc, steps, maxSize)
    push!(yvals, results["energy"][end])
end
println(yvals)
#=p = scatter(yvals)=#
#=display(p)=#
