using Plots
include("../src/Base.jl")
include("../src/Models.jl")
include("../src/Eigen.jl")
include("../src/Correlations.jl")

gsEnergy = []
totSpinSqCouplingVals = 0.6:0.1:10
critCoupling = 1
offCritCoupling = 1.5
vneMax = []
numSpins = 2:5
for N in numSpins
    basis = BasisStates(2 * N)
    H = LiebMattis(N; totSpinSqCoupling=critCoupling)
    #=E, X = Spectrum(H, basis)=#
    @time E, X = Spectrum(H, basis; symmetries="Z", symmetrySector=[0])
    minSector = E |> keys |> collect .|> abs |> minimum
    E = E[minSector]
    X = X[minSector]
    degen = count(e -> e ≤ minimum(E) + 1e-8, E)
    println(degen)
    vne = []
    for Xi in X[E .≤ minimum(E) + 1e-8]
        rhoA = ReducedDensityMatrix(Xi, collect(1:N))
        push!(vne, EntEntropy(rhoA))
    end
    #=println(vne)=#
    push!(vneMax, maximum(vne))
    #=display(X[sortperm(vne)][1])=#
    push!(gsEnergy, E[1])
end
#=println(vneMax)=#
scatter(numSpins, vneMax)
#=scatter(totSpinSqCouplingVals, vneMax)=#
