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
numSpins = [3]
N = 3
for offCritCoupling in [1.3,]
    basis = BasisStates(2 * N)
    H = LiebMattis(N; totSpinSqCoupling=offCritCoupling)
    E, X = Spectrum(H, basis)
    #=@time E, X = Spectrum(H, basis; symmetries="Z", symmetrySector=[0])=#
    #=minSector = E |> keys |> collect .|> abs |> minimum=#
    #=E = E[minSector]=#
    #=X = X[minSector]=#
    println(E[1])
    display(sort(values(X[1]) |> collect))
    degen = count(e -> e ≤ minimum(E) + 1e-8, E)
    println(degen)
    vne = []
    for Xi in X[E .≤ minimum(E) + 1e-8]
        rhoA = ReducedDensityMatrixFast(Xi, collect(1:N))
        push!(vne, EntEntropy(rhoA))
    end
    println(vne)
    push!(vneMax, maximum(vne))
    push!(gsEnergy, E[1])
end
#=println(vneMax)=#
scatter(numSpins, vneMax)
#=scatter(totSpinSqCouplingVals, vneMax)=#
