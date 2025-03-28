using Plots
include("storage/programmingProjects/Spins/src/Base.jl")
include("storage/programmingProjects/Spins/src/Models.jl")
include("storage/programmingProjects/Spins/src/Eigen.jl")
include("storage/programmingProjects/Spins/src/Correlations.jl")

Es = []
vals = 0:10:100
N = 3
basis = BasisStates(2 * N)
for intraCoupling in vals
   H = LiebMattis(N, 0.5 * intraCoupling; globalField=1e-4)
   E, X = Spectrum(H, basis, symmetries="")
   S = 0
   degen = count(e -> e ≤ minimum(E) + 1e-8, E)
   println(degen)
   for Xi in X[E .≤ minimum(E) + 1e-8]
       rhoA = ReducedDensityMatrixFast(Xi, collect(N+1:2*N))
       #=eigVals, _ = eigen(rhoA)=#
       #=println(histogram(eigVals))=#
       S = ifelse(S > EntEntropy(rhoA), S, EntEntropy(rhoA))
       #=S = ifelse(S < EntEntropy(rhoA), S, EntEntropy(rhoA))=#
   end
   println(S)
   push!(Es, S)
end
scatter(vals, Es)
