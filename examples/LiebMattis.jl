using Plots
include("../src/Base.jl")
include("../src/Models.jl")
include("../src/Eigen.jl")
include("../src/Correlations.jl")

Es = []
vals = [0.5]
N = 3
basis = BasisStates(2 * N)
S = []
for intraCoupling in vals
   H = LiebMattis(N, intraCoupling; globalField=0.)
   E, X = Spectrum(H, basis, symmetries="")
   #=E = E[0]=#
   println(E[1])
   #=X = X[0]=#
   #=println(X[3])=#
   degen = count(e -> e ≤ minimum(E) + 1e-8, E)
   println(degen)
   for Xi in X[E .≤ minimum(E) + 1e-8]
       rhoA = ReducedDensityMatrix(Xi, collect(N+1:2*N))
       push!(S, EntEntropy(rhoA))
   end
   println(sort(S))
   push!(Es, E[1])
end
scatter(vals, Es)
