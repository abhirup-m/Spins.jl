include("../src/Base.jl")
include("../src/Models.jl")
include("../src/Eigen.jl")
include("../src/Correlations.jl")

N = 2
basis = BasisStates(N)
model = Heisenberg(N);
E, X = Spectrum(model, basis)
println(E)
