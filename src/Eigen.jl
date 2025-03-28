using LinearAlgebra

function Spectrum(
        hamiltonian::Vector{Tuple{String,Vector{Int64},Float64}},
        basisStates::Vector{Dict{BitVector,Float64}};
        symmetries::String="",
        tolerance::Float64=1e-10,
    )
    if isempty(symmetries)
        hamiltonianMatrix = OperatorMatrix(basisStates, hamiltonian)
        eigenValues, X = eigen(Hermitian(hamiltonianMatrix))
        eigenStates = eltype(basisStates)[] 
        for eigenState in eachcol(X)
            eigenState = mergewith(+, [Dict(k => v * cj for (k,v) in basisStates[j]) for (j, cj) in enumerate(eigenState) if abs(cj) > tolerance]...)
            push!(eigenStates, eigenState)
        end
    elseif symmetries == "Z"
        numSpins = basisStates[1] |> keys |> collect |> first |> length
        totSzValues = range(-numSpins, numSpins, step=2) |> collect
        classifiedBasis = Dict{Int64, typeof(basisStates)}(m => eltype(basisStates)[] for m in totSzValues)
        for state in basisStates
            totSz = 2 * count(==(1), state |> keys |> collect |> first) - numSpins
            push!(classifiedBasis[totSz], state)
        end
        eigenStates = Dict{Int64, Vector{Vector{Float64}}}(m => Vector{Float64}[] for m in totSzValues)
        eigenValues = Dict{Int64, Vector{Float64}}(m => Float64[] for m in totSzValues)
        for totSz in totSzValues
            hamiltonianMatrix = OperatorMatrix(classifiedBasis[totSz], hamiltonian)
            eigenValues[totSz], X = eigen(Hermitian(hamiltonianMatrix))
            eigenStates[totSz] = [Xi for Xi in eachcol(X)]
        end
    end
    return eigenValues, eigenStates
end
export Spectrum
