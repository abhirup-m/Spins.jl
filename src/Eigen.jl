using LinearAlgebra

function Spectrum(
        hamiltonian::Vector{Tuple{String,Vector{Int64},Float64}},
        basisStates::Vector{Dict{BitVector,Float64}};
        symmetries::String="",
        symmetrySector::Vector{Int64}=Int64[],
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
        if !isempty(symmetrySector)
            totSzValues = intersect(totSzValues, symmetrySector)
        end
        classifiedBasis = Dict{Int64, typeof(basisStates)}(m => eltype(basisStates)[] for m in totSzValues)
        for state in basisStates
            totSz = 2 * count(==(1), state |> keys |> collect |> first) - numSpins
            if totSz in totSzValues
                push!(classifiedBasis[totSz], state)
            end
        end
        eigenValues = Dict{Int64, Vector{Float64}}(m => Float64[] for m in totSzValues)
        eigenStates = Dict{Int64, typeof(basisStates)}(m => eltype(basisStates)[] for m in totSzValues)
        for totSz in totSzValues
            hamiltonianMatrix = OperatorMatrix(classifiedBasis[totSz], hamiltonian)
            eigenValues[totSz], X = eigen(Hermitian(hamiltonianMatrix))
            for eigenState in eachcol(X)
                eigenState = mergewith(+, [Dict(k => v * cj for (k,v) in classifiedBasis[totSz][j]) for (j, cj) in enumerate(eigenState) if abs(cj) > tolerance]...)
                push!(eigenStates[totSz], eigenState)
            end
        end
    end
    return eigenValues, eigenStates
end
export Spectrum
