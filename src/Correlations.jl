using LinearAlgebra

function timeEvolvedMatrixElement(basisStates::Dict{Float64,Vector{BitArray}}, observeState::Dict{BitVector,Float64}, operator::Dict{Tuple{String,Vector{Int64}},Float64}, hamiltonian::Dict{Tuple{String,Vector{Int64}},Float64}, timeMax::Float64, deltaTime::Float64)
    hamiltonianMatrix = generalOperatorMatrix(basisStates, hamiltonian)
    propagatorSectors = Dict(sector => (I + (hamiltonianMatrix[sector] .* 1im .* deltaTime ./ 2)) / (I - (hamiltonianMatrix[sector] .* 1im .* deltaTime ./ 2)) for (sector, matrix) in hamiltonianMatrix)
    matrixElement = Float64[]
    observeStateClassified = Dict()
    for (sector, states) in basisStates
        observeStateClassified[sector] = Dict()
        weights = [state ∈ keys(observeState) ? observeState[state] : 0 for state in states]
        observeStateClassified[sector] = weights
    end
    operatorMatrix = generalOperatorMatrix(basisStates, operator)
    @showprogress for time in range(0, stop=timeMax, step=deltaTime)
        push!(matrixElement, 0)
        for sector in keys(observeStateClassified)
            matrixElement[end] += real(observeStateClassified[sector]' * operatorMatrix[sector] * observeStateClassified[sector])
        end
        for sector in keys(operatorMatrix)
            operatorMatrix[sector] = propagatorSectors[sector] * operatorMatrix[sector] * propagatorSectors[sector]'
        end
    end
    return matrixElement
end

