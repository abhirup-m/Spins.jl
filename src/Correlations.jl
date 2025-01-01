using LinearAlgebra, ProgressMeter

function ExpecValue(
        state::Dict{BitVector,Float64}, 
        operator::Vector{Tuple{String, Vector{Int64}, Float64}},
    )
    modifiedState = ApplyOperator(operator, state)
    return StateOverlap(state, modifiedState)
end


function timeEvolvedMatrixElement(basisStates::Dict{Float64,Vector{BitArray}}, observeState::Dict{BitVector, Float64}, operator::Dict{Tuple{String,Vector{Int64}},Float64}, hamiltonian::Dict{Tuple{String,Vector{Int64}},Float64}, timeMax::Float64, deltaTime::Float64)
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


function ReducedDensityMatrix(
        state::Dict{BitVector,Float64}, 
        nonTracedIndices::Vector{Int64},
    )
    stateNorm = sum(values(state) .^ 2 )^0.5
    normalisedState = deepcopy(state)
    map!(v -> v / stateNorm, values(normalisedState))
    nonTracedBasis = BasisStates(length(nonTracedIndices))
    reducedDensityMatrix = zeros(length(nonTracedBasis), length(nonTracedBasis))
    tracedIndices = [i for i in normalisedState |> keys |> first |> length if i ∉ nonTracedIndices]
    for (bstate1, val1) in normalisedState
        for (bstate2, val2) in normalisedState
            if bstate1[tracedIndices] == bstate2[tracedIndices]
                nonTracedState1 = Dict(bstate1[nonTracedIndices] => 1.)
                nonTracedState2 = Dict(bstate2[nonTracedIndices] => 1.)
                reducedDensityMatrix[findfirst(==(nonTracedState1), nonTracedBasis), 
                                     findfirst(==(nonTracedState2), nonTracedBasis),
                                    ] += val1' * val2
            end
        end
    end
    return reducedDensityMatrix
end
