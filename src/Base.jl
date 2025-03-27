function BasisStates(
        numSpins::Int64
    )
    basis = Dict{BitVector,Float64}[]
    for decimalNum in 0:2^numSpins-1
        config = digits(decimalNum, base=2, pad=numSpins) |> reverse
        push!(basis, Dict(BitVector(config) => 1.0))
    end
    return basis
end
export BasisStates


function TransformBit(
        qubit::Bool, 
        operator::Char
    )
    @assert operator in ('z', '+', '-')
    if operator == 'z'
        return qubit, 2 * (qubit - 0.5)
    elseif (operator == '+' && qubit == 0) || (operator == '-' && qubit == 1)
        return 1 - qubit, 1
    else
        return qubit, 0
    end
end
export TransformBit


function ApplyOperatorChunk(
        opType::String,
        opMembers::Vector{Int64},
        opStrength::Float64,
        incomingState::Dict{BitVector,Float64};
        tolerance::Float64=1e-16
    )
    outgoingState = Dict{BitVector,Float64}()
    for i in eachindex(opMembers)[end:-1:2]
        if opMembers[i] ∈ opMembers[i+1:end]
            continue
        end
        if opType[i] == '+'
            filter!(p -> p[1][opMembers[i]] == 0, incomingState)
        elseif opType[i] == '-'
            filter!(p -> p[1][opMembers[i]] == 1, incomingState)
        end
    end
    for (incomingBasisState, coefficient) in incomingState

        newCoefficient = coefficient
        outgoingBasisState = copy(incomingBasisState)

        # for each basis state, obtain a modified state after applying the operator tuple
        for (siteIndex, operator) in zip(reverse(opMembers), reverse(opType))
            newQubit, factor = TransformBit(outgoingBasisState[siteIndex], operator)
            if factor == 0
                newCoefficient = 0
                break
            end
            outgoingBasisState[siteIndex] = newQubit
            newCoefficient *= factor
        end

        if abs(newCoefficient) > tolerance
            if haskey(outgoingState, outgoingBasisState)
                outgoingState[outgoingBasisState] += opStrength * newCoefficient
            else
                outgoingState[outgoingBasisState] = opStrength * newCoefficient
            end
        end
    end
    return outgoingState
end
export ApplyOperatorChunk


function ApplyOperator(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        incomingState::Dict{BitVector,Float64};
        tolerance::Float64=1e-16
    )
    @assert !isempty(operator)
    @assert maximum([maximum(positions) for (_, positions, _) in operator]) ≤ length.(keys(incomingState))[1]

    return mergewith(+, fetch.([Threads.@spawn ApplyOperatorChunk(opType, opMembers, opStrength, copy(incomingState); tolerance=tolerance) 
                                for (opType, opMembers, opStrength) in operator])...)

    return outgoingState
end
export ApplyOperator


function OperatorMatrix(
        basisStates::Vector{Dict{BitVector,Float64}},
        operator::Vector{Tuple{String,Vector{Int64},Float64}};
        tolerance::Float64=1e-16,
    )
    operatorMatrix = zeros(length(basisStates), length(basisStates))
    newStates = fetch.([Threads.@spawn ApplyOperator(operator, incomingState; tolerance=tolerance)
                        for incomingState in basisStates])
    Threads.@threads for incomingIndex in findall(!isempty, newStates)
        Threads.@threads for outgoingIndex in eachindex(basisStates)
            operatorMatrix[outgoingIndex, incomingIndex] = StateOverlap(basisStates[outgoingIndex], newStates[incomingIndex])
        end
    end
    return operatorMatrix
end
export OperatorMatrix


function Commutator(
        basisStates::Vector{Dict{BitVector,Float64}},
        operatorLeft::Vector{Tuple{String,Vector{Int64},Float64}},
        operatorRight::Vector{Tuple{String,Vector{Int64},Float64}},
    )
    matrixLeft = OperatorMatrix(basisStates, operatorLeft)
    matrixRight = OperatorMatrix(basisStates, operatorRight)
    return matrixLeft * matrixRight - matrixRight * matrixLeft
end
export Commutator


function StateOverlap(
        stateBra::Dict{BitVector,Float64}, 
        stateKet::Dict{BitVector,Float64}
    )
    overlap = 0.
    keys2 = keys(stateBra)
    for (key, val) in stateKet
        if key ∈ keys2
            overlap += val * stateBra[key]'
        end
    end
    return overlap
end
export StateOverlap
