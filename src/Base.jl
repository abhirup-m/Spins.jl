function getBasisStates(numSpins::Int64)
    basisStates = Dict{Float64, Vector{BitArray}}()
    allSequences = digits.(0:2^numSpins-1, base=2, pad=numSpins) |> reverse
    allSequencesTotSigmaz = sum.([sequence .- 0.5 for sequence in allSequences]) .* 2
    for sigmaz in -numSpins:2:numSpins
        basisStates[sigmaz] = allSequences[allSequencesTotSigmaz .== sigmaz]
    end
    return basisStates
end


function TransformBit(qubit::Bool, operator::Char)
    @assert operator in ('z', '+', '-')
    if operator == 'z'
        return qubit, 2 * (qubit - 0.5)
    elseif (operator == '+' && qubit == 0) || (operator == '-' && qubit == 1)
        return 1 - qubit, 1
    else
        return qubit, 0
    end
end


function applyOperatorOnState(stateDict::Dict{BitVector,Float64}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    @assert maximum([maximum(opMembers) for (_, opMembers) in keys(operatorList)]) ≤ length(collect(keys(stateDict))[1])

    # define a dictionary for the final state obtained after applying 
    # the operators on all basis states
    completeOutputState = Dict{BitVector,Float64}()

    # loop over all operator tuples within operatorList
    for ((opType, opMembers), opStrength) in pairs(operatorList)
        @assert length(opMembers) == length(opType) == length(unique(opMembers))

        for (state, coefficient) in stateDict
            # println(state[opMembers], state[opMembers] .== 0, opType)
            # if '-' ∈ opType[state[opMembers] .== 0] || '+' ∈ opType[state[opMembers] .== 1]
            #     continue
            # end

            newState = copy(state)
            newCoefficient = coefficient
            for (operator, spinIndex) in zip(opType, opMembers)
                newQubit, factor = TransformBit(state[spinIndex], operator)
                newState[spinIndex] = newQubit
                newCoefficient *= factor
                if newCoefficient == 0
                    break
                end
            end

            if newCoefficient != 0
                if newState in keys(completeOutputState)
                    completeOutputState[newState] += opStrength * newCoefficient
                else
                    completeOutputState[newState] = opStrength * newCoefficient
                end
            end
        end
    end
    return completeOutputState
end


function generalOperatorMatrix(basisStates::Dict{Float64,Vector{BitArray}}, operatorList::Dict{Tuple{String,Vector{Int64}},Float64})
    operatorFullMatrix = Dict(key => zeros(length(value), length(value)) for (key, value) in basisStates)

    # loop over each symmetry sector
    for (key, bstates) in collect(basisStates)

        # loop over the basis states in the current symmetry sector
        for (index, state) in collect(enumerate(bstates))

            # loop over the new states generated upon applying
            # the operator to the current basis state
            for (nState, coeff) in applyOperatorOnState(Dict(state => 1.0), operatorList)
                operatorFullMatrix[key][index, bstates.==[nState]] .= coeff
            end
        end
    end
    return operatorFullMatrix
end
