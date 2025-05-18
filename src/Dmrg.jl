using LinearAlgebra, ProgressMeter

"""
Slice a Hamiltonian into multiple parts in order to
constitute a Hamiltonian flow that can be studied
using iterative diagonalisation. The slicing happens
at the indices provided in indexParitions; for each
index i, we retain all terms in the Hamiltonian that
act on i. If this is the first index, we retain all
terms that act only on indices upto i. For example, given
a Hamiltonian H = c^†_1 c_2 + c^†_2 c_3 + c^†_3 c_4 + h.c.
and indexPartitions = [2, 4], this function returns two
Hamiltonians, one for the first step H1 = c^†_1 c_2 and
one for the second step H2 = c^†_2 c_3 + c^†_3 c_4.
"""
function MinceHamiltonian(
        hamiltonian::Vector{Tuple{String, Vector{Int64}, Float64}},
        indexPartitions::Vector{Int64},
    )
    
    # first create the subspaces by defining start and stop
    # indices for each subspace.
    subspaces = NTuple{2, Int64}[]
    leadingIndex = 0
    for index in indexPartitions
        push!(subspaces, (leadingIndex+1, index))
        leadingIndex = index
    end

    # check that the subspaces are mutually exclusive
    @assert issorted(vcat(subspaces...))

    # allocate the Hamiltonian for each step of partitioning
    hamFlow = Vector{Tuple{String, Vector{Int64}, Float64}}[[] for _ in indexPartitions]
    for (opType, members, coupling) in hamiltonian

        # for each term in the Hamiltonian, check in which subspace
        # its farthest index lies.
        subspaceIndex = findfirst(s -> s[1] ≤ maximum(members) ≤ s[2], subspaces)
        push!(hamFlow[subspaceIndex], (opType, members, coupling))
    end
    return hamFlow
end
export MinceHamiltonian


function Operators(
        numSites::Int64;
        rev::Bool=false,
        startFrom::Int64=1,
        padBy::Int64=1,
    )
    operators = Dict{Tuple{Int64, Char}, Matrix{Float64}}()
    for site in 1:numSites
        sigmaPlus = [[0, 1] [0, 0]]
        leftPadding = I(2^(site-1))
        rightPadding = I(2^(numSites - site))
        if rev
            operators[(site + startFrom - 1, '+')] = kron(kron(I(padBy), rightPadding), sigmaPlus, leftPadding)
        else
            operators[(site + startFrom - 1, '+')] = kron(leftPadding, sigmaPlus, kron(rightPadding, I(padBy)))
        end
        operators[(site + startFrom - 1, '-')] = operators[(site + startFrom - 1, '+')]'
        operators[(site + startFrom - 1, 'z')] = 2 .* (operators[(site + startFrom - 1, '+')] * operators[(site + startFrom - 1, '-')] .- 0.5 .* I(operators[(site + startFrom - 1, '-')] |> size |> first))
    end
    return operators
end


function GetMatrix(
        hamiltonian::Vector{Tuple{String, Vector{Int64}, Float64}},
        operators::Dict{Tuple{Int64, Char}, Matrix{Float64}},
    )
    matrix = zeros(size(operators |> values |> collect |> first)...)
    for (action, members, coupling) in hamiltonian
        matrix += coupling * prod([operators[chunk] for chunk in zip(members, action)])
    end
    return matrix
end


function GetMatrix(
        hamiltonian::Vector{Tuple{String, Vector{Int64}, Float64}},
        leftOperators::Dict{Tuple{Int64, Char}, Matrix{Float64}},
        rightOperators::Dict{Tuple{Int64, Char}, Matrix{Float64}},
        partition::Int64,
    )
    matrix = zeros((size(leftOperators |> values |> collect |> first) .* size(rightOperators |> values |> collect |> first))...)
    for (action, members, coupling) in hamiltonian
        term = I(matrix |> size |> first)
        for chunk in zip(members, action)
            if chunk[1] ≤ partition
                term *= kron(leftOperators[chunk], I(leftOperators[chunk] |> size |> first)) 
            else
                location = 2 * partition - chunk[1] + 1
                term *= kron(I(rightOperators[(location, chunk[2])] |> size |> first), rightOperators[(location, chunk[2])])
            end
        end
        matrix += coupling * term
    end
    return matrix
end


function ReflectHamiltonian(
        hamiltonian::Vector{Tuple{String, Vector{Int64}, Float64}},
        lastSite::Int64,
    )
    reflectedHamiltonian = Tuple{String, Vector{Int64}, Float64}[]
    for (action, members, coupling) in hamiltonian
        reflectMembers = lastSite .- members 
        push!(reflectedHamiltonian, (action, reflectMembers, coupling))
    end
    return reflectedHamiltonian
end


function RotationMatrices(
        evecs::Matrix{Float64},
        maxSize::Int64,
    )
    totalDimension = trunc.(Int, size(evecs).^0.5)
    gstateTensor = reshape(evecs[:, 1], totalDimension)
    rhoSystem = gstateTensor * gstateTensor'
    evals, basis = eigen(Hermitian(rhoSystem)) 
    if length(evals) < maxSize
        maxSize = length(evals)
    end
    rotationSystem = basis[:, sortperm(evals, rev=true)[1:maxSize]]

    rhoEnv = gstateTensor' * gstateTensor
    evals, basis = eigen(Hermitian(rhoEnv)) 
    rotationEnv = basis[:, sortperm(evals, rev=true)[1:maxSize]]
    return rotationSystem, rotationEnv
end


function InfiniteDMRG(
        initHamiltonian::Vector{Tuple{String, Vector{Int64}, Float64}},
        bondHamiltonian::Function,
        locHamiltonian::Function,
        numSteps::Int64,
        maxSize::Int64,
    )
    currentSites = unique(vcat([members for (_, members, _) in initHamiltonian]...))
    leader = maximum(currentSites)
    addedHamiltonian = vcat(locHamiltonian(leader), bondHamiltonian(leader))
    bondSites = unique(vcat([members for (_, members, _) in addedHamiltonian]...))
    bondEnd = maximum(bondSites)
    sysOperators = Operators(bondEnd)
    envOperators = Operators(bondEnd; rev=true)

    step = 1
    system = copy(initHamiltonian)
    display(GetMatrix(system, Operators(8)))
    E, _ = eigen(Hermitian(GetMatrix(system, Operators(8))))
    println(E[1])
    sysExpand = vcat(system, addedHamiltonian)
    system = GetMatrix(sysExpand, sysOperators)
    environment = GetMatrix(sysExpand, envOperators)
    results = Dict("energy" => Float64[])
    @showprogress for step in 1:numSteps
        sysEnvCoupling = GetMatrix(bondHamiltonian(bondEnd), sysOperators, envOperators, bondEnd)
        superHamiltonian = kron(system, I(environment |> size |> first)) + kron(I(system |> size |> first), environment) + sysEnvCoupling
        evals, evecs = eigen(Hermitian(superHamiltonian))
        push!(results["energy"], evals[1] / (2 * bondEnd))

        rotationSys, rotationEnv = RotationMatrices(evecs, maxSize)

        leader = bondEnd

        oldSize = size(rotationSys' * system * rotationSys)

        addedHamiltonian = vcat(locHamiltonian(leader), bondHamiltonian(leader))
        bondSites = filter(>(leader), vcat([members for (_, members, _) in bondHamiltonian(leader)]...)) |> unique
        bondEnd = maximum(bondSites)
        for ((k1, v1), (k2, v2)) in zip(sysOperators, envOperators)
            sysOperators[k1] = kron(rotationSys' * v1 * rotationSys, I(2^length(bondSites)))
            envOperators[k2] = kron(I(2^length(bondSites)), rotationEnv' * v2 * rotationEnv)
        end

        merge!(sysOperators, Operators(length(bondSites); startFrom=leader+1, padBy=oldSize[1]))
        merge!(envOperators, Operators(length(bondSites); startFrom=leader+1, rev=true, padBy=oldSize[1]))

        addedMatrixSys = GetMatrix(addedHamiltonian, sysOperators)
        addedMatrixEnv = GetMatrix(addedHamiltonian, envOperators)

        system = kron(rotationSys' * system * rotationSys, I(2^length(bondSites)))
        environment = kron(I(2^length(bondSites)), rotationEnv' * environment * rotationEnv)
        system .+= addedMatrixSys
        environment .+= addedMatrixEnv
    end
    return results
end
export InfiniteDMRG
