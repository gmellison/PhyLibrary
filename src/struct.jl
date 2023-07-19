##########################################################
# struct
############################################################
struct SEQdata
    taxaname::Vector{String}
    sequence::Vector{String}
    gene::Matrix{String}

    function SEQdata(taxaname::Vector{String}, sequence::Vector{String}, gene::Matrix{String})
        new(taxaname, sequence, gene)
    end

    function SEQdata(taxaname::Vector{String}, sequence::Vector{String})
        new(taxaname, sequence)
    end
end

mutable struct Node 
    father::Int64
    nson::Int64
    sons::Vector{Int64}
    taxaname::AbstractString
    brlens::Float64
    theta::Float64
    label::String

    function Node()
        new(-9, 0, [-9,-9], "", -9.0, -9.0,"")
    end

    function Node(rootnode::Bool)
        if rootnode
            new(-9, 0, [-9,-9,-9], "", -9.0, -9.0,"")
        else
            new(-9, 0, [-9,-9], "", -9.0, -9.0,"")
        end
    end

    function Node(father::Int64, nson::Int64, sons::Vector{Int64}, taxaname::AbstractString, brlens::Float64, theta::Float64, label::String)
        new(father, nson, sons, taxaname, brlens, theta, label)
    end
end

mutable struct Tree
    rootnode::Int64
    ntaxa::Int64 
    nodes::Array{Node,1}
    rooted::Bool
    delnodes::Vector{Int64}
    taxaname::Vector{String}

    function Tree()
        new(-9, 0, Node[], false, Int64[], String[])
    end

    function Tree(ntaxa::Int64, rooted::Bool)
        nodes=Node[]
        rootnode = ntaxa+1
        if rooted
        nnodes = 2*ntaxa-1
        else 
        nnodes = 2*ntaxa-2
        end

        if rooted
        for i in 1:nnodes
            push!(nodes, Node())
        end
        else
        for i in 1:nnodes
            if i == ntaxa+1
            push!(nodes, Node(true))
            else
            push!(nodes, Node())
            end
        end
        end
        taxaname = fill("", ntaxa)
        new(rootnode, ntaxa, nodes, rooted, Int64[], taxaname)
    end

        function Tree(rootnode::Int64, ntaxa::Int64, nodes::Array{Node,1}, rooted::Bool)
        new(rootnode, ntaxa, nodes, rooted, Int64[], String[])
    end
end

mutable struct ModelParameters
    nstates::Int64
    states::Vector{Char}
    conditionprob::Matrix{Float64}
    qmatrix::Matrix{Float64}
    qmatrix_eigvecs::Matrix{Float64}
    qmatrix_eigvecs_inv::Matrix{Float64}
    qmatrix_eigvals::Vector{Float64}
    basefreq::Vector{Float64}
  
    function ModelParameters(tree::Tree, states::Vector{Char})
        nstates = length(states)
        if tree.rooted
            new(nstates, states, fillzeros(2*tree.ntaxa-1,nstates), zeros(nstates,nstates), zeros(nstates,nstates), zeros(nstates,nstates), zeros(nstates), zeros(nstates))
        else 
            new(nstates, states, zeros(2*tree.ntaxa-2,nstates), zeros(nstates,nstates), zeros(nstates,nstates), zeros(nstates,nstates), zeros(nstates), zeros(nstates))
        end
    end
  
    function ModelParameters(states::Vector{Char}, conditionprob::Matrix{Float64})
        nstates = length(states)
        new(nstates, states, conditionprob, zeros(nstates,nstates), zeros(nstates,nstates), zeros(nstates,nstates), zeros(nstates), zeros(nstates))
    end
  
    function ModelParameters(states::Vector{Char}, conditionprob::Matrix{Float64}, qmatrix::Matrix{Float64}, basefreq::Vector{Float64})
        nstates = length(states)
        qmatrix_eigvals = exp.(eigvals(qmatrix))
        qmatrix_eigvecs = eigvecs(qmatrix)
        qmatrix_eigvecs_inv = inv(qmatrix_eigvecs) 
        
        if sum(basefreq) != 1.0
            error("the sum of base frequencies is not equal to 1")
        end
        if length(basefreq) != nstates
            error("the number of states is not equal to the length of basefreq")
        end

        new(nstates, states, conditionprob, qmatrix, qmatrix_eigvecs, qmatrix_eigvecs_inv, qmatrix_eigvals, basefreq)
    end
end
  
  
