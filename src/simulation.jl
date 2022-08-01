
##################################################
# simulation
####################################################
function simNucleotideSingleBranch(states::Vector{Char}, seqlength::Int64, basefreq::Vector{Float64}, rates::Vector{Float64}, brlens::Float64)
    nstates = length(states)
    if length(basefreq) != nstates
        error("the length of basefreq != the number of states")
    end

    if sum(basefreq) != 1.0
        error("the sum of basefreq must be equal to 1.0")
    end

    #create the Q matrix
    qmatrix = submodelGTR_rateMatrix(basefreq, rates)
    
    #find the transition probability matrix
    pmatrix = exp(qmatrix*brlens)
  
    #find the probabilities of 16 doublets
    doubleP = zeros(nstates,nstates)
    for i in 1:nstates 
      doubleP[i,:] = basefreq[i]*pmatrix[i,:]
    end
  
    #generate random nucleotides
    x = fill(0,(2,seqlength))
    x[1,:] = sample(1:nstates, Weights(basefreq), seqlength)
    for i in 1:seqlength 
      x[2,i] = sample(1:nstates, Weights(pmatrix[x[1,i],:]))
    end
  
    data = repeat([""],2)
    data[1] = join(transpose(x[1,:]))
    data[2] = join(transpose(x[2,:]))
  
    #data = replace.(y, r"[ \]\[]"=>"")
    for i in 1:nstates
      data = replace.(data,"$i"=>states[i])
    end
  
    return data
end
  
function simNucleotide(states::Vector{Char}, tree::Tree, basefreq::Vector{Float64}, rates::Vector{Float64}, seqlength::Int64)
    nstates = length(states)
    if length(basefreq) != nstates
        error("the length of basefreq != the number of states")
    end

    if sum(basefreq) != 1.0
        error("the sum of basefreq must be equal to 1.0")
    end

    qmatrix = submodelGTR_rateMatrix(basefreq, rates)
    traversal=Int64[]
    preorderTraversal(tree.rootnode, tree, traversal)
  
    for i in 1:length(traversal)
        if traversal[i] == tree.rootnode
            tree.nodes[traversal[i]].label = join(sample(1:nstates, Weights(basefreq), seqlength))
        else
            brlens = tree.nodes[traversal[i]].brlens
            pmatrix = exp(qmatrix*brlens)
            father = tree.nodes[traversal[i]].father
            seq = parse.(Int64, split(tree.nodes[father].label,""))
            x = collect(1:seqlength)
            for j in 1:length(seq)
                x[j] = sample(collect(1:nstates), Weights(pmatrix[seq[j],:]),1)[1]
            end
            tree.nodes[traversal[i]].label = join(x)
        end
    end
  
    sequence = fill("",tree.ntaxa)
    for i in 1:tree.ntaxa
        x = tree.nodes[i].label
        for j in 1:nstates
            x = replace(x,"$j"=>states[j])
        end
        sequence[i] = x
    end
    return SEQdata(tree.taxaname, sequence)
end
  