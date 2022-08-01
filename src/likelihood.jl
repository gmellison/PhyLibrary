function preorderTraversal(inode::Int64, tree::Tree, traversal::Vector{Int64})
    push!(traversal, inode)
    if tree.nodes[inode].nson > 0
        for i in 1:tree.nodes[inode].nson
            son = tree.nodes[inode].sons[i]
            preorderTraversal(son, tree, traversal)
        end
    end
end
  
function conditionProbAllNode(tree::Tree, traversal::Vector{Int64}, modelpara::ModelParameters)
    for i in 1:length(traversal)
        if traversal[i] == tree.rootnode
            son1 = tree.nodes[traversal[i]].sons[1]
            son2 = tree.nodes[traversal[i]].sons[2]
            son3 = tree.nodes[traversal[i]].sons[3]
            son1prob = modelpara.conditionprob[son1,:]
            son2prob = modelpara.conditionprob[son2,:]
            son3prob = modelpara.conditionprob[son3,:]
            brlens1 = tree.nodes[son1].brlens
            brlens2 = tree.nodes[son2].brlens
            brlens3 = tree.nodes[son3].brlens
   
            pmatrix1 = modelpara.qmatrix_eigvecs * Diagonal(modelpara.qmatrix_eigvals.^brlens1) * modelpara.qmatrix_eigvecs_inv
            pmatrix2 = modelpara.qmatrix_eigvecs * Diagonal(modelpara.qmatrix_eigvals.^brlens2) * modelpara.qmatrix_eigvecs_inv
            pmatrix3 = modelpara.qmatrix_eigvecs * Diagonal(modelpara.qmatrix_eigvals.^brlens3) * modelpara.qmatrix_eigvecs_inv
    
            nstates = length(son1prob)
            prob = zeros(nstates)
            for j in 1:nstates
                for k in 1:nstates
                    for l in 1:nstates
                        for m in 1:nstates
                            prob[j] += pmatrix1[j,k]*son1prob[k]*pmatrix2[j,l]*son2prob[l]*pmatrix3[j,m]*son3prob[m]
                        end
                    end
                end
            end
            modelpara.conditionprob[traversal[i],:] = prob
        elseif traversal[i] > tree.ntaxa
            son1 = tree.nodes[traversal[i]].sons[1]
            son2 = tree.nodes[traversal[i]].sons[2]
            son1prob = modelpara.conditionprob[son1,:]
            son2prob = modelpara.conditionprob[son2,:]
            brlens1 = tree.nodes[son1].brlens
            brlens2 = tree.nodes[son2].brlens
  
            pmatrix1 = modelpara.qmatrix_eigvecs * Diagonal(modelpara.qmatrix_eigvals.^brlens1) * modelpara.qmatrix_eigvecs_inv
            pmatrix2 = modelpara.qmatrix_eigvecs * Diagonal(modelpara.qmatrix_eigvals.^brlens2) * modelpara.qmatrix_eigvecs_inv
  
            modelpara.conditionprob[traversal[i],:] = conditionProbSingleNode(son1prob, son2prob, pmatrix1, pmatrix2)
        end
    end
end
    
    
function conditionProbSingleNode(son1prob::Vector{Float64}, son2prob::Vector{Float64}, pmatrix1::Matrix{Float64}, pmatrix2::Matrix{Float64})
    nstates = length(son1prob)
    prob = zeros(nstates)
    for i in 1:nstates
        for j in 1:nstates
            for k in 1:nstates
                prob[i] += pmatrix1[i,j]*son1prob[j]*pmatrix2[i,k]*son2prob[k]
            end
        end
    end
    return prob
end
  
function LoglikelihoodDNA(patternfreq::Vector{Int64}, pattern::Vector{String}, tree::Tree, modelpara::ModelParameters)
    loglike = 0.0
    for i in 1:length(pattern)
      loglike += (patternfreq[i]*LoglikelihoodSingleDNA(pattern[i], tree, modelpara))
    end
    return loglike
end
  
function LoglikelihoodSingleDNA(nuc::String, tree::Tree, modelpara::ModelParameters)
    for i in 1:tree.ntaxa
        modelpara.conditionprob[i,:] = zeros(modelpara.nstates)
        for j in 1:modelpara.nstates
            if nuc[i] == modelpara.states[j]
                modelpara.conditionprob[i,j] = 1
                break
            end
        end
    end
  
    traversal=Int64[]
    preorderTraversal(tree.rootnode, tree, traversal)
    conditionProbAllNode(tree, reverse(traversal), modelpara)
    return log(transpose(modelpara.conditionprob[tree.rootnode,:]) * modelpara.basefreq)
end

function LoglikelihoodDNApairwise(doubletfreq::Matrix{Int64}, tree::Tree, modelpara::ModelParameters)
    distance = findAllPathDistance(tree)
    dis = distance[upper_tri(distance)] 
    npair = length(dis)
  
    loglike = 0.0
    for i in 1:npair
        brlens = dis[i]
        #pmatrix = exp(qmatrix*brlens)
        pmatrix = modelpara.qmatrix_eigvecs * Diagonal(modelpara.qmatrix_eigvals.^brlens) * modelpara.qmatrix_eigvecs_inv
        for j in 1:4
            pmatrix[j,:] = pmatrix[j,:] * modelpara.basefreq[j]
        end
        logprob = log.(vec(pmatrix))
        loglike += transpose(doubletfreq[i,:]) * logprob
    end
    return loglike
end
        
  