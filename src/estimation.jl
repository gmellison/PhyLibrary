function submodelGTR_rateMatrix_estimation(data::Vector{String}, states::Vector{Char})
  nstates = length(states)
  basefreq = sequence_base_freq(data, states)
  doublefreq = sequence_doublet_freq(data, states)
  doublefreq = reshape(doublefreq, nstates, nstates)
  doublefreq = (doublefreq + transpose(doublefreq))/2

  pmatrix_hat = zeros(nstates, nstates)
  for i in 1:nstates 
    pmatrix_hat[i,:] = doublefreq[i,:]/basefreq[i]
  end

  eig = eigen(pmatrix_hat)
  eva = zeros(nstates, nstates)
  eva[diag(eva)] = log.(eig.values)
  eve = eig.vectors
  qt = eve * eva * inv(eve)

  brlens = - transpose(basefreq) * qt[diag(qt)]
  qmatrix_hat = qt/brlens
  qmatrix_hat = qmatrix_hat * diag(1 ./ basefreq)

  #making the last rate [nstates-1, nstate] = 1.0
  qmatrix_hat = qmatrix_hat/qmatrix_hat[nstates-1, nstates]

  rates_hat = qmatrix_hat[upper_tri(qmatrix_hat)]
  result = Dict([("basefreq", basefreq),("rates", abs.(rates_hat)),("brlens",brlens)])

  return result
end

function treeNewtonRaphsonSingleBrlensFull(tree::Tree, inode::Int64, patternfreq::Vector{Int64}, pattern::Vector{String}, modelpara::ModelParameters)
  minimumbrlens = 0.0001
  a1 = tree.nodes[inode].brlens
  delta = 0.00001
  a2 = a1+delta
  a3 = a1-delta
  diff = 1.0

  while diff > ERROR
      tree.nodes[inode].brlens = a1
      x = LoglikelihoodDNA(patternfreq, pattern, tree, modelpara)
      tree.nodes[inode].brlens = a2
      y = LoglikelihoodDNA(patternfreq, pattern, tree, modelpara)
      tree.nodes[inode].brlens = a3
      z = LoglikelihoodDNA(patternfreq, pattern, tree, modelpara)

      der1 = (y-x)/delta
      der2 = (y-2*x+z)/(delta^2)

      if abs(der2) < ERROR
        println("der2 is too small\n\n")
        break
      else
        w = der1/der2
      end

      #update brlens a1, a2, a3
      a1 = a1 - w
      if a1 < minimumbrlens
        a1 = minimumbrlens
        break
      end
      a2 = a1+delta
      a3 = a1-delta
      diff = abs(w)
  end
end

function treeNewtonRaphsonSingleBrlensPairwise(tree::Tree, inode::Int64, doubletfreq::Matrix{Int64}, modelpara::ModelParameters)
  minimumbrlens = 0.0001
  a1 = tree.nodes[inode].brlens
  delta = 0.00001
  a2 = a1+delta
  a3 = a1-delta
  diff = 1.0

  while diff > ERROR
      tree.nodes[inode].brlens = a1
      x = LoglikelihoodDNApairwise(doubletfreq, tree, modelpara)
      tree.nodes[inode].brlens = a2
      y = LoglikelihoodDNApairwise(doubletfreq, tree, modelpara)
      tree.nodes[inode].brlens = a3
      z = LoglikelihoodDNApairwise(doubletfreq, tree, modelpara)

      der1 = (y-x)/delta
      der2 = (y-2*x+z)/(delta^2)

      if abs(der2) < ERROR
        println("der2 is too small\n\n")
        break
      else
        w = der1/der2
      end

      #update brlens a1, a2, a3
      a1 = a1 - w
      if a1 < minimumbrlens
        a1 = minimumbrlens
        break
      end
      a2 = a1+delta
      a3 = a1-delta
      diff = abs(w)
  end
end


function treeNewtonRaphsonBrlens(tree::Tree, patternfreq::Vector{Int64}, pattern::Vector{String}, modelpara::ModelParameters)
  diff = 1.0
  nnodes = length(tree.nodes)
  threshold = ERROR * nnodes
  niteration = 1
  while diff > threshold && niteration < 100
      diff = 0.0
      for i in 1:nnodes
          x = tree.nodes[i].brlens
          if x >= 0.0
              treeNewtonRaphsonSingleBrlensFull(tree, i, patternfreq, pattern, modelpara)
              diff += abs(tree.nodes[i].brlens-x)
          end
      end
      niteration += 1
  end
end
 
function treeMaximumLikelihood(data::SEQdata, states::Vector{Char}, starttree::Tree, likelihood="full") 
  pattern_freq = alignmentPattern(data.sequence)
  pattern = names(pattern_freq)[1]
  patternfreq = Vector(pattern_freq)

  doublet = sequence_doublet_freq_pair(data,states)
  doubletfreq = doublet["doubletfreq"]
  
  rate_estimates = submodelGTR_rateMatrix_estimation(data.sequence, states)
  rates = rate_estimates["rates"]
  basefreq = rate_estimates["basefreq"] 
  basefreq[length(states)] = 1.0 - sum(basefreq[1:(length(states)-1)])
  println("The base frequencies: ", basefreq, "\nThe sum is ", sum(basefreq))
  qmatrix = submodelGTR_rateMatrix(basefreq, rates)

  #starttree = treeGenerator(data.taxaname,"")
  print("\n\n\n       Julia: Maximum Likelihood Estimation of Phylogenetic Trees\n\n                          Liang Liu\n\n")
  println("       The initial tree:  ", writeTree(starttree,"text",true),"\n")
  modelpara = ModelParameters(states,zeros(2*starttree.ntaxa-2,4),qmatrix,basefreq)

  if likelihood == "full"
    oldloglike = LoglikelihoodDNA(patternfreq, pattern, starttree, modelpara)
  else
    oldloglike = LoglikelihoodDNApairwise(doubletfreq, starttree, modelpara)
  end

  niteration = 0
  oldtree = deepcopy(starttree)
  newtree = deepcopy(starttree)
  nomove = 0

  while nomove < 20 && niteration < 10000
      println("       ", niteration,"   ------  ", oldloglike)
      treeNNI(newtree)
      if likelihood == "full"
        newloglike = LoglikelihoodDNA(patternfreq, pattern, newtree, modelpara)
      else
        newloglike = LoglikelihoodDNApairwise(doubletfreq, newtree, modelpara)
      end
      #newloglike = LoglikelihoodDNA(patternfreq, pattern, newtree, modelpara)
      #@time newloglike = LoglikelihoodDNApairwise(doubletfreq, newtree, qmatrix, basefreq)
      if newloglike > oldloglike
          oldtree = deepcopy(newtree)
          oldloglike = newloglike
          nomove = 0
      else
          newtree = deepcopy(oldtree)
          nomove += 1
      end
      niteration += 1

      if mod(niteration, 100) == 0
        println("\n\n       Updating the branch length ...")
        @time for i in 1:(2*newtree.ntaxa-2)
          x = newtree.nodes[i].brlens
          if x >= 0.0
            if likelihood == "full"
              treeNewtonRaphsonSingleBrlensFull(newtree, i, patternfreq, pattern, modelpara)
            else
              treeNewtonRaphsonSingleBrlensPairwise(newtree, i, doubletfreq, modelpara)
            end
          end
        end
      end
  end

  println("\n\n       Finalize branch length estimation")
  for i in 1:(2*newtree.ntaxa-2)
    println("       Updating the branch length ", i)
    x = newtree.nodes[i].brlens
    if x >= 0.0
      if likelihood == "full"
        treeNewtonRaphsonSingleBrlensFull(newtree, i, patternfreq, pattern, modelpara)
      else
        treeNewtonRaphsonSingleBrlensPairwise(newtree, i, doubletfreq, modelpara)
      end
    end
  end

  if likelihood == "full"
    maxloglike = LoglikelihoodDNA(patternfreq, pattern, newtree, modelpara)
  else
    maxloglike = LoglikelihoodDNApairwise(doubletfreq, newtree, modelpara)
  end

  println("\n       The maximum loglikelihood score: ", maxloglike)
  println("\n       The final tree:", writeTree(newtree,"text",true))

  return newtree, maxloglike

end

  
  