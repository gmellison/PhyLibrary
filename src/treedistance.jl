function findAllPathDistance(tree::Tree)
    x = findAllAncestorDistance(tree)
    result = zeros(tree.ntaxa, tree.ntaxa)
    for i in 1:(tree.ntaxa-1)
      for j in (i+1):tree.ntaxa
        anci = get(x[i], "ancestor", 3)
        ancj = get(x[j], "ancestor", 3)
        disi = get(x[i], "distance", 3)
        disj = get(x[j], "distance", 3)
        d1 = 0.0
        d2 = 0.0
        for k in 1:length(anci)
          index = (anci[k] .== ancj)
          if sum(index) > 0
              d1 = disi[k]
              d2 = disj[index]
              break
          end
        end
        result[i,j] = round(d1+d2[1],digits=10)
      end
    end
    return result
end
  
function findPathDistance(tree::Tree, inode::Int64, jnode::Int64)
    d1 = -1.0
    d2 = -1.0
    if inode > tree.ntaxa || jnode > tree.ntaxa
        error("this function calculates distance between leaves")
    end
  
    x1 = findAncestorDistance(tree,inode)
    x2 = findAncestorDistance(tree,jnode)
  
    anc1 = get(x1, "ancestor", 3)
    anc2 = get(x2, "ancestor", 3)
    dis1 = get(x1, "distance", 3)
    dis2 = get(x2, "distance", 3)
  
    for i in 1:length(anc1)
        index = (anc1[i] .== anc2)
        if sum(index) > 0
            d1 = dis1[i]
            d2 = dis2[index]
            break
        end
    end
  
    if d1 == -1.0 || d2 == -1.0
        error("could not find the common ancestor of two nodes")
    end
  
    return round(d1+d2[1],digits=10)
end
  
function findAncestorDistance(tree::Tree, inode::Int64)
    anc = Int64[]
    distance = Float64[]
  
    #initialize anc and distance
    push!(anc,inode)
    push!(distance,0.0)
  
    #update anc and distance
    father = inode
    dis = 0.0
    while father != tree.rootnode
        if tree.nodes[father].brlens > 0.0
            dis = dis + tree.nodes[father].brlens
        end
        father = tree.nodes[father].father
        push!(anc,father)
        push!(distance,dis)
    end
    result = Dict([("ancestor", anc), ("distance", round.(distance,digits=10))])
    return result
end
  
  
function findAllAncestorDistance(tree::Tree)
    result = Dict[]
    for i in 1:tree.ntaxa
      x = findAncestorDistance(tree, i)
      push!(result, x)
    end
    return result
end
  
function treeSinglePartition(tree::Tree, inode::Int64)
    if inode == tree.rootnode
        error("the node cannot be the root node")
    end
  
    candidates = fill(false, tree.ntaxa)
    leaves = Int64[]
    candidates[findDownLeaves(tree, inode, leaves)] .= true
  
    return candidates
  end
  
function treeAllPartition(tree::Tree)
    result = fill(false, (2*tree.ntaxa-3, tree.ntaxa))
    if tree.rooted
        candidates = collect(1:(2*tree.ntaxa-1))
        delnode = cat(tree.rootnode,tree.nodes[tree.rootnode].sons[1],dims=1)
        deleteat!(candidates, sort(delnode))
        for i in 1:length(candidates)
            result[i,:] = treeSinglePartition(tree, candidates[i])
        end
    else
        candidates = collect(1:(2*tree.ntaxa-2))
        deleteat!(candidates, tree.rootnode)
        for i in 1:length(candidates)
            result[i,:] = treeSinglePartition(tree, candidates[i])
        end
    end
    return result
end
        
  
function treeRFdistance(tree1::Tree, tree2::Tree)
    partition1 = treeAllPartition(tree1)
    partition1 = partition1[:, sortperm(tree1.taxaname)]
    partition2 = treeAllPartition(tree2) 
    partition2 = partition2[:, sortperm(tree2.taxaname)]
  
    distance = 2*(tree1.ntaxa - 3)
    for i in (tree1.ntaxa+1):(2*tree1.ntaxa-3)
        for j in (tree2.ntaxa+1):(2*tree2.ntaxa-3)
            if sum(partition1[i,:]+partition2[j,:] .== 1) == 0 || sum(partition1[i,:]+partition2[j,:] .== 1) == tree1.ntaxa
                distance -= 2
                break
            end
        end
    end
    return distance
  end
  
function findAllTriple(tree::Tree)
    anc = findAllAncestorDistance(tree)
    orderedtaxa = sort(tree.taxaname)
    ntriple = Int64(tree.ntaxa*(tree.ntaxa-1)*(tree.ntaxa-2)/6)
    triple = fill("", ntriple)
    index = 1
    for i in 1:(tree.ntaxa-2)
        for j in (i+1):(tree.ntaxa-1)
            for k in (j+1):(tree.ntaxa)
                triple[index] = findSingleTriple(tree, orderedtaxa[i], orderedtaxa[j], orderedtaxa[k], anc)
                index += 1
            end
        end
    end
    return triple
end
  
function findSingleTriple(tree::Tree, taxon1::String, taxon2::String, taxon3::String, ancestors::Vector{Dict})
    if !tree.rooted
      error("It must be a rooted tree")
    end
  
    orderedtaxa = sort([taxon1,taxon2,taxon3])
    leaf1 = ((1:tree.ntaxa)[orderedtaxa[1] .== tree.taxaname])[1]
    leaf2 = ((1:tree.ntaxa)[orderedtaxa[2] .== tree.taxaname])[1]
    leaf3 = ((1:tree.ntaxa)[orderedtaxa[3] .== tree.taxaname])[1]
    x1 = get!(ancestors[leaf1],"ancestor",3)
    x2 = get!(ancestors[leaf2],"ancestor",3)
    x3 = get!(ancestors[leaf3],"ancestor",3)
  
    #find the position of the common ancestor
    anc1 = findCommonAncestor(x1, x2)[2]
    anc2 = findCommonAncestor(x1, x3)[2]
  
    #three triples are determined by two positions anc1 and anc2.
    if anc1 < anc2
        return string("((", orderedtaxa[1], ",", orderedtaxa[2], ")", ",", orderedtaxa[3],")")
    elseif anc1 > anc2
        return string("((", orderedtaxa[1], ",", orderedtaxa[3], ")", ",", orderedtaxa[2],")")
    else
        return string("((", orderedtaxa[2], ",", orderedtaxa[3], ")", ",", orderedtaxa[1],")")
    end
end
  
function findCommonAncestor(anc1::Vector{Int64}, anc2::Vector{Int64})
    result = 0
    for i in 1:length(anc1)
        index = (anc1[i] .== anc2)
        if sum(index) > 0
            result = i
            break
        end
    end
    return cat(anc1[result], result, dims=1)
end
 
function tripleFrequancy(treefile, taxaname::Vector{String})
    tree = readTree(treefile)
    ntree = length(tree)
    ntaxa = length(taxaname)
    if ntaxa < 3
        error("triple involves at least 3 species")
    end
    ntriples = Int64(ntaxa*(ntaxa-1)*(ntaxa-2)/6)

    taxontriple = fill("",(ntriples,3))
    index = 1
    for i in 1:(ntaxa-2)
        for j in (i+1):(ntaxa-1)
            for k in (j+1):ntaxa
                taxontriple[index,:] = sort([taxaname[i], taxaname[j], taxaname[k]])
                index += 1
            end
        end
    end

    triple = fill("", ntree, ntriples)
    for w in 1:ntree
        anc = findAllAncestorDistance(tree[w])
        index = 1
        for i in 1:(ntaxa-2)
            for j in (i+1):(ntaxa-1)
                for k in (j+1):ntaxa
                    if sum(taxaname[i] .== tree[w].taxaname) == 0 || sum(taxaname[j] .== tree[w].taxaname) == 0 || sum(taxaname[k] .== tree[w].taxaname) == 0 
                        index += 1
                    else
                        triple[w,index] = findSingleTriple(tree[w], taxaname[i], taxaname[j], taxaname[k], anc)
                        index += 1
                    end
                end
            end
        end
    end

    result = fill(1,(ntriples,4))
    for i in 1:ntriples
        x1 = string("((", taxontriple[i,1], ",", taxontriple[i,2], ")", ",", taxontriple[i,3],")")
        x2 = string("((", taxontriple[i,1], ",", taxontriple[i,3], ")", ",", taxontriple[i,2],")")
        x3 = string("((", taxontriple[i,2], ",", taxontriple[i,3], ")", ",", taxontriple[i,1],")")
        result[i,1] = sum("" .== triple[:,i])
        result[i,2] = sum(x1 .== triple[:,i])
        result[i,3] = sum(x2 .== triple[:,i])
        result[i,4] = sum(x3 .== triple[:,i])
    end

    triplefreq = Tables.table(hcat(taxontriple, result),header = ["taxon1","taxon2","taxon3","missing","triple 12|3","triple 13|2","triple 23|1"])
    return triplefreq
end

function treeTripleDistance(tree1::Tree, tree2::Tree)
    triple1 = findAllTriple(tree1)
    triple2 = findAllTriple(tree2)
  
    distance = 0
    for i in 1:length(triple1)
        if sum(triple1[i] .== triple2) == 0
            distance += 1
        end
    end
    return distance
end
  
function findSingleQuartet(tree::Tree, taxon1::String, taxon2::String, taxon3::String, taxon4::String, ancestors::Vector{Dict})
    if tree.rooted
      error("It must be an unrooted tree")
    end
  
    orderedtaxa = sort([taxon1,taxon2,taxon3,taxon4])
    leaf1 = ((1:tree.ntaxa)[orderedtaxa[1] .== tree.taxaname])[1]
    leaf2 = ((1:tree.ntaxa)[orderedtaxa[2] .== tree.taxaname])[1]
    leaf3 = ((1:tree.ntaxa)[orderedtaxa[3] .== tree.taxaname])[1]
    leaf4 = ((1:tree.ntaxa)[orderedtaxa[4] .== tree.taxaname])[1]
    x1 = get!(ancestors[leaf1],"ancestor",3)
    x2 = get!(ancestors[leaf2],"ancestor",3)
    x3 = get!(ancestors[leaf3],"ancestor",3)
    x4 = get!(ancestors[leaf4],"ancestor",3)
  
    #find the position of the common ancestor
    anc1 = findCommonAncestor(x1, x2)[2]
    anc2 = findCommonAncestor(x1, x3)[2]
    anc3 = findCommonAncestor(x1, x4)[2]
    anc4 = findCommonAncestor(x2, x3)[2]
    anc5 = findCommonAncestor(x2, x4)[2]
  
    #three triples are determined by two positions anc1 and anc2.
    if anc1 < anc2
        return string("((", orderedtaxa[1], ",", orderedtaxa[2], ")", ",(", orderedtaxa[3],",", orderedtaxa[4], "))")
    elseif anc2 < anc3
        return string("((", orderedtaxa[1], ",", orderedtaxa[3], ")", ",(", orderedtaxa[2],",", orderedtaxa[4], "))")
    elseif anc3 < anc1
        return string("((", orderedtaxa[1], ",", orderedtaxa[4], ")", ",(", orderedtaxa[2],",", orderedtaxa[3], "))")
    elseif anc4 < anc5
        return string("((", orderedtaxa[1], ",", orderedtaxa[4], ")", ",(", orderedtaxa[2],",", orderedtaxa[3], "))")
    elseif anc5 < anc4
        return string("((", orderedtaxa[1], ",", orderedtaxa[3], ")", ",(", orderedtaxa[2],",", orderedtaxa[4], "))")
    else
        error("this is a polynomy node")
    end
end
  
function findAllQuartet(tree::Tree)
    anc = findAllAncestorDistance(tree)
    orderedtaxa = sort(tree.taxaname)
    nquartet = Int64(tree.ntaxa*(tree.ntaxa-1)*(tree.ntaxa-2)*(tree.ntaxa-3)/24)
    quartet = fill("", nquartet)
    index = 1
    for i in 1:(tree.ntaxa-3)
        for j in (i+1):(tree.ntaxa-2)
            for k in (j+1):(tree.ntaxa-1)
                for l in (k+1):(tree.ntaxa)
                    quartet[index] = findSingleQuartet(tree, orderedtaxa[i], orderedtaxa[j], orderedtaxa[k], orderedtaxa[l], anc)
                    index += 1
                end
            end
        end
    end
    return quartet
end
  
function treeQuartetDistance(tree1::Tree, tree2::Tree)
    quartet1 = findAllQuartet(tree1)
    quartet2 = findAllQuartet(tree2)
  
    distance = 0
    for i in 1:length(quartet1)
        if sum(quartet1[i] .== quartet2) == 0
            distance += 1
        end
    end
    return distance
end

function treeSpeciesFrequency(tree::Vector{Tree})
    spname=Vector{String}
    for i in 1:length(tree)
        if i == 1
            spname = copy(tree[i].taxaname)
        else
            spname = vcat(spname,tree[i].taxaname)
        end
    end
    result = freqtable(spname)
    return result
end

function treeDistribution(tree::Vector{Tree})
    nalltree = length(tree)
    ntree = length(tree)
    treestr = String[]
    treeprob = Float64[]
    while ntree > 0
        push!(treestr, writeTree(tree[1],"text",true))
        index = fill(false,ntree)
        for j in 2:ntree
            if treeRFdistance(tree[1],tree[j]) > 0
                index[j] = true
            end
        end
        tree = tree[index]
        ntree = sum(index)
        push!(treeprob, (length(index)-ntree)/nalltree)
    end
    return treestr, treeprob
end