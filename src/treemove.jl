function findDownLeaves(tree::Tree, inode::Int64, leaves::Vector{Int64})
    if inode <= tree.ntaxa
      push!(leaves,inode)
    else
      if inode == tree.rootnode && !tree.rooted
        son1 = tree.nodes[inode].sons[1]
        son2 = tree.nodes[inode].sons[2]
        son3 = tree.nodes[inode].sons[3]
        findDownLeaves(tree,son1,leaves)
        findDownLeaves(tree,son2,leaves)
        findDownLeaves(tree,son3,leaves)
      else
        son1 = tree.nodes[inode].sons[1]
        son2 = tree.nodes[inode].sons[2]
        findDownLeaves(tree,son1,leaves)
        findDownLeaves(tree,son2,leaves)
      end
    end
end
  
function findDownNodes(tree::Tree, inode::Int64, nodes::Vector{Int64})
    push!(nodes,inode)
    if inode > tree.ntaxa
      if inode == tree.rootnode && !tree.rooted
        son1 = tree.nodes[inode].sons[1]
        son2 = tree.nodes[inode].sons[2]
        son3 = tree.nodes[inode].sons[3]
        findDownNodes(tree,son1,nodes)
        findDownNodes(tree,son2,nodes)
        findDownNodes(tree,son3,nodes)
      else
        son1 = tree.nodes[inode].sons[1]
        son2 = tree.nodes[inode].sons[2]
        findDownNodes(tree,son1,nodes)
        findDownNodes(tree,son2,nodes)
      end
    end
end
  
function deleteNode(tree::Tree, inode::Int64)
      if inode == tree.rootnode 
      error("cannot delete the tree root")
    end
  
    if sum(inode .== tree.delnodes) > 0
      error("the node has already been deleted")
    end
  
    if tree.ntaxa-sum(tree.delnodes .< tree.ntaxa) < 4
      error("cannot delete nodes if # of taxa < 4")
    end
  
      father = tree.nodes[inode].father
  
      if father == tree.rootnode
      if tree.rooted
        if tree.nodes[father].sons[1] == inode 
          brother = tree.nodes[father].sons[2]
        else 
          brother = tree.nodes[father].sons[1]
        end
        tree.rootnode = brother
      else
        if tree.nodes[father].sons[1] == inode 
          if tree.nodes[father].sons[2] > tree.ntaxa
            rootnode = tree.nodes[father].sons[2]
            brother = tree.nodes[father].sons[3]
          elseif tree.nodes[father].sons[3] > tree.ntaxa
            rootnode = tree.nodes[father].sons[3]
            brother = tree.nodes[father].sons[2]
          else
            #if both sons are leaves, this is a rooted tree of 2 taxa
            deleteat!(tree.nodes[father].sons,1)
            tree.nodes[father].nson = 2
            tree.rooted = true
          end
        elseif tree.nodes[father].sons[2] == inode
          if tree.nodes[father].sons[1] > tree.ntaxa
            rootnode = tree.nodes[father].sons[1]
            brother = tree.nodes[father].sons[3]
          elseif tree.nodes[father].sons[3] > tree.ntaxa
            rootnode = tree.nodes[father].sons[3]
            brother = tree.nodes[father].sons[1]
          else
            #if both sons are leaves, this is a rooted tree of 2 taxa
            deleteat!(tree.nodes[father].sons,2)
            tree.nodes[father].nson = 2
            tree.rooted = true
          end
        else 
          if tree.nodes[father].sons[1] > tree.ntaxa
            rootnode = tree.nodes[father].sons[1]
            brother = tree.nodes[father].sons[2]
          elseif tree.nodes[father].sons[2] > tree.ntaxa
            rootnode = tree.nodes[father].sons[2]
            brother = tree.nodes[father].sons[1]
          else
            #if both sons are leaves, this is a rooted tree of 2 taxa
            deleteat!(tree.nodes[father].sons,3)
            tree.nodes[father].nson = 2
            tree.rooted = true
          end  
        end
        if !tree.rooted
          tree.nodes[father].nson = 2
          deleteat!(tree.nodes[father].sons,3)
          tree.rootnode = rootnode
          tree.nodes[tree.rootnode].nson = 3
          push!(tree.nodes[rootnode].sons, brother)
          tree.nodes[brother].father = rootnode
        end
      end
      else
      if tree.nodes[father].sons[1] == inode 
        brother = tree.nodes[father].sons[2]
      else 
        brother = tree.nodes[father].sons[1]
      end
          grandfather = tree.nodes[father].father
          if tree.nodes[grandfather].sons[1] == father
              tree.nodes[grandfather].sons[1] = brother
          else
              tree.nodes[grandfather].sons[2] = brother
      end
          tree.nodes[brother].father = grandfather
    end
    tree.nodes[tree.rootnode].father = -9
    nodes = Int64[]
    findDownNodes(tree, inode, nodes)
    tree.delnodes = cat(tree.delnodes,father, nodes,dims=1)
  
end  
  
  
function addNode(tree::Tree, fromnode::Int64, tonode::Int64)
    fromfather = tree.nodes[fromnode].father
    if sum(tree.delnodes .== fromnode) + sum(tree.delnodes .== fromfather) < 2
      error("the added node is already in the tree")
    end
    if sum(tree.delnodes .== tonode) > 0 
      error("tonode has been deleted")
    end
  
    #update tree.delnodes
    index = fill(true,length(tree.delnodes))
    nodes = Int64[]
    findDownNodes(tree, Int64(fromnode), nodes)
    push!(nodes,fromfather)
    for i in 1:length(tree.delnodes)
      if sum(nodes .== tree.delnodes[i]) > 0
        index[i] = false
      end
    end
    tree.delnodes = tree.delnodes[index]
  
    #update ancestor of delnodes
    if tree.nodes[fromfather].sons[1] == fromnode
      brother = tree.nodes[fromfather].sons[2]
    else
      brother = tree.nodes[fromfather].sons[1]
    end
    if sum(tree.delnodes .== brother) > 0
      tree.nodes[brother].father = tree.nodes[fromfather].father
    end
  
    #if fromfather is rootnode, fromnode is attached to rootnode, changing rooted to unrooted 
    if tree.rooted && fromfather == tree.rootnode
      push!(tree.nodes[tree.rootnode].sons, fromnode)
      tree.nodes[tree.rootnode].nson = 3
      tree.rooted = false
      return tree
    end
  
    if tonode == tree.rootnode
      if tree.rooted
        #update rootnode
        tree.rootnode = fromfather
  
        #update ancestral node of tonode
        tree.nodes[tonode].father = fromfather
  
        #update fromfather
        tree.nodes[fromfather].sons[1] = fromnode
        tree.nodes[fromfather].sons[2] = tonode
        tree.nodes[fromfather].father = -9
  
        #update fromnode
        tree.nodes[fromnode].father = fromfather
      else
        #update rootnode
        tree.rootnode = fromfather
  
        #update fromfather
        tree.nodes[fromfather].sons[1] = fromnode
        tree.nodes[fromfather].sons[2] = tonode
        push!(tree.nodes[fromfather].sons, tree.nodes[tonode].sons[1])
        tree.nodes[fromfather].nson = 3
        tree.nodes[fromfather].father = -9
  
        #update tree.nodes[tonode].sons[1]
        tree.nodes[tree.nodes[tonode].sons[1]].father = fromfather
  
        #update ancestral node of tonode
        tree.nodes[tonode].father = fromfather
        tree.nodes[tonode].nson = 2
        deleteat!(tree.nodes[tonode].sons, 1)
  
        #update fromnode
        tree.nodes[fromnode].father = fromfather
      end
    else
      if !tree.rooted && tree.nodes[tonode].father == tree.rootnode
        tofather = tree.nodes[tonode].father
  
        #update the ancestral information of tonode
        tree.nodes[tonode].father = fromfather; 
  
        #update tofather
        if tree.nodes[tofather].sons[1] == tonode
          tree.nodes[tofather].sons[1] = fromfather
        elseif tree.nodes[tofather].sons[2] == tonode
          tree.nodes[tofather].sons[2] = fromfather
        else
          tree.nodes[tofather].sons[3] = fromfather
        end
        
        #update fromnode
        tree.nodes[fromnode].father = fromfather
  
        #update the ancestral information of fromfather
        tree.nodes[fromfather].sons[1] = tonode
        tree.nodes[fromfather].sons[2] = fromnode
        tree.nodes[fromfather].father = tofather
      else
        tofather = tree.nodes[tonode].father
        
        #update the ancestral information of tonode
        tree.nodes[tonode].father = fromfather; 
  
        #update tofather
        if tree.nodes[tofather].sons[1] == tonode
          tree.nodes[tofather].sons[1] = fromfather
        else
          tree.nodes[tofather].sons[2] = fromfather
        end
        
        #update fromnode
        tree.nodes[fromnode].father = fromfather
  
        #update the ancestral information of fromfather
        tree.nodes[fromfather].sons[1] = tonode
        tree.nodes[fromfather].sons[2] = fromnode
        tree.nodes[fromfather].father = tofather
      end
    end	
end
  
function swapNodes(tree::Tree, inode::Int64, jnode::Int64)
    if sum(tree.delnodes .== inode) + sum(tree.delnodes .== jnode) > 0 
      error("at least one node has been deleted")
    end
  
    if inode > tree.ntaxa || jnode > tree.ntaxa
      error("only swap terminal nodes (leaves)")
    end
    
    ifather = tree.nodes[inode].father;
    jfather = tree.nodes[jnode].father;
    tree.nodes[inode].father = jfather;
    tree.nodes[jnode].father = ifather;
  
    if tree.rooted
      if tree.nodes[ifather].sons[1] == inode
        tree.nodes[ifather].sons[1] = jnode
      else
        tree.nodes[ifather].sons[2] = jnode
      end
      if(tree.nodes[jfather].sons[1] == jnode)
        tree.nodes[jfather].sons[1] = inode
      else
        tree.nodes[jfather].sons[2] = inode
      end
    else
      if ifather == tree.rootnode
        if tree.nodes[ifather].sons[1] == inode
          tree.nodes[ifather].sons[1] = jnode
        elseif tree.nodes[ifather].sons[2] == inode
          tree.nodes[ifather].sons[2] = jnode
        else
          tree.nodes[ifather].sons[3] = jnode
        end
        if(tree.nodes[jfather].sons[1] == jnode)
          tree.nodes[jfather].sons[1] = inode
        else
          tree.nodes[jfather].sons[2] = inode
        end
      elseif jfather == tree.rootnode
        if(tree.nodes[jfather].sons[1] == jnode)
          tree.nodes[jfather].sons[1] = inode
        elseif tree.nodes[jfather].sons[2] == jnode
          tree.nodes[jfather].sons[2] = inode
        else
          tree.nodes[jfather].sons[3] = inode
        end
        if tree.nodes[ifather].sons[1] == inode
          tree.nodes[ifather].sons[1] = jnode
        else
          tree.nodes[ifather].sons[2] = jnode
        end
      else
        if tree.nodes[ifather].sons[1] == inode
          tree.nodes[ifather].sons[1] = jnode
        else
          tree.nodes[ifather].sons[2] = jnode
        end
        if(tree.nodes[jfather].sons[1] == jnode)
          tree.nodes[jfather].sons[1] = inode
        else
          tree.nodes[jfather].sons[2] = inode
        end
      end
    end
end
  
  
##################################################
# tree rearrangement NNI
##################################################
function deleteNodeNNI(tree::Tree, inode::Int64)
    father = tree.nodes[inode].father
    grandfather = tree.nodes[father].father
  
    #updtae brother
    if tree.nodes[father].sons[1] == inode 
      brother = tree.nodes[father].sons[2]
    else 
      brother = tree.nodes[father].sons[1]
    end
    tree.nodes[brother].father = grandfather
    
    #update grandfather
    if grandfather == tree.rootnode
      if tree.nodes[grandfather].sons[1] == father
        tree.nodes[grandfather].sons[1] = brother
      elseif tree.nodes[grandfather].sons[2] == father
        tree.nodes[grandfather].sons[2] = brother
      else
        tree.nodes[grandfather].sons[3] = brother
      end
    else
      if tree.nodes[grandfather].sons[1] == father
        tree.nodes[grandfather].sons[1] = brother
      else
        tree.nodes[grandfather].sons[2] = brother
      end
    end
end  
  
function addNodeNNI(tree::Tree, fromnode::Int64, tonode::Int64)
    fromfather = tree.nodes[fromnode].father
    tofather = tree.nodes[tonode].father
  
    if !tree.rooted && tofather == tree.rootnode
      #update tonode
      tree.nodes[tonode].father = fromfather
  
      #update tofather
      if tree.nodes[tofather].sons[1] == tonode
        tree.nodes[tofather].sons[1] = fromfather
      elseif tree.nodes[tofather].sons[2] == tonode
        tree.nodes[tofather].sons[2] = fromfather
      else
        tree.nodes[tofather].sons[3] = fromfather
      end
      
      #update fromfather
      tree.nodes[fromfather].sons[1] = tonode
      tree.nodes[fromfather].sons[2] = fromnode
      tree.nodes[fromfather].father = tofather
    else   
      #update tonode
      tree.nodes[tonode].father = fromfather; 
  
      #update tofather
      if tree.nodes[tofather].sons[1] == tonode
        tree.nodes[tofather].sons[1] = fromfather
      else
        tree.nodes[tofather].sons[2] = fromfather
      end
  
      #update fromfather
      tree.nodes[fromfather].sons[1] = tonode
      tree.nodes[fromfather].sons[2] = fromnode
      tree.nodes[fromfather].father = tofather
    end
end
  
function treeNNI(tree::Tree)
    if tree.rooted
      #select at random an internal branch to perform NNI
      candidates = collect((tree.ntaxa+1):(2*tree.ntaxa-1))
      deleteat!(candidates,[tree.rootnode-tree.ntaxa])  
      if length(candidates) == 1
        internode = candidates
      else
        internode = sample(candidates)
      end
  
      #find the 4 nodes (a1,a2,b1,b2) involved in NNI
      a1 = tree.nodes[internode].sons[1]
      a2 = tree.nodes[internode].sons[2]
  
      father = tree.nodes[internode].father 
      if father == tree.rootnode
        if tree.nodes[father].sons[1] == internode
          b1 = b2 = tree.nodes[father].sons[2]
        else
          b1 = b2 = tree.nodes[father].sons[1]
        end
      else
        b1 = father
        if tree.nodes[father].sons[1] == internode
          b2 = tree.nodes[father].sons[2]
        else
          b2 = tree.nodes[father].sons[1]
        end
      end
      #randomly select two nodes for NNI
      node1 = sample([a1,a2])
      node2 = sample([b1,b2])
    else
      candidates = collect((tree.ntaxa+1):(2*tree.ntaxa-2))
      #delete the rootnode
      deleteat!(candidates,[tree.rootnode-tree.ntaxa])
      #randomly select an internal branch
      if length(candidates) == 1
        internode = candidates
      else
        internode = sample(candidates)
      end
  
      #find the 4 nodes (a1,a2,b1,b2) involved in NNI
      a1 = tree.nodes[internode].sons[1]
      a2 = tree.nodes[internode].sons[2]
  
      father = tree.nodes[internode].father 
      if father == tree.rootnode
        if tree.nodes[father].sons[1] == internode
          b1 = tree.nodes[father].sons[2]
          b2 = tree.nodes[father].sons[3]
        elseif tree.nodes[father].sons[2] == internode
          b1 = tree.nodes[father].sons[1]
          b2 = tree.nodes[father].sons[3]
        else
          b1 = tree.nodes[father].sons[1]
          b2 = tree.nodes[father].sons[2]
        end
      else
        b1 = father
        if tree.nodes[father].sons[1] == internode
          b2 = tree.nodes[father].sons[2]
        else
          b2 = tree.nodes[father].sons[1]
        end
      end
      #randomly select two nodes for NNI
      node1 = sample([a1,a2])
      node2 = sample([b1,b2])
    end
  
    deleteNodeNNI(tree, Int64(node1))
    addNodeNNI(tree, Int64(node1), Int64(node2))
end
        
function treeSPR(tree::Tree)
    if tree.rooted
        #find the candidate node1 to prune
        candidates = collect(1:(2*tree.ntaxa-1))
        son1 = tree.nodes[tree.rootnode].sons[1]
        son2 = tree.nodes[tree.rootnode].sons[2]
        delnode = [tree.rootnode, son1, son2]
        deleteat!(candidates, sort(delnode))  
        node1 = sample(candidates)
  
        #find the candidate node2 to insert node1
        candidates = collect(1:(2*tree.ntaxa-1))
        father = tree.nodes[node1].father
        if tree.nodes[father].sons[1] == node1
            brother = tree.nodes[father].sons[2]
        else
            brother = tree.nodes[father].sons[1]
        end
        nodes = Int64[]
        findDownNodes(tree, Int64(node1), nodes)
        nodes = cat(nodes,brother,father,tree.rootnode,dims=1)
        deleteat!(candidates, sort(nodes))
        node2 = sample(candidates)
    else
        #find the candidate node1 to prune
        candidates = collect(1:(2*tree.ntaxa-2))
        son1 = tree.nodes[tree.rootnode].sons[1]
        son2 = tree.nodes[tree.rootnode].sons[2]
        son3 = tree.nodes[tree.rootnode].sons[3]
        delnode = [tree.rootnode, son1, son2, son3]
        deleteat!(candidates, sort(delnode))  
        node1 = sample(candidates)
  
        #find the candidate node2 to insert node1
        candidates = collect(1:(2*tree.ntaxa-2))
        father = tree.nodes[node1].father
        if tree.nodes[father].sons[1] == node1
            brother = tree.nodes[father].sons[2]
        else
            brother = tree.nodes[father].sons[1]
        end
        nodes = Int64[]
        findDownNodes(tree, Int64(node1), nodes)
        nodes = cat(nodes,brother,father, tree.rootnode,dims=1)
        deleteat!(candidates, sort(nodes))
        node2 = sample(candidates)
    end
    deleteNodeNNI(tree, Int64(node1))
    addNodeNNI(tree, Int64(node1), Int64(node2))
  end
        
function treeReroot(tree::Tree, inode::Int64)
    if inode == tree.rootnode
        error("cannot root a tree using rootnode")
    end
  
    if tree.rooted
        if inode == tree.nodes[tree.rootnode].sons[1] || inode == tree.nodes[tree.rootnode].sons[2]
            println("the root does not change")
            return tree
        end
  
        #find the key nodes to update
        inodefather = tree.nodes[inode].father
        updatefather = tree.nodes[inodefather].father
        updatenode_brlens = tree.nodes[inodefather].brlens
        original_rootson1 = tree.nodes[tree.rootnode].sons[1]
        original_rootson2 = tree.nodes[tree.rootnode].sons[2]
        combinedbrlens = tree.nodes[original_rootson1].brlens+tree.nodes[original_rootson2].brlens
  
        #update rootnode
        tree.nodes[tree.rootnode].sons[1] = inode
        tree.nodes[tree.rootnode].sons[2] = inodefather 
        tree.nodes[inode].father = tree.rootnode
        tree.nodes[inode].brlens = tree.nodes[inode].brlens/2
        tree.nodes[inodefather].brlens = tree.nodes[inode].brlens
        tree.nodes[inodefather].father = tree.rootnode
  
        #update other nodes
        updatenode = inodefather
        sonnode = inode
  
        if updatefather == tree.rootnode
          if updatenode == original_rootson1
              if tree.nodes[original_rootson1].sons[1] == inode
                  tree.nodes[original_rootson1].sons[1] = original_rootson2
              else
                  tree.nodes[original_rootson1].sons[2] = original_rootson2
              end
              tree.nodes[original_rootson2].brlens = combinedbrlens
              tree.nodes[original_rootson2].father = original_rootson1
          else
              if tree.nodes[original_rootson2].sons[1] == inode
                  tree.nodes[original_rootson2].sons[1] = original_rootson1
              else
                  tree.nodes[original_rootson2].sons[2] = original_rootson1
              end
              tree.nodes[original_rootson1].brlens = combinedbrlens
              tree.nodes[original_rootson1].father = original_rootson2
          end
        else
          while updatefather != tree.rootnode   
              if tree.nodes[updatenode].sons[1] == sonnode
                  tree.nodes[updatenode].sons[1] = updatefather
              else
                  tree.nodes[updatenode].sons[2] = updatefather
              end
              x = tree.nodes[updatefather].father
              y = tree.nodes[updatefather].brlens
              tree.nodes[updatefather].father = updatenode
              tree.nodes[updatefather].brlens = updatenode_brlens
              sonnode = updatenode
              updatenode = updatefather 
              updatenode_brlens = y     
              updatefather = x       
          end
          if updatenode == original_rootson1
              if tree.nodes[updatenode].sons[1] == tree.nodes[updatenode].father
                  tree.nodes[updatenode].sons[1] = original_rootson2
              else
                  tree.nodes[updatenode].sons[2] = original_rootson2
              end
              tree.nodes[original_rootson2].father = original_rootson1
              tree.nodes[original_rootson2].brlens = combinedbrlens
          elseif updatenode == original_rootson2
              if tree.nodes[updatenode].sons[1] == tree.nodes[updatenode].father
                  tree.nodes[updatenode].sons[1] = original_rootson1
              else
                  tree.nodes[updatenode].sons[2] = original_rootson1
              end
              tree.nodes[original_rootson1].father = original_rootson2
              tree.nodes[original_rootson1].brlens = combinedbrlens
          else
              error("bug in the algorithm")
          end
        end
    else
        push!(tree.nodes,Node())
        newroot = 2*tree.ntaxa-1
        tree.nodes[newroot].nson = 2
        tree.nodes[tree.rootnode].nson = 2
        tree.rooted = true
        if inode == tree.nodes[tree.rootnode].sons[1] || inode == tree.nodes[tree.rootnode].sons[2] || inode == tree.nodes[tree.rootnode].sons[3]
            tree.nodes[newroot].sons[1] = tree.rootnode
            tree.nodes[newroot].sons[2] = inode
            tree.nodes[inode].father = newroot
            tree.nodes[inode].brlens = tree.nodes[inode].brlens/2
            tree.nodes[tree.rootnode].father = newroot
            tree.nodes[tree.rootnode].brlens = tree.nodes[inode].brlens
            deleteat!(tree.nodes[tree.rootnode].sons,findfirst(inode .== tree.nodes[tree.rootnode].sons))
            tree.rootnode = newroot
        else
            #find the key nodes to update
            inodefather = tree.nodes[inode].father
            original_rootson1 = tree.nodes[tree.rootnode].sons[1]
            original_rootson2 = tree.nodes[tree.rootnode].sons[2]
            original_rootson3 = tree.nodes[tree.rootnode].sons[3]
  
            #update rootnode
            tree.nodes[newroot].sons[1] = inode
            tree.nodes[newroot].sons[2] = inodefather
            tree.nodes[inode].father = newroot
            tree.nodes[inode].brlens = tree.nodes[inode].brlens/2
            updatefather = tree.nodes[inodefather].father
            updatenode_brlens = tree.nodes[inodefather].brlens
            tree.nodes[inodefather].father = newroot
            tree.nodes[inodefather].brlens = tree.nodes[inode].brlens
  
            #update other nodes
            updatenode = inodefather
            sonnode = inode
  
            if updatefather == tree.rootnode   
              if tree.nodes[updatenode].sons[1] == inode
                  tree.nodes[updatenode].sons[1] = tree.rootnode
              else
                  tree.nodes[updatenode].sons[2] = tree.rootnode
              end
              tree.nodes[tree.rootnode].father = updatenode 
              tree.nodes[tree.rootnode].brlens = updatenode_brlens
              deleteat!(tree.nodes[tree.rootnode].sons,findfirst(updatenode .== tree.nodes[tree.rootnode].sons))
              tree.rootnode = newroot              
            else
              while updatefather != tree.rootnode   
                  if tree.nodes[updatenode].sons[1] == sonnode
                      tree.nodes[updatenode].sons[1] = updatefather
                  else
                      tree.nodes[updatenode].sons[2] = updatefather
                  end
                  x = tree.nodes[updatefather].father
                  y = tree.nodes[updatefather].brlens
                  tree.nodes[updatefather].father = updatenode
                  tree.nodes[updatefather].brlens = updatenode_brlens
                  sonnode = updatenode
                  updatenode = updatefather
                  updatenode_brlens = y    
                  updatefather = x       
              end           
              if tree.nodes[updatenode].sons[1] == tree.nodes[updatenode].father
                  tree.nodes[updatenode].sons[1] = tree.rootnode
              else
                  tree.nodes[updatenode].sons[2] = tree.rootnode
              end
              tree.nodes[tree.rootnode].father = updatenode 
              tree.nodes[tree.rootnode].brlens = updatenode_brlens
              deleteat!(tree.nodes[tree.rootnode].sons,findfirst(updatenode .== tree.nodes[tree.rootnode].sons))
              tree.rootnode = newroot  
            end      
        end
    end
end
  
function treeRerootSubtree(tree::Tree, rootnode::Int64, inode::Int64)
    if inode == tree.nodes[rootnode].sons[1] || inode == tree.nodes[rootnode].sons[2]
        return tree
    end
  
    #nodes = Int64[]
    #findDownNodes(tree, rootnode, nodes)
    #if sum(inode .== nodes) == 0
    #  error("inode must be a descendant node of rootnode")
    #end

    #find the key nodes to update
    inodefather = tree.nodes[inode].father
    updatefather = tree.nodes[inodefather].father
    updatenode_brlens = tree.nodes[inodefather].brlens
    original_rootson1 = tree.nodes[rootnode].sons[1]
    original_rootson2 = tree.nodes[rootnode].sons[2]
    combinedbrlens = tree.nodes[original_rootson1].brlens+tree.nodes[original_rootson2].brlens
  
    #update rootnode
    tree.nodes[rootnode].sons[1] = inode
    tree.nodes[rootnode].sons[2] = inodefather
    tree.nodes[inode].father = rootnode
    tree.nodes[inode].brlens = tree.nodes[inode].brlens/2
    tree.nodes[inodefather].father = rootnode
    tree.nodes[inodefather].brlens = tree.nodes[inode].brlens
  
    #update other nodes
    updatenode = inodefather
    sonnode = inode
  
    if updatefather == rootnode
        if updatenode == original_rootson1
            if tree.nodes[original_rootson1].sons[1] == inode
                tree.nodes[original_rootson1].sons[1] = original_rootson2
            else
                tree.nodes[original_rootson1].sons[2] = original_rootson2
            end
            tree.nodes[original_rootson2].father = original_rootson1
            tree.nodes[original_rootson2].brlens = combinedbrlens
        else
            if tree.nodes[original_rootson2].sons[1] == inode
                tree.nodes[original_rootson2].sons[1] = original_rootson1
            else
                tree.nodes[original_rootson2].sons[2] = original_rootson1
            end
            tree.nodes[original_rootson1].father = original_rootson2
            tree.nodes[original_rootson1].brlens = combinedbrlens
        end
    else
        while updatefather != rootnode   
            if tree.nodes[updatenode].sons[1] == sonnode
                tree.nodes[updatenode].sons[1] = updatefather
            else
                tree.nodes[updatenode].sons[2] = updatefather
            end
            x = tree.nodes[updatefather].father
            y = tree.nodes[updatefather].brlens
            tree.nodes[updatefather].father = updatenode
            tree.nodes[updatefather].brlens = updatenode_brlens
            sonnode = updatenode
            updatenode = updatefather
            updatenode_brlens = y     
            updatefather = x       
        end
        if updatenode == original_rootson1
            if tree.nodes[original_rootson1].sons[1] == tree.nodes[original_rootson1].father
                tree.nodes[original_rootson1].sons[1] = original_rootson2
            else
                tree.nodes[original_rootson1].sons[2] = original_rootson2
            end
            tree.nodes[original_rootson2].father = original_rootson1
            tree.nodes[original_rootson2].brlens = combinedbrlens
        elseif updatenode == original_rootson2
            if tree.nodes[original_rootson2].sons[1] == tree.nodes[original_rootson2].father
                tree.nodes[original_rootson2].sons[1] = original_rootson1
            else
                tree.nodes[original_rootson2].sons[2] = original_rootson1
            end
            tree.nodes[original_rootson1].father = original_rootson2
            tree.nodes[original_rootson1].brlens = combinedbrlens
        else
            error("bug in the algorithm")
        end
    end
end

function treeReorderTaxa(tree::Tree, newtaxaname::Vector{String})
  index = collect(1:tree.ntaxa)
  for i in 1:tree.ntaxa
      index[i] = findall(x-> x == newtaxaname[i], tree.taxaname)[1]
  end

  tree.nodes[1:tree.ntaxa] = tree.nodes[index]

  index = collect(1:tree.ntaxa)
  for i in 1:tree.ntaxa
      index[i] = findall(x-> x == tree.taxaname[i], newtaxaname)[1]
  end
  for i in (tree.ntaxa+1):length(tree.nodes)
      for j in 1:tree.nodes[i].nson
          if tree.nodes[i].sons[j] < tree.ntaxa+1
              tree.nodes[i].sons[j] = index[tree.nodes[i].sons[j]]
          end
      end
  end
  tree.taxaname = newtaxaname
end

function treeRerootUnrooted(tree::Tree, inode::Int64)
    if inode == tree.rootnode
        return tree
    end

    if inode <= tree.ntaxa
      error("the new root must be an internal node")
    end

    if tree.rooted
      error("the input tree must be unrooted")
    end

    #inode is the new root
    tree.nodes[inode].nson = 3
    tree.nodes[tree.rootnode].nson = 2
    if inode == tree.nodes[tree.rootnode].sons[1] || inode == tree.nodes[tree.rootnode].sons[2] || inode == tree.nodes[tree.rootnode].sons[3]
        push!(tree.nodes[inode].sons, tree.rootnode)
        tree.nodes[tree.rootnode].father = inode
        tree.nodes[tree.rootnode].brlens = tree.nodes[inode].brlens
        deleteat!(tree.nodes[tree.rootnode].sons,findfirst(inode .== tree.nodes[tree.rootnode].sons))
        tree.rootnode = inode
        tree.nodes[inode].father = -9
        tree.nodes[inode].brlens = -9
    else
        #find the key nodes to update
        inodefather = tree.nodes[inode].father
        updatefather = tree.nodes[inodefather].father
        updatenode_brlens = tree.nodes[inodefather].brlens
        original_rootson1 = tree.nodes[tree.rootnode].sons[1]
        original_rootson2 = tree.nodes[tree.rootnode].sons[2]
        original_rootson3 = tree.nodes[tree.rootnode].sons[3]

        #update new root inode
        push!(tree.nodes[inode].sons, inodefather)
        tree.nodes[inodefather].father = inode
        tree.nodes[inodefather].brlens = tree.nodes[inode].brlens
        tree.nodes[inode].father = -9
        tree.nodes[inode].brlens = -9

        #update other nodes
        updatenode = inodefather
        sonnode = inode

        if updatefather == tree.rootnode   
          if tree.nodes[updatenode].sons[1] == inode
              tree.nodes[updatenode].sons[1] = tree.rootnode
          else
              tree.nodes[updatenode].sons[2] = tree.rootnode
          end
          tree.nodes[tree.rootnode].father = updatenode 
          tree.nodes[tree.rootnode].brlens = updatenode_brlens
          deleteat!(tree.nodes[tree.rootnode].sons,findfirst(updatenode .== tree.nodes[tree.rootnode].sons))
          tree.rootnode = inode            
        else
          while updatefather != tree.rootnode   
              if tree.nodes[updatenode].sons[1] == sonnode
                  tree.nodes[updatenode].sons[1] = updatefather
              else
                  tree.nodes[updatenode].sons[2] = updatefather
              end
              x = tree.nodes[updatefather].father
              y = tree.nodes[updatefather].brlens
              tree.nodes[updatefather].father = updatenode
              tree.nodes[updatefather].brlens = updatenode_brlens
              sonnode = updatenode
              updatenode = updatefather
              updatenode_brlens = y    
              updatefather = x       
          end           
          if tree.nodes[updatenode].sons[1] == tree.nodes[updatenode].father
              tree.nodes[updatenode].sons[1] = tree.rootnode
          else
              tree.nodes[updatenode].sons[2] = tree.rootnode
          end
          tree.nodes[tree.rootnode].father = updatenode 
          tree.nodes[tree.rootnode].brlens = updatenode_brlens
          deleteat!(tree.nodes[tree.rootnode].sons,findfirst(updatenode .== tree.nodes[tree.rootnode].sons))
          tree.rootnode = inode 
        end      
    end
end

  
function treeTBR(tree::Tree)
    if tree.rooted
        #find the candidate node1 to prune
        candidates = collect((tree.ntaxa+1):(2*tree.ntaxa-1))
        son1 = tree.nodes[tree.rootnode].sons[1]
        son2 = tree.nodes[tree.rootnode].sons[2]
        delnode = [tree.rootnode, son1, son2]
        deleteat!(candidates, sort(delnode[delnode .> tree.ntaxa] .- tree.ntaxa)) 
        if length(candidates) == 0
            error("cannot do TBR if ntaxa < 5")
        elseif length(candidates) == 1
            node1 = candidates
        else
            node1 = sample(candidates)
        end
  
        #find the candidate node2 to insert node1
        candidates = collect(1:(2*tree.ntaxa-1))
        father = tree.nodes[node1].father
        if tree.nodes[father].sons[1] == node1
            brother = tree.nodes[father].sons[2]
        else
            brother = tree.nodes[father].sons[1]
        end
        nodes = Int64[]
        findDownNodes(tree, Int64(node1), nodes)
        nodes = cat(nodes,brother,father,tree.rootnode,dims=1)
        deleteat!(candidates, sort(nodes))
        node2 = sample(candidates)
  
        #find a branch to insert the subtree
        nodes = Int64[]
        findDownNodes(tree, Int64(node1), nodes)
        if length(nodes) > 3
            inode = sample(nodes[.!(nodes .== node1)])
            treeRerootSubtree(tree, Int64(node1), Int64(inode))
        end
    else
        #find the candidate node1 to prune
        candidates = collect((tree.ntaxa+1):(2*tree.ntaxa-2))
        son1 = tree.nodes[tree.rootnode].sons[1]
        son2 = tree.nodes[tree.rootnode].sons[2]
        son3 = tree.nodes[tree.rootnode].sons[3]
        delnode = [tree.rootnode, son1, son2, son3]
        deleteat!(candidates, sort(delnode[delnode .> tree.ntaxa] .- tree.ntaxa))  
        node1 = sample(candidates)
  
        #find the candidate node2 to insert node1
        candidates = collect(1:(2*tree.ntaxa-2))
        father = tree.nodes[node1].father
        if tree.nodes[father].sons[1] == node1
            brother = tree.nodes[father].sons[2]
        else
            brother = tree.nodes[father].sons[1]
        end
        nodes = Int64[]
        findDownNodes(tree, Int64(node1), nodes)
        nodes = cat(nodes,brother,father, tree.rootnode,dims=1)
        deleteat!(candidates, sort(nodes))
        node2 = sample(candidates)
  
        #find a branch to insert the subtree
        nodes = Int64[]
        findDownNodes(tree, Int64(node1), nodes)
        if length(nodes) > 3
            inode = sample(nodes[.!(nodes .== node1)])
            treeRerootSubtree(tree, Int64(node1), Int64(inode))
        end
    end
    deleteNodeNNI(tree, Int64(node1))
    addNodeNNI(tree, Int64(node1), Int64(node2))
end
 
function treeGenerator(taxanames::Vector{String}, type::String)
  ntaxa = length(taxanames)
  tree = Tree(ntaxa, false)
  brlens = 0.01

  if ntaxa < 4
      error("the number of taxa should be greater than 3")
  end

  #define root
  tree.rootnode = ntaxa+1
  tree.nodes[tree.rootnode].father = -9
  tree.nodes[tree.rootnode].nson = 3
  tree.nodes[tree.rootnode].sons[1] = ntaxa
  tree.nodes[tree.rootnode].sons[2] = ntaxa-1
  tree.nodes[tree.rootnode].sons[3] = tree.rootnode+1

  tree.nodes[ntaxa].father = tree.rootnode
  tree.nodes[ntaxa].brlens = brlens
  tree.nodes[ntaxa-1].father = tree.rootnode
  tree.nodes[ntaxa-1].brlens = brlens

  #define other internodes
  for i in (ntaxa+2) : (2*ntaxa-3)
      tree.nodes[i].father = i-1
      tree.nodes[i].nson = 2
      tree.nodes[i].sons[1] = i+1
      tree.nodes[i].sons[2] = 2*ntaxa-i
      tree.nodes[i].brlens = brlens

      tree.nodes[2*ntaxa-i].father = i
      tree.nodes[2*ntaxa-i].brlens = brlens
  end

  #define the last internode
  tree.nodes[2*ntaxa-2].father = 2*ntaxa-3
  tree.nodes[2*ntaxa-2].nson = 2
  tree.nodes[2*ntaxa-2].sons[1] = 1
  tree.nodes[2*ntaxa-2].sons[2] = 2
  tree.nodes[2*ntaxa-2].brlens = brlens

  tree.nodes[1].father = 2*ntaxa-2
  tree.nodes[1].brlens = brlens
  tree.nodes[2].father = 2*ntaxa-2
  tree.nodes[2].brlens = brlens

  tree.taxaname = taxanames
  for i in 1:ntaxa
      tree.nodes[i].taxaname = taxanames[i]
  end
  return tree
end


      