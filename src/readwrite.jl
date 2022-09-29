function readDNA(file, format="nexus")
    X = readlines(file)
    X = X[X .!= ""]
    X = replace.(X, "\t"=>" ")
    X = replace.(X, r"\[.*\]"=>"")

    if format == "nexus"
        i = occursin.(r"dimension"i, X)
        str = X[i][1]
        str = replace(str, " "=>"")
        str = collect(eachmatch(r"=[0-9]+", str))
        ntaxa = parse(Int, replace(str[1].match, "="=>""))
        seqlength = parse(Int, replace(str[2].match, "="=>""))

        i = findall(==(1),occursin.(r"matrix"i, X))[1]
        str = X[i:end]
        j = findall(==(1),occursin.(r";", str))[1]
        str = str[2:(j-1)]
    
        taxaname = fill("", ntaxa)
        sequence = fill("", ntaxa)

        m = match.(r"[-?\w]+", str[1:ntaxa])
        for i in 1:ntaxa
            taxaname[i] = m[i].match
        end

        for i in eachindex(str) 
            if i <= ntaxa
                x = replace(str[i], taxaname[i]=>"")
                x = replace(x,r" "=>"")
                sequence[i] = x
            else
                index = i%ntaxa
                if index == 0
                    index = ntaxa
                end
                x = replace(str[i], taxaname[index]=>"")
                x = replace(x, r" "=>"")
                sequence[index] = sequence[index] * x
            end
        end
        
        if length(sequence[1]) != seqlength
            error("the actual sequence length != nchar")
        end

        index = findall(==(1),occursin.(r"charset"i, X))
        if length(index) > 0
            str = X[index]
            str = replace.(str, r"charset"i=>"")
            str = replace.(str, r"\s"=>"")
            genename = collect.(eachmatch.(r"\w+=", str))
            location_a = collect.(eachmatch.(r"\w+-", str))
            location_b = collect.(eachmatch.(r"\w+;", str))
            gene = fill("", (length(index),3))
            for i in 1:length(index)
                gene[i,1] = genename[i][1].match
                gene[i,2] = location_a[i][1].match
                gene[i,3] = location_b[i][1].match
            end
            gene = replace.(gene,r"[-;=]"=>"")
        else 
            gene = missing
        end   
    elseif format == "phylip" 
        str = collect(eachmatch(r"[0-9]+",X[1]))
        ntaxa = parse(Int, str[1].match)
        seqlength = parse(Int, str[2].match)
    
        taxaname = fill("", ntaxa)
        sequence = fill("", ntaxa)      
        str = X[2:end]
        m = match.(r"[-?\w]+", str[1:ntaxa])
        for i in 1:ntaxa
            taxaname[i] = m[i].match
        end

        for i in 1:length(str) 
            if i <= ntaxa
                x = replace(str[i], taxaname[i]=>"")
                x = replace(x,r" "=>"")
                sequence[i] = x
            else
                index = i%ntaxa
                if index == 0
                    index = ntaxa
                end
                x = replace(str[i], taxaname[index]=>"")
                x = replace(x, r" "=>"")
                sequence[index] = sequence[index] * x
            end
        end
        
        if length(sequence[1]) != seqlength
            error("the actual sequence length != nchar")
        end
        gene = missing
    elseif format == "fasta"
        index = findall(==(1),occursin.(r">"i, X))
        ntaxa = length(index)
        taxaname = replace.(X[index],">"=>"")
        sequence = fill("", ntaxa)
        for i in 1:ntaxa
            a = index[i]+1
            if i < ntaxa
                b = index[i+1]-1
            else
                b = length(X)
            end
            sequence[i] = string(X[a:b])
        end
        sequence = replace.(sequence, r"[\"\[\s\],]"=>"")   
        gene = missing
    else
        error("format must be nexus or phylip or fasta")
    end

    if ismissing(gene) 
      result = SEQdata(taxaname, uppercase.(sequence))
    else 
      result = SEQdata(taxaname, uppercase.(sequence), gene)
    end

    return result
end

function writeDNA(data::SEQdata; file::String, format="nexus")
    f = open(file, "w")
    format = lowercase(format)

    if format == "fasta"
        for i in 1:length(data.taxaname)
            write(f, ">", data.taxaname[i], "\n", data.sequence[i],"\n")
        end
    elseif format == "phylip"
        ntaxa = length(data.taxaname)
        seqlength = length(data.sequence[1])
        write(f, "$ntaxa  $seqlength\n")
        for i in 1:length(data.taxaname)
            write(f, data.taxaname[i], "     ", data.sequence[i],"\n")
        end
    elseif format == "nexus"
        ntaxa = length(data.taxaname)
        seqlength = length(data.sequence[1])
        x = "#NEXUS\nBEGIN DATA;\nDIMENSIONS  NTAX=$ntaxa NCHAR=$seqlength;\nFORMAT DATATYPE=DNA  INTERLEAVE=YES MISSING=? GAP=-;\nMATRIX\n" 
        write(f, x)
        for i in 1:length(data.taxaname)
            write(f, data.taxaname[i], "     ", data.sequence[i],"\n")
        end
        write(f,";\nEND;")
    else
        error("The output format must be fasta or phylip or nexus")
    end

    close(f)
end

function readTree(str)
    if typeof(str) == Vector{String}
      tree = Tree[]
      for i in 1:length(str)
        push!(tree, readAtree(str[i]))
      end
    elseif typeof(str) == String
      if str[1] == '('
        tree = readAtree(str)
      else
        treestr = readlines(str)
        if length(treestr) == 1
          tree = readAtree(treestr[1])
        else
          tree = Tree[]
          for i in 1:length(treestr)
            push!(tree, readAtree(treestr[i]))
          end
        end
      end
    end
    return tree
end
  
function readAtree(treestr::String)
  treestr = replace(treestr, r"\s"=>"");
  ntaxa = count(j->j==',', treestr)+1
  leftp = count(j->j=='(', treestr)
  rightp = count(j->j==')', treestr)

  if leftp != rightp
    error("the number of ( != the number of )")
  end

  if treestr[end] != ';'
    error("the tree must end with ;")
  end

  if rightp == ntaxa-1
    rooted = true
  elseif rightp == ntaxa - 2
    rooted = false
  else 
    rooted = false
    #error("the number of ) is wrong")
  end

  tree = Tree(ntaxa, rooted)
  nnode = ntaxa
  cfather = 0
  taxa = 1
  inodeb = 0
  i = 1

  while treestr[i] != ';'
    if treestr[i] == '('
      cnode = nnode+1
      nnode += 1            
      if cfather >= 1 
        tree.nodes[cfather].nson += 1
        tree.nodes[cfather].sons[tree.nodes[cfather].nson] = cnode;
        tree.nodes[cnode].father=cfather;
      else
        tree.rootnode=cnode;
        cfather=cnode;
      end
      cfather = cnode
      i += 1
    elseif treestr[i] == ')'  
      inodeb = cfather
      cfather = tree.nodes[cfather].father
      j = i
      while treestr[i+1] != ',' && treestr[i+1] != ':' && treestr[i+1] != ';' && treestr[i+1] != ')'
        i += 1
      end
      tree.nodes[inodeb].label = treestr[(j+1):i]
      i += 1
    elseif treestr[i] == ':'
      i = i+1
      j = i
      while isdigit(treestr[i]) || treestr[i] == '.' 
        i += 1
      end
      tree.nodes[inodeb].brlens = parse(Float64, treestr[j:(i-1)])
    elseif treestr[i] == '#'
      i = i+1
      j = i
      while isdigit(treestr[i]) || treestr[i] == '.' 
        i += 1
      end
      tree.nodes[inodeb].theta = parse(Float64, treestr[j:(i-1)])
    elseif treestr[i] == '['
      i = i+1
      j = i
      while treestr[i] != ']' 
        i += 1
      end
      tree.nodes[inodeb].label = treestr[j:(i-1)]
    elseif isletter(treestr[i])
      j = i
      while treestr[i] != ',' && treestr[i] != ':' && treestr[i] != ')'
        i += 1
      end
      taxaname = treestr[j:(i-1)]
      tree.nodes[taxa].taxaname = taxaname
      tree.taxaname[taxa] = taxaname
      tree.nodes[taxa].father = cfather
      tree.nodes[cfather].nson += 1
      tree.nodes[cfather].sons[tree.nodes[cfather].nson] = taxa
      inodeb = taxa;
      taxa += 1;
    else
      i += 1
    end
  end
  return tree
end

  
function writeSubtree(inode::Int64, tree::Tree, printBrlens::Bool, printTheta::Bool, printLabel::Bool)
    if tree.nodes[inode].sons[1] > 0
        if inode == tree.rootnode && !tree.rooted
          son1 = tree.nodes[inode].sons[1]
          son2 = tree.nodes[inode].sons[2]
          son3 = tree.nodes[inode].sons[3]
          
          str1 = writeSubtree(son1, tree, printBrlens, printTheta, printLabel)
          x = "($str1"
          brlens = tree.nodes[son1].brlens
          theta = tree.nodes[son1].theta
          label = tree.nodes[son1].label
          if printLabel && label != ""
            x = "$x$label"
          end
          if printBrlens && brlens>0
            x = "$x:$brlens"
          end
          if printTheta && theta>0
            x = "$x#$theta"
          end
        
          str2 = writeSubtree(son2, tree, printBrlens, printTheta, printLabel)
          x = "$x,$str2"
          brlens = tree.nodes[son2].brlens
          theta = tree.nodes[son2].theta
          label = tree.nodes[son2].label
          if printLabel && label != ""
            x = "$x$label"
          end
          if printBrlens && brlens>0
            x = "$x:$brlens"
          end
          if printTheta && theta>0
            x = "$x#$theta"
          end
          
          str3 = writeSubtree(son3, tree, printBrlens, printTheta, printLabel) 
          x = "$x,$str3"
          brlens = tree.nodes[son3].brlens
          theta = tree.nodes[son3].theta
          label = tree.nodes[son3].label
          if printLabel && label != ""
            x = "$x$label"
          end
          if printBrlens && brlens>0
            x = "$x:$brlens"
          end
          if printTheta && theta>0
            x = "$x#$theta"
          end
          x = "$x)"
        else
          son1 = tree.nodes[inode].sons[1]
          son2 = tree.nodes[inode].sons[2]
  
          str1 = writeSubtree(son1, tree, printBrlens, printTheta, printLabel)
          x = "($str1"
          brlens = tree.nodes[son1].brlens
          theta = tree.nodes[son1].theta
          label = tree.nodes[son1].label
          if printLabel && label != ""
            x = "$x$label"
          end
          if printBrlens && brlens>0
            x = "$x:$brlens"
          end
          if printTheta && theta>0
            x = "$x#$theta"
          end
          
          str2 = writeSubtree(son2, tree, printBrlens, printTheta, printLabel)
          x = "$x,$str2"
          brlens = tree.nodes[son2].brlens
          theta = tree.nodes[son2].theta
          label = tree.nodes[son2].label
          if printLabel && label != ""
            x = "$x$label"
          end
          if printBrlens && brlens>0
            x = "$x:$brlens"
          end
          if printTheta && theta>0
            x = "$x#$theta"
          end
          x = "$x)"
        end
    else 
      x = tree.nodes[inode].taxaname
    end
  
    if inode == tree.rootnode
        x = "$x;"
    end
    
    return x
end
  
function writeTree(tree=Tree[], file="text", printBrlens=false, printTheta=false, printLabel=false)
    if typeof(tree) == Tree
        ntree = 1
    elseif typeof(tree) == Vector{Tree}
        ntree = length(tree)
    else
        error("The type of input trees is not Tree")
    end

    if ntree == 1
        rootnode = tree.rootnode
        treestr = writeSubtree(rootnode, tree, printBrlens, printTheta, printLabel)
    end

    if ntree > 1
        treestr = fill("", ntree)
        for i in 1:ntree
            rootnode = tree[i].rootnode
            treestr[i] = writeSubtree(rootnode, tree[i], printBrlens, printTheta, printLabel)
        end
    end

    if file == "text"
        return treestr
    else
        if ntree == 1 
            write(file, treestr)
        else
            write(file, join(treestr,"\n"))
        end
    end
end
  
  
