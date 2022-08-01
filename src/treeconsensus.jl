
function treeConsensus(multree::Vector{Tree})
    ntree = length(multree)
    ntaxa = multree[1].ntaxa
    taxaname = sort(multree[1].taxaname) 
    taxabrlens = Float64[]
    parlens = Float64[]

    #do we need to calculate branch lengths?
    if multree[1].nodes[1].brlens > 0
        calbrlens = true
    else
        calbrlens = false
    end

    if calbrlens
        brlens = zeros(2*ntaxa-3, ntree)
    end

    #find all partitions across trees
    npartition = 2*ntaxa-3
    str = String[]

    for j in 1:ntree
        partition = treeAllPartition(multree[j])
        partition = partition[:, sortperm(multree[j].taxaname)]
        Threads.@spawn for i in (ntaxa+1):npartition
            push!(str, join(Int.(partition[i,:])))
            push!(str, join(Int.(.!partition[i,:])))
        end
        #branch length
        if calbrlens > 0
            indicator = 1
            for i in 1:(2*ntaxa-2)
                if i <= ntaxa
                    brlens[(1:ntaxa)[partition[i,:]][1],j] = multree[j].nodes[i].brlens
                    indicator += 1
                elseif i == multree[j].rootnode
                    continue
                else
                    brlens[indicator,j] = multree[j].nodes[i].brlens
                    indicator += 1
                end
            end
        end
    end

    #sort partitions by frequency
    freq = freqtable(str)
    partition = sort(freq,rev=true)
    parstr = names(partition)[1]

    #calculate the branch lengths of sorted partitions. the number of internal branches = ntaxa-3
    if calbrlens
        taxabrlens = zeros(ntaxa)
        for i in 1:ntaxa
            taxabrlens[i] = sum(brlens[i,:])/ntree
        end

        brlens = brlens[(ntaxa+1):npartition,:]
        parlens = zeros(length(parstr))
        for i in 1:length(parstr)
            x = 1:length(str)
            y = fld.(x[parstr[i] .== str] .- 1,2) .+ 1
            parlens[i] = sum(brlens[y])/length(y)
        end
    end
    
    ##############################################################
    #find the most frequent partitions of the consensus tree
    ##############################################################
    parmatrix = fill(2,(length(parstr),ntaxa))
    for i in 1:length(parstr)
        parmatrix[i,:] = parse.(Int64,split(parstr[i],""))
    end
    parmatrix = hcat(parmatrix,reshape(partition,(length(partition),1)))
    npartition = size(parmatrix)[1]
    indicator = 1
    while npartition != ntaxa - 3 && indicator < ntaxa
        index = fill(true, npartition)
        for i in (indicator+1):npartition
            x = parmatrix[i,1:ntaxa] - parmatrix[indicator,1:ntaxa]
            y = parmatrix[i,1:ntaxa] + parmatrix[indicator,1:ntaxa]
            #find the enclosed or exclusive partitions 
            if (sum(x .== 1) > 0 && sum(x .== -1) > 0 && sum(y .== 2) > 0) || (sum(y .== 1) == ntaxa)
                index[i] = false
            end
        end
        indicator += 1
        parmatrix = parmatrix[index,:]
        if calbrlens
            parlens = parlens[index]
        end
        npartition = sum(index)
    end

    treestr = partition2tree(parmatrix, taxaname, taxabrlens, parlens)

    return treestr
end


function partition2tree(parmatrix::Matrix{Int64}, taxaname::Vector{String}, taxabrlens::Vector{Float64}, parbrlens::Vector{Float64})
    n = size(parmatrix)
    ntaxa = length(taxaname)
    npartition = n[1]
    if n[2] != ntaxa+1
        error("The partition matrix does not match taxa names")
    end

    if length(parbrlens) > 0
        calbrlens = true
    else
        calbrlens = false
    end

    if calbrlens
        if length(parbrlens) != npartition
            error("The partition matrix does not match branch lengths")
        end
        taxaname = string.(taxaname,":", taxabrlens[1:ntaxa])
    end

    while sum(parmatrix[:,1:ntaxa]) > npartition
        for i in 1:npartition
            if sum(parmatrix[i,1:ntaxa]) == 2
                x = (1:ntaxa)[parmatrix[i,1:ntaxa] .== 1]
                
                if calbrlens
                    taxaname[x[1]] = join(["(",taxaname[x[1]],",",taxaname[x[2]],")",parmatrix[i,ntaxa+1],":",parbrlens[i]])
                else
                    taxaname[x[1]] = join(["(",taxaname[x[1]],",",taxaname[x[2]],")",parmatrix[i,ntaxa+1]])
                end

                taxaname[x[2]] = ""
                parmatrix[:,x[2]] .= 0  
            end
        end
    end

    treestr = join(["(",join(taxaname[taxaname .!= ""],","),")100;"])
    return treestr
end



###################################################################################
# the consensus tree generated by RAxML does not have a proper format
# this function can correct the format and open the tree file in figtree
###################################################################################
function convertRAxMLconsensusTree(treestr::String)
    treestr1 = collect(treestr)
    for m in eachmatch(r":[[:digit:]]*.[[:digit:]]*\[[[:digit:]]*\]", treestr)
        i = m.offset
        j = i + length(m.match) - 1
        str = treestr[i:j]
        w = match(r"\[",str).offset
        str = str[w:end] * str[1:(w-1)]
        treestr1[i:j] = collect(str)
    end
    return join(treestr1)
end
