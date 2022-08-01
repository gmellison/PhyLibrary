function sequence_doublet_freq(data::Vector{String}, states::Vector{Char})
    nstates = length(states)
    it = Iterators.product(ntuple(_ -> states, 2)...)
    doubletall = vec(permutedims(join.(collect(it))))
    ndoublet = length(doubletall)

    str = split.(data,"")
    n = length(data)
    doublet_freq = zeros(ndoublet)
    for i in 1:(n-1)
        for j in (i+1):n
            doublet = str[i] .* str[j]
            for k in 1:ndoublet
                doublet_freq[k] = doublet_freq[k] + sum(doublet .== doubletall[k])
            end
        end
    end
    return doublet_freq/sum(doublet_freq)
end

function sequence_base_freq(data::Vector{String}, states::Vector{Char})
    nstates = length(states)
    seq = string(data)
    basefreq_hat = zeros(nstates)
    for i in 1:nstates
        basefreq_hat[i] = count(j->j==states[i], seq)
    end
    basefreq_hat = basefreq_hat/sum(basefreq_hat)
end

function alignmentPattern(sequence::Vector{String})
    nchar = length(sequence[1])
    ntaxa = length(sequence)
    dna = fill("", ntaxa, nchar)
    for i in 1:ntaxa
        dna[i,:] = split(sequence[i],"")
    end
    
    site = fill("",nchar)
    for i in 1:nchar
        site[i] = join(dna[:,i])
    end
    return freqtable(site)
end

function sequence_doublet_freq_pair(data::SEQdata, states::Vector{Char})
    ntaxa = length(data.taxaname)
    npair = floor(Int, ntaxa*(ntaxa-1)/2)
    pairtaxa = fill("", (npair, 2))
    str = split.(data.sequence, "")

    nstates = length(states)
    it = Iterators.product(ntuple(_ -> states, 2)...)
    doubletall = vec(permutedims(join.(collect(it))))
    ndoublet = nstates^2
  
    doublet_freq = fill(0, npair, ndoublet)
    index = 1
    for i in 1:(ntaxa-1)
      for j in (i+1):ntaxa
        doublet = str[i] .* str[j]
        for k in 1:ndoublet
          doublet_freq[index, k] = sum(doublet .== doubletall[k])
        end
        pairtaxa[index, 1] = data.taxaname[i]
        pairtaxa[index, 2] = data.taxaname[j]
        index = index + 1
      end
    end
    result = Dict([("doubletfreq", doublet_freq),("pair", pairtaxa)])
    return result
end
  
