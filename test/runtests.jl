using PhyLibrary
using Distributions

####################################################
# simulation
####################################################
basefreq = [1/4,1/4,1/4,1/4]
rates = [0.1,0.2,0.3,0.4,0.5,1.0]
brlens = 0.01

states = ['A','C','G','T']
data = simNucleotideSingleBranch(states, 1000, basefreq, rates, brlens)

treestr = "((((S1:0.01,S4:0.1):0.1,S3:0.21):0.09,S2:0.1):0.1,S5:0.09,S6:0.2);"
tree = readTree(treestr)
data = simNucleotide(states, tree, basefreq, rates, 1000)

#####################################################
# estimation
#####################################################
estimates = submodelGTR_rateMatrix_estimation(data.sequence, states)

##################################################
# read DNA sequences
##################################################
nexusfile = joinpath(@__DIR__, "..", "data", "test.nex")
fastafile = joinpath(@__DIR__, "..", "data", "test.fasta")
phyfile = joinpath(@__DIR__, "..", "data", "test.phy")

w=readDNA(fastafile, "fasta")
w=readDNA(nexusfile);
w=readDNA(phyfile,"phylip");
writeDNA(w, file="testwrite.nex", format="nexus")

states = ['A','C','G','T']

sequence_base_freq(w.sequence, states)
sequence_doublet_freq(w.sequence, states)

###########################################################
# read and write a tree
#############################################################
treestr1 = "(((((a:1.0,b:1.0):2.0,c:3.0):4.0,d:5.0):0.1,e:0.1):0.1,f:0.7);"
treestr2 = "(((S1:0.1,S2:0.1#0.345)6.45:0.1,S3:0.1)[3.4]:0.09,S4:0.1,S5:0.09);"
tree1 = readTree(treestr1)
tree2 = readTree(treestr2)
writeTree(tree1, "testwritetree.tre", true)

#########################################################
# delete or add a node in tree
#########################################################
treestr = "(((S1:0.1,S2:0.1#0.345)6.45:0.1,S3:0.1)[3.4]:0.09,S4:0.1,S5:0.09);"
tree = readTree(treestr)
leaf = Int64[]
#findDownLeaves(tree, Int64(9),leaf)
#findDownNodes(tree, Int64(9),leaf)

deleteNode(tree, Int64(7))
writeTree(tree)

addNode(tree,Int64(8),Int64(5))
writeTree(tree)

swapNodes(tree, Int64(4),Int64(1))
writeTree(tree)

#############################################################
# tree rearrangement
#############################################################
treestr = "(((((a:1.0,b:1.0):2.0,c:3.0):4.0,d:5.0):0.1,e:0.1):0.1,f:0.7);"
tree = readTree(treestr)
treeRerootSubtree(tree,9,1)
writeTree(tree,"text",true)

treestr = "((((S1:0.1,S2:0.1#0.345)6.45:0.01,S3:0.15):0.09,S4:0.1):0.021,S5:0.109,S6:0.011);"
tree = readTree(treestr)
treeReroot(tree,Int64(1))
writeTree(tree,"text",true)

treestr = "((((S1:0.1,S2:0.1#0.345)6.45:0.01,S3:0.15):0.09,S4:0.1):0.021,S5:0.109,S6:0.011);"
tree = readTree(treestr)
treeRerootUnrooted(tree, 9)
writeTree(tree,"text",true)

treestr = "((((S1:0.1,S2:0.1#0.345)6.45:0.01,S3:0.15):0.09,S4:0.1):0.021,S5:0.109,S6:0.011);"
tree = readTree(treestr)
treeNNI(tree)
writeTree(tree,"text",true)
treeTBR(tree)
writeTree(tree,"text",true)

for i in 1:100
    treeTBR(tree)
    println(writeTree(tree))
end

##################################################################
# tree distance
##################################################################
treestr1 = "((((S1:0.1,S2:0.1):0.1,S3:0.1):0.09,S4:0.1):0.1,S5:0.09,S6);"
treestr2 = "((((S1:0.1,S4:0.1):0.1,S3:0.1):0.09,S2:0.1):0.1,S5:0.09,S6);"
tree1 = readTree(treestr1)
tree2 = readTree(treestr2)

findAncestorDistance(tree2,1)
anc = findAllAncestorDistance(tree2)
findSingleQuartet(tree2, "S5","S2","S3","S4", anc)
findAllQuartet(tree2)
treeQuartetDistance(tree1,tree2)

treestr1 = "(((((a:1.0,b:1.0):2.0,c:3.0):4.0,d:5.0):0.1,e:0.1):0.1,f:0.7);"
treestr2 = "(((((a:1.0,d:1.0):2.0,c:3.0):4.0,b:5.0):0.1,e:0.1):0.1,f:0.7);"

tree1 = readTree(treestr1)
tree2 = readTree(treestr2)

treeSinglePartition(tree1, 9)
treeRFdistance(tree1,tree2)

findPathDistance(tree1,1,3)
findAllPathDistance(tree1)

anc = findAllAncestorDistance(tree2)
findSingleTriple(tree1, "a","b","d", anc)
findAllTriple(tree1)
treeTripleDistance(tree1,tree2)

######################################################################################
# likelihood calculation for 4 states
######################################################################################
basefreq = [1/4,1/4,1/4,1/4]
rates = [1,1,1,1,1.0,1.0]
qmatrix = submodelGTR_rateMatrix(basefreq, rates)

treestr = "((((S1:0.01,S4:0.1):0.1,S3:0.21):0.09,S2:0.1):0.1,S5:0.09,S6:0.2);"
tree = readTree(treestr)
states = ['A','C','G','T']
modelpara = ModelParameters(states, zeros(2*tree.ntaxa-2,length(states)),qmatrix,basefreq)

it = Iterators.product(ntuple(_ -> states, 6)...)
data = collect(it);
npattern = length(data)
loglike = zeros(npattern)
for i in 1:npattern
    loglike[i] = LoglikelihoodSingleDNA(join(data[i]),tree,modelpara)
end

sum(exp.(loglike))

#################################################################################
# likelihood calculation for 3 states
################################################################################
basefreq = [0.25,0.25,0.5]
rates = [1.5,0.7,1.0]
qmatrix = submodelGTR_rateMatrix(basefreq, rates)
treestr = "((((S1:0.01,S4:0.1):0.1,S3:0.21):0.09,S2:0.1):0.1,S5:0.09,S6:0.2);"
tree = readTree(treestr)
states = ['0','1','2']
modelpara = ModelParameters(states, zeros(2*tree.ntaxa-2,length(states)),qmatrix,basefreq)

it = Iterators.product(ntuple(_ -> modelpara.states, 6)...)
data1 = collect(it);
npattern = length(data1)
loglike = zeros(npattern)
for i in 1:npattern
    loglike[i] = LoglikelihoodSingleDNA(join(data1[i]),tree,modelpara)
end
sum(exp.(loglike))

####################################################################################
# an example
######################################################################################
phyfile = joinpath(@__DIR__, "..", "data", "sim.phy")
data = readDNA(phyfile,"phylip");
states = ['A','C','G','T']

doublet = sequence_doublet_freq_pair(data, states)
doubletfreq = doublet["doubletfreq"]
pattern_freq = alignmentPattern(data.sequence)
pattern = names(pattern_freq)[1]
patternfreq = Vector(pattern_freq)

basefreq = [1/4,1/4,1/4,1/4]
rates = [1,1,1,1,1.0,1.0]
qmatrix = submodelGTR_rateMatrix(basefreq, rates)

treestr = "((((S1:0.01,S4:0.1):0.1,S3:0.21):0.09,S2:0.1):0.1,S5:0.09,S6:0.2);"
tree = readTree(treestr)
treeReorderTaxa(tree, data.taxaname)
modelpara = ModelParameters(states, zeros(2*tree.ntaxa-2,4), qmatrix, basefreq)

LoglikelihoodDNA(patternfreq, pattern, tree, modelpara)
LoglikelihoodDNApairwise(doubletfreq, tree, modelpara)

submodelGTR_rateMatrix_estimation(data.sequence, states)

log1=zeros(20)
log2=zeros(20)
for i in 1:20
    treeNNI(tree)
    log1[i] = LoglikelihoodDNA(patternfreq, pattern, tree, modelpara)
    log2[i] = LoglikelihoodDNApairwise(doubletfreq, tree, modelpara)
end


br = collect(1:20)
log1=zeros(20)
log2=zeros(20)
for i in 1:20
    tree.nodes[8].brlens = br[i]
    log1[i] = LoglikelihoodDNA(patternfreq, pattern, tree, modelpara)
    log2[i] = LoglikelihoodDNApairwise(doubletfreq, tree, modelpara)
end



####################################################################################
## simulate and calculate ML for DiMethyl model
####################################################################################

states = ['1','2','3']
treestr = "((((S1:0.01,S4:0.1):0.1,S3:0.21):0.09,S2:0.1):0.1,S5:0.09,S6:0.2);"
tree = readTree(treestr)
data = simDiMethyl(states, tree, 1.0, 2.0, 1000)
treeMaximumLikelihood2(data, states, treeGenerator(data.taxaname,""),"DiMethyl")


#########################################################################
# estimation brlens
#########################################################################
phyfile = joinpath(@__DIR__, "..", "data", "sim.phy")
data = readDNA(phyfile,"phylip");
states = ['A','C','G','T']
doublet = sequence_doublet_freq_pair(data, states)
doubletfreq = doublet["doubletfreq"]
pattern_freq = alignmentPattern(data.sequence)
pattern = names(pattern_freq)[1]
patternfreq = Vector(pattern_freq)

basefreq = [1/4,1/4,1/4,1/4]
rates = [1,1,1,1,1.0,1.0]
qmatrix = submodelGTR_rateMatrix(basefreq, rates)

treestr = "((((S1:0.01,S4:0.1):0.1,S3:0.21):0.09,S2:0.1):0.1,S5:0.09,S6:0.2);"
tree = readTree(treestr)
modelpara = ModelParameters(states, zeros(2*tree.ntaxa-2,4), qmatrix, basefreq)

#calculate likelihood
LoglikelihoodDNA(patternfreq, pattern, tree, modelpara)

#change the tree root, make sure that the likleihood remains the same
treeRerootUnrooted(tree,9)
LoglikelihoodDNA(patternfreq, pattern, tree, modelpara)

#optimize branch lengths using Newton-Raphson 
for i in 1:tree.ntaxa
    treeNewtonRaphsonSingleBrlensFull(tree, i, patternfreq, pattern, modelpara)
end

treeNewtonRaphsonBrlens(tree, patternfreq, pattern, modelpara)

#######################################################################
# maximum likelihood phylogenetic tree
#######################################################################
phyfile = joinpath(@__DIR__, "..", "data", "sim.phy")
data = readDNA(phyfile,"phylip");
states = ['A','C','G','T']
doublet = sequence_doublet_freq_pair(data, states)
doubletfreq = doublet["doubletfreq"]

treeMaximumLikelihood(data, states, treeGenerator(data.taxaname,""))

#################################################
# triplet frequency
#################################################
treefile = joinpath(@__DIR__, "..", "data", "intron_root_mle_outlier0.7.tre")
taxon = ["cuclscnrs","chlmydtsm","trcrythrl","ptrclsgtt","ptgnsfsct","antrstmsc","atlntsrgr","brhnsbstr","pdcpscrst"]
tripleFrequancy(treefile, taxon)

tree = readTree(treefile)
treeSpeciesFrequency(tree)

treestr = ["((a,b),c,d);","((a,c),b,d);","((a,b),c,d);","((a,c),b,d);","((a,c),b,d);","((a,d),b,c);"]
tree = readTree(treestr)
treeDistribution(tree)

####################################################
# tree HypothesisTests
######################################################
using CSV
using DataFrames

file = joinpath(@__DIR__, "..", "data", "cds.csv")
data = CSV.read(file, DataFrame)
ntriple = size(data)[1]
result = zeros(ntriple,2)
for i in 1:ntriple
    x = Vector(data[i,5:7])
    y = tripleTest(x)
    result[i,1] = y[1]
    result[i,2] = y[2]
end
data[!,:MajorTest] = result[:,1]
data[!,:MinorTest] = result[:,2]
    

#########################################################
# consensus
#########################################################
treestr = "((((a,b),c),d),e,f);"
tree = readTree(treestr)
multree = Tree[]
push!(multree,deepcopy(tree))
for i in 1:99
    treeNNI(tree)
    push!(multree,deepcopy(tree))
end

treeConsensus(multree)

######################################################
# tree plot
#######################################################
treestr = ["((a,b),c,d);","((a,c),b,d);","((a,b),c,d);","((a,c),b,d);","((a,c),b,d);","((a,d),b,c);"]
plotFigtree(treestr)

###################################################
#
###################################################
file = joinpath(@__DIR__, "..", "data", "RAxML_MajorityRuleExtendedConsensusTree.S.con.tre")
treestr = readline(file);
newtree = convertRAxMLconsensusTree(treestr)
