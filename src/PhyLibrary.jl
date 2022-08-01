module PhyLibrary

    using StatsBase
    using LinearAlgebra
    using Distributions
    using FreqTables
    using Tables
    using HypothesisTests
    
    const ERROR = 1e-5

    export
        #data and tree definition    
        Tree,
        Node,
        SEQdata,
        ModelParameters,

        #estimation
        submodelGTR_rateMatrix_estimation,
        treeNewtonRaphsonSingleBrlensFull,
        treeNewtonRaphsonSingleBrlensPairwise,
        treeNewtonRaphsonBrlens,
        treeMaximumLikelihood,

        #likelihood
        preorderTraversal,
        LoglikelihoodDNA,
        LoglikelihoodSingleDNA,
        LoglikelihoodDNApairwise,

        #matrix 
        lower_tri,
        upper_tri,
        diag,

        #read and write
        readDNA,
        readAtree,
        readTree,
        writeDNA,
        writeSubtree,
        writeTree,

        #sequence
        sequence_doublet_freq,
        sequence_base_freq,
        alignmentPattern,
        sequence_doublet_freq_pair,

        #simulation
        simNucleotideSingleBranch,
        simNucleotide,

        #substitution models
        submodelGTR_rateMatrix,

        #tree consensus
        treeConsensus,
        convertRAxMLconsensusTree,

        #tree distance
        findAncestorDistance,
        findAllAncestorDistance,
        findPathDistance,
        findAllPathDistance,
        findPathDistance,
        treeSinglePartition,
        treeAllPartition,
        treeRFdistance,
        findSingleTriple,
        findAllTriple,
        tripleFrequancy,
        treeTripleDistance,
        findSingleQuartet,
        findAllQuartet,
        treeQuartetDistance,
        treeSpeciesFrequency,
        treeDistribution,

        #tree plot
        plotFigtree,

        #tree rearrangement
        deleteNode,
        addNode,
        swapNodes,
        treeNNI,
        treeSPR,
        treeTBR,
        treeReorderTaxa,
        treeReroot,
        treeRerootSubtree,
        treeRerootUnrooted,
        treeGenerator,

        #tree HypothesisTests
        tripleTest

    include("struct.jl")
    include("estimation.jl")
    include("likelihood.jl")
    include("matrix.jl")
    include("readwrite.jl")
    include("sequence.jl")
    include("simulation.jl")
    include("substitution.jl")
    include("treeconsensus.jl")
    include("treedistance.jl")
    include("treemove.jl")
    include("treeplot.jl")
    include("treetest.jl")

end # module
