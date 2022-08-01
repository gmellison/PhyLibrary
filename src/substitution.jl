###########################################################################
# substitution models
###########################################################################
function submodelGTR_rateMatrix(basefreq::Vector{Float64}, rates::Vector{Float64})
    nstates = length(basefreq)
    Qmatrix = zeros(nstates, nstates)
    Qmatrix[upper_tri(Qmatrix)] = rates
    Qmatrix = Qmatrix+transpose(Qmatrix)
    for i in 1:nstates
        Qmatrix[:,i] = Qmatrix[:,i] .* basefreq[i]
    end

    scale = 0
    for i in 1:nstates
        Qmatrix[i,i] = -sum(Qmatrix[i,:])
        scale += Qmatrix[i,i] 
    end

    # rescale Qmatrix such that the diagonal sum is -1
    delta = - transpose(Qmatrix[diag(Qmatrix)]) * basefreq

    return Qmatrix/delta
end