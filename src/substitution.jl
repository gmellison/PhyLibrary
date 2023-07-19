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


function submodelDiMethyl_basefreq(a::Float64, b::Float64)

    return [b^2,2*a*b,a^2]./((a+b)^2)

end


function submodelDiMethyl_rateMatrix(a::Float64, b::Float64)

    Qmatrix = zeros(3, 3)
    Qmatrix[upper_tri(Qmatrix)] = [2*a,0,a] 
    Qmatrix[lower_tri(Qmatrix)] = [b,0,2*b] 
    Qmatrix[diag(Qmatrix)] = [-2*a, -(b+a), -2*b]

    basefreq = submodelDiMethyl_basefreq(a,b)

    #Dpi = zeros(3, 3)
    #Dpi[diag(Dpi)] = basefreq

    # rescale Qmatrix such that the diagonal sum is -1
    delta = -transpose(Qmatrix[diag(Qmatrix)]) * basefreq

    return Qmatrix/delta # should this return Qmatrix/delta * Dpi?
end



