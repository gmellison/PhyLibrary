function tripleTest(triplefreq::Vector{Int64})
    x = sort(triplefreq)
    n = x[3]+x[2]
    major_pvalue = pvalue(BinomialTest(x[3],x[3]+x[2],0.5))/2
    minor_pvalue = pvalue(BinomialTest(x[1],x[1]+x[2],0.5))
    return major_pvalue, minor_pvalue
end