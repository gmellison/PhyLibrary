function lower_tri(x::Matrix{Float64})
    (n,m) = size(x)
    if n != m 
        error("It must be a square matrix")
    end

    v = collect(Int, 1:(n*(n-1)/2))
    index = 1
    for i in 1:(n-1)
        for j in (i+1):n
        v[index] = (i-1)*n+j
        index = index + 1
        end
    end
    return v
end
  
function upper_tri(x::Matrix{Float64})
    (n,m) = size(x)
    if n != m 
        error("It must be a square matrix")
    end

    v = collect(Int, 1:(n*(n-1)/2))
    index = 1
    for i in 1:(n-1)
        for j in (i+1):n
        v[index] = (j-1)*n+i
        index = index + 1
        end
    end
    return v
end
  
function diag(x)
    dim = size(x)
    if length(dim) == 0
        x = floor(Int, x)
        if x<1
            error("the dimension must be > 1")
        end
        result = zeros(x,x)
        for i in 1:x
            result[i,i] = 1
        end
    elseif length(dim) == 1
        n = length(x)
        result = zeros(n,n)
        for i in 1:n
            result[i,i] = x[i]
        end
    elseif length(dim) == 2
        (n,m) = size(x)
        if n != m 
            error("It must be a square matrix")
        end
        result = collect(1:(n+1):(n^2))
    else
        error("the input must be a number, a vector, or a matrix")
    end
    return result
end
