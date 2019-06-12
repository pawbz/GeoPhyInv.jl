using Ipopt
using Statistics
using StatsBase


# Smoothing function for tampering
function smoothstep(x, N)
    if x >= 1
        return(1)
    end
    
    if x <= 0
        return(0)
    end
    
    f = 0

    for i in 0:N
        f = f .+ binomial(N + i, i) * binomial(2 * N + 1, N - i) * (- x) ^ i
    end
    
    f = x ^ (N+1) * f 
    
    if 0 < x < 1
        return(f)
    end
end
    
# Tampering of the properties
function tamp!(Q,n)
    nz = size(Q,1)
    nx = size(Q,2)
    
    for i in 1:nz
        for j in 1:nx
            lat = max( n - i, i - (nx - n) )
            vert = max( n - j, j - (nz - n) )
            
            dist = max(lat, vert)
                      
            Q[i,j] = Q[i,j] * smoothstep(1 - (dist+1) / n, 3)
        end
    end
end
