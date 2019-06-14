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


# Generic objective function evaluation
function function_eval(x, paE, x0, c_y, Cinv; layered=false, log_normal_param=false)
    nz = size(paE.mgrid[1],1)
    nx = size(paE.mgrid[2],1)
    
#     If layered model, then build 2D model
    if layered
        X = zeros(nz,nx)
    
        for i in 1:nz
            X[i,:] .= x[i]
        end
    else
        X = copy(x)
    end
    
#     If log normal parametrization, get physical model
    if log_normal_param
        X = exp.(X)
    end
    
#     If non log normal parametrization, nothing more to do
    if !log_normal_param
       
    end
    
    data_loss = 1 / c_y * func(X,paE)
        
    if layered
        prior_loss = (x .- x0)' * Cinv * (x .- x0)
    else
        prior_loss = (reshape(x,nz*nx) .- reshape(x0,nz*nx))' * Cinv * (reshape(x,nz*nx) .- reshape(x0,nz*nx))
    end
    
    return (data_loss + prior_loss)
end


# Generic gradient eval
function gradient_eval(x, grad_f, paE, x0, c_y, Cinv; layered=true, log_normal_param=false)
    nz = size(paE.mgrid[1],1)
    nx = size(paE.mgrid[2],1)
    
#     If layered model, then build 2D model
    if layered
        X = zeros(nz,nx)
    
        for i in 1:nz
            X[i,:] .= x[i]
        end
    else
        X = copy(x)
    end
    
#     If log normal parametrization, get physical model
    if log_normal_param
        X = exp.(X)
    end    
    
    mod!(paE,reshape(X,nz*nx),FGσ())
    
    grad_tot = 1 / c_y * reshape(paE.g, (nz * nx,1))[:,1]
    
    if log_normal_param
        grad_tot = grad_tot .* reshape(X, nz * nx)
    end
    
    if layered
        proj = zeros(nz, nz * nx)
    
        for i in 1:nz
            proj_l = zeros(nz, nx)
            proj_l[i,:] .= 1.

            proj[i,:] .= reshape(proj_l, nz * nx)
        end
        
        grad_tot = proj * grad_tot .+ 2 * Cinv * (x .- x0)
    else
        grad_tot = grad_tot .+ 2 * Cinv * (reshape(x,nz*nx) .- reshape(x0,nz*nx)) 
    end
    
    grad_f .= grad_tot
end
                

# Optimization routine
function solve_inverse_problem(paE, x0, c_y, Cinv; layered=true, log_normal_param=false)
    nz = size(paE.mgrid[1],1)
    nx = size(paE.mgrid[2],1)
    
    if layered
        x = copy(paE.σ[:,1])
    else
        x = copy(paE.σ)                    
    end
    
    f(x) = function_eval(x, paE, x0, c_y, Cinv; layered=layered, log_normal_param=log_normal_param) 
    gradf(x,grad_f) = gradient_eval(x, grad_f, paE, x0, c_y, Cinv; layered=layered, log_normal_param=log_normal_param)
                        
    function void_g(x, g)
    
    end

    function void_g_jac(x, mode, rows, cols, values)
    
    end
                        
    if log_normal_param
        prob = createProblem(nz, -10. .* ones(nz), 10. .* ones(nz), 0, 
                    Array{Float64}(undef,0), Array{Float64}(undef,0), 0, 0,
                    f, void_g, gradf, void_g_jac, nothing)
    else
        prob = createProblem(nz, 0. .* ones(nz), 100. .* ones(nz), 0, 
                    Array{Float64}(undef,0), Array{Float64}(undef,0), 0, 0,
                    f, void_g, gradf, void_g_jac, nothing)
    end
                            
    prob.x = copy(x0)
                            
    addOption(prob, "hessian_approximation", "limited-memory")
    addOption(prob, "max_iter", 30)
    addOption(prob, "limited_memory_max_history", 30)
                            
    solveProblem(prob)
                            
    return prob.x
end