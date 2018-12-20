# World of testing
#Testing the forward-mode of the Poisson solver module Poisson.jl using combinations of random eigenfunctions as a source
#Developed by:
#Niels Grobbe, Massachusetts Institute of Technology, USA.
#In collaboration with: Aimé Fournier & Laurent Demanet, Massachusetts Institute of Technology, USA.
#Date: October, 2017
#Contact: ngrobbe@gmail.com


using JuMIT
using Polynomials
using Test
using StatsBase


############################################    
# for manual random or uniform sigma testing
# sigma=ones(rand(16:64),rand(16:64)); # random size
# # sigma=ones(40,40);    
############################################
    
# function input random source alterations and sigma
#ntst = 8;#  how many tests
ntst=32; # how many tests: larger number, more more gridpoints (finer grid)
wavs = [rand(1:8,2,5); 8*randn(1,5)]; # pre-set random wavenumbers --> how many complex source terms we want to add..
t = ["Computed potential" "True potential" "interior" "boundary"];
m=["or","vk"];

err=zeros(ntst,4);  
for i = 1 : ntst
    sigma=ones(floor(Int,2^(4+2.5*(i - 1)/(ntst - 3))),
	       floor(Int,2^(4+2.5*(i - 1)/(ntst - 3))));

    global fields # to make it a global variable

# Get model size (also in module, but here local for creating the input for the module)
    nz=size(sigma,1);
    nx=size(sigma,2);
    dz= 1.0*inv(nz - 1);
    dx= 1.0*inv(nx - 1);

    z=(0 : nz - 1)*dz;
    x=(0:nx-1)*dx;

# For random tests
    if size(wavs, 1) == 3
   	 l = wavs[1 : 2,:];
   	 j = wavs[3,:]'; # why 3? (take 3rd row which is amplitudes..) # --> make row vector for c definition below
    else
# for just a single test (if random tests not specified)
    	l = [rand(1:ceil(nz/3),1,3) ;rand(1:ceil(nx/3),1,3)];  #% random wavevector ... #% for eigenfunction source
    	j = 8*randn(1, size(l, 2));
    end

    l = pi*l; 
# preallocate source and reference electrical potential
    src = zeros(nz,nx);
    psi_ref = src;

    for k = 1 : size(l,2)
        c = j[k].*cos.(l[1,k].*z)*cos.(l[2,k].*x'); # result is size nz*nx
        src = src + c; # add eigenfunction to source (to make source more complex)
        psi_ref = psi_ref - c./(l[:,k]'*l[:,k]);
    end

    psi_ref = psi_ref./sigma;
println("Solving Poisson equation, testrun number:",i)

    psi = JuMIT.Poisson.solve(src,sigma,1) # call forward solver

# subtract mean for DC
    psi_ref=psi_ref.-StatsBase.mean(psi_ref);
    psi=psi.-StatsBase.mean(psi);

    fields=(psi,psi_ref); #make tuple

#println("testing...")
# Normalized RSME errors for interior points only, and interior+boundary points
    err[i,1] = size(fields[1],1);
    err[i,2] = size(fields[1],2);
    err[i,3] = sqrt(sum(sum((fields[1][2 : end - 1, 2 : end - 1] - # diff recovered-true only for inner points without boundaries
                fields[2][2 : end - 1, 2 : end - 1]).^2))/sum(sum(fields[2][2 : end - 1, 2 : end - 1].^2)));
    err[i, 4] = sqrt(sum(sum([(fields[1][1,1 : end - 1]-
                fields[2][1,1 : end - 1])' fields[1][1 : end - 1,end]'-
                fields[2][1 : end - 1,end]' (fields[1][end,end : -1 : 2] -
                fields[2][end,end : -1 : 2])' fields[1][end : -1: 2,1]' - fields[2][end : -1: 2,1]'].^2))/
    sum(sum([fields[2][1,1 : end - 1]' fields[2][1 : end - 1,end]' fields[2][end,end : -1 : 2]' fields[2][end : -1: 2,1]'].^2)));

# test to study error decay for interior only, and interior+boundary points
 if i > 1
        x = minimum(err[:,1 : 2], dims=2); # error over smallest dimension
      for k = 3:4
           global p,chk
           p = polyfit(log.(x[1:i]), log.(err[1:i,k]),1);  # fit polynomial of first degree using least squares fit
#           println(p)
           chk=round(p[1],4)
       end
   end

end # end i-loop

# Test passed?

# first pass interior points, then also include boundary points
@test chk<=-0.8
@test (err[end,3]/err[1,3]) < 0.2
@test (err[end,4]/err[1,4]) < 0.2 

#if(chk<=-0.8 && err[end,3]/err[1,3] < 0.2 && err[end,4]/err[1,4] < 0.2 ) i
#println("tests successfully passed!")
#println("log-log slope <=-0.8? slope=",chk)
#println("Normalized RMSE interior < 0.2? err[end,3]/err[1,3]=",err[end,3]/err[1,3])
#println("Normalized RMSE int.+boundary < 0.2? err[end,4]/err[1,4]=",err[end,4]/err[1,4])
#
#println("log-log slope 'chk' <=-0.8? slope=",chk)
#println("iteration=",collect(1:length(err[:,3])));
#println("Normalized RMSE interior < 0.5? err[:,3]=",err[:,3])
#println("Normalized RMSE int.+boundary points < 0.5? err[:,4]=",err[:,4])
#error("ERROR: RMSE and/or log-log slope too large: try another run, or try increasing amount of gridpoints via ntst")
