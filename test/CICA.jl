using JuMIT
using Distributions
using Base.Test
using MultivariateStats

nt = 2048
recv_n=2

uni=Uniform(-1, 1)

s1=rand(uni, nt)
s2=randn(nt)
unip = Uniform(-pi, pi)


s1 = s1 .* exp.(im*rand(unip,nt))
s2 = s2 .* exp.(im*rand(unip,nt))

S = hcat(s1, s2)';
A = complex.(randn(recv_n, 2), randn(recv_n, 2));
X = A * S;

ica=JuMIT.CICA.ICA(X, 2, nbins=10)
@time Sout=JuMIT.CICA.fastica!(ica, A=A);
err = Sout[2]
@test all(err.<1e-1)
##
# REAL ICA
s1=rand(uni, nt)
s2=randn(nt)
S = hcat(s1, s2)';
A = randn(recv_n, 2);
X = A * S;

ica=JuMIT.CICA.ICA(X, 2)
M = fit(ICA, X, 2; do_whiten=true, winit=ica.W,fun=icagfun(:gaus)); W = M.W;
@time Sout=JuMIT.CICA.fastica!(ica, A=A);

# using fastica from multivariate stats
err=JuMIT.CICA.error_unmixing(W',A)
#err = Sout[2]
@test (err<1e-1)


nbins=4
nt=1001;
splits = [round(Int, s) for s in linspace(1,nt,nbins+1)];
k=Array{UnitRange{Int64}}(nbins)
for ib in 1:nbins
           k[ib]=splits[ib]:splits[ib+1]
end
1:size(q,1), splits[idx]+1:splits[idx+1]



#####

function fix_scaling!(W, bins, magic_recv)
        nbins=length(bins)
        for ib in 1:nbins
                WW = transpose(W[:,:,ib]);
                A = inv(WW);
                A[:,1] /= A[magic_recv,1];
                A[:,2] /= A[magic_recv,2];
                W[:,:,ib] = transpose(inv(A));
        end
end


function error_unmixing(ica, A)
        nbins=length(ica.bins)
        Q=ica.Q
        W=ica.W
        for ib=1:nbins
                 QQ=view(Q,:,:,ib)
                 WW=view(W,:,:,ib)
                 SSE[ib] = error_scale_perm(WW'*QQ', A)
        end
end
