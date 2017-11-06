using JuMIT
using Distributions
using Base.Test
using MultivariateStats

nt = 10000
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

for nbins in [1,2,3,4]
        ica=JuMIT.CICA.ICA(X, 2, nbins=nbins)
        @time Sout=JuMIT.CICA.fastica!(ica, A=A);
        err = Sout[2]
        @test all(err.<1e-1)
        serr=JuMIT.CICA.error_after_ica(ica, S)
        @test all(serr.<0.05)
end


# REAL ICA
nt = 5000
s1=rand(uni, nt)
s2=randn(nt)
S = hcat(s1, s2)';
A = randn(recv_n, 2);
X = A * S;

ica=JuMIT.CICA.ICA(X, 2)
#M = fit(ICA, X, 2; do_whiten=true, winit=ica.W[:,:,1],fun=icagfun(:gaus)); W = M.W;
@time Sout=JuMIT.CICA.fastica!(ica, A=A);
err = Sout[2]
@test all(err.<1e-1)
# check error in the recovered components
serr=JuMIT.CICA.error_after_ica(ica, S)
@test all(serr.<0.05)

ica=JuMIT.CICA.ICA(X, 2, nbins=2)
#M = fit(ICA, X, 2; do_whiten=true, winit=ica.W[:,:,1],fun=icagfun(:gaus)); W = M.W;
@time Sout=JuMIT.CICA.fastica!(ica, A=A);
err = Sout[2]
@test all(err.<1e-1)
