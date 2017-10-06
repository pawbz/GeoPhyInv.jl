using JuMIT
using Distributions
using Base.Test

nt = 1024
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
Sout = JuMIT.CICA.fastica!(X, 2,A=A)

err = Sout[4][end]

@test (err<1e-1)
