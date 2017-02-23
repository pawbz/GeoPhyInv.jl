using MultivariateStats
using PyPlot
using Base.Test

n = 1000
m = 8
noise_perc = 0.01

t = linspace(0.0, 60.0, n)
s1 = sin(t * 2)
s2 = 1.0 - 2.0 * Bool[isodd(floor(Int, _ / 3)) for _ in t]
s3 = Float64[mod(_, 5.0) for _ in t]

s1 = s1 - mean(s1);
s2 = s2 - mean(s2);
s3 = s3 - mean(s3);

s1 += noise_perc * randn(n)
s2 += noise_perc * randn(n)
s3 += noise_perc * randn(n)
# just random vectors
s4 =  randn(n)
s5 =  randn(n)



##

println(mean(s2.*s1))

plot(s5); show()
close("all")
println(mean(s1)*mean(s2))

S = hcat(s1, s4)'

k = size(S,1)
figure(4);
plot(S')
title("three source signals")
show()


A = randn(m, k)

X = A * S
mv = vec(mean(X,2))

figure(2)
plot(X'); show()

M = fit(ICA, X, k; do_whiten=true)
figure(3)
W = M.W;
λ = transform(M,X);
plot(λ'); show()



@test_approx_eq transform(M, X) W' * (X .- mv)
E = W' * (X .- mv)

#@test_approx_eq W'W eye(k)

print(mean(abs((s4.-vec(E[2,:])))))


##

close("all")
