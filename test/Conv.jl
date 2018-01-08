using Revise
using JuMIT
using Base.Test

np=1000
n=777
n2=100
x = randn(n,n2); xa = similar(x)
z = zeros(np,n2);

# cover all lags for serial mode
@time for i in [0, 2, 4, n-1]
	JuMIT.Conv.pad_truncate!(x, z, n-i-1,i,np,1)
	JuMIT.Conv.pad_truncate!(xa, z, n-i-1,i,np,-1)
	@test x ≈ xa
end

np=1000
n=777
x = randn(n); xa = similar(x)
z = zeros(np);
@time for i in [0, 2, 4, n-1]
	JuMIT.Conv.pad_truncate!(x, z, n-i-1,i,np,1)
	JuMIT.Conv.pad_truncate!(xa, z, n-i-1,i,np,-1)
	@test x ≈ xa
end


# fast_filt
# note that this test works for only some delta functions
for nw in [7, 8], nr in [20, 10], ns in [11, 10]

	r = zeros(nr); w = zeros(nw); s = zeros(ns);
	ra = similar(r); wa = similar(w); sa = similar(s);

	r[5]=1.; w[3]=1.;

	JuMIT.Conv.mod!(s,r,w,:d)
	JuMIT.Conv.mod!(s,r,wa,:wav)
	@test w ≈ wa

	JuMIT.Conv.mod!(s,ra,w,:gf)
	@test ra ≈ r

	JuMIT.Conv.mod!(sa,ra,wa,:d)
	@test sa ≈ s
end



## dot product test for fast filt

function filt_loop(func, n2; )
	nwvec=[101, 100]
	nrvec=[1000, 1001]
	nsvec=[1500, 901]
	np2=2024;
	for nw in nwvec, nr in nrvec, ns in nsvec

		r=randn(nr, n2...)
		s=randn(ns, n2...)
		w=randn(nw, n2...)
		func(s,r,w,:d)

		sa=randn(ns, n2...);
		ra=similar(r)
		func(sa,ra,w,:gf)

		# dot product test
		@test dot(s, sa) ≈ dot(ra, r)

		r=randn(nr, n2...)
		s=randn(ns, n2...)
		w=zeros(nw, n2...)
		func(s,r,w,:wav)


		wa=randn(nw, n2...);
		ra=similar(r);
		func(s,ra,wa,:gf)

		# dot product test
		@test dot(r, ra) ≈ dot(wa, w)
	end
end

n2=128
@time filt_loop(JuMIT.Conv.mod!, n2)
@time filt_loop(JuMIT.Conv.mod!, n2)


using BenchmarkTools
# check if mod! will not have any allocations
n=1000
println("====================")
s=zeros(n,n)
r=zeros(n,n)
w=zeros(n,n)
pa=JuMIT.Conv.Param(d=s,gf=r,wav=w,dims=(n,), ntgf=n, ntd=n, ntwav=n);

@btime JuMIT.Conv.mod!(pa, :d)
@btime JuMIT.Conv.mod!(pa, :gf)
@btime JuMIT.Conv.mod!(pa, :wav)
