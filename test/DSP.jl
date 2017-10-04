addprocs(4)
using JuMIT
using Base.Test

np=16000
n=7777
n2=100
x = randn(n,n2); xa = similar(x)
z = complex.(zeros(np,n2),zeros(np,n2));

# cover all lags for serial mode
@time for i in [0, 2, 4, n-1]
	JuMIT.DSP.nlag_npow2_pad_truncate!(x, z, n-i-1,i,np,1)
	JuMIT.DSP.nlag_npow2_pad_truncate!(xa, z, n-i-1,i,np,-1)
	@test x ≈ xa
end

# fast_filt
# note that this test works for only some delta functions
for nw in [7, 8], nr in [20, 10], ns in [11, 10]

	r = zeros(nr); w = zeros(nw); s = zeros(ns);
	ra = similar(r); wa = similar(w); sa = similar(s);

	r[5]=1.; w[3]=1.;

	JuMIT.DSP.fast_filt!(s,r,w,:s)
	JuMIT.DSP.fast_filt!(s,r,wa,:w)
	@test w ≈ wa

	JuMIT.DSP.fast_filt!(s,ra,w,:r)
	@test ra ≈ r

	JuMIT.DSP.fast_filt!(sa,ra,wa,:s)
	@test sa ≈ s
end



## dot product test for fast filt

function filt_loop(func, n2; inplaceflag=false)
	if(inplaceflag)
		nwvec=[1024,1024]
		nrvec=[1024,1024]
		nsvec=[1024,1024]
	else
		nwvec=[101, 100]
		nrvec=[1000, 1001]
		nsvec=[1500, 901]
	end
	np2=1024;
	for nw in nwvec, nr in nrvec, ns in nsvec, nrp in [0, 500], nsp in [0, 501], nwp in [0, 50]

		r=randn(nr, n2...)
		s=randn(ns, n2...)
		w=randn(nw, n2...)
		func(s,r,w,:s, nrplags=nrp, nsplags=nsp, nwplags=nwp, np2=np2)

		sa=randn(ns, n2...);
		ra=similar(r)
		func(sa,ra,w,:r, nrplags=nrp, nsplags=nsp, nwplags=nwp, np2=np2)

		# dot product test
		@test dot(s, sa) ≈ dot(ra, r)

		r=randn(nr, n2...)
		s=randn(ns, n2...)
		w=zeros(nw, n2...)
		func(s,r,w,:w, nrplags=nrp, nsplags=nsp, nwplags=nwp, np2=np2)


		wa=randn(nw, n2...);
		ra=similar(r);
		func(s,ra,wa,:r, nrplags=nrp, nsplags=nsp, nwplags=nwp, np2=np2)

		# dot product test
		@test dot(r, ra) ≈ dot(wa, w)
	end
end

n2=128
@time filt_loop(JuMIT.DSP.fast_filt!, n2)
@time filt_loop(JuMIT.DSP.fast_filt!, n2)
#@time filt_loop(JuMIT.DSP.fast_filt_parallel!, n2)
#@time filt_loop(JuMIT.DSP.fast_filt_parallel!, n2)

#n2=(128,4)
#@time filt_loop(JuMIT.DSP.fast_filt!, n2)
#@time filt_loop(JuMIT.DSP.fast_filt!, n2)
#@time filt_loop(JuMIT.DSP.fast_filt_parallel!, n2)
#@time filt_loop(JuMIT.DSP.fast_filt_parallel!, n2)
