using JuMIT
using Base.Test

np=16;
n=7
x = randn(n); xa = similar(x)
z = complex.(zeros(np),zeros(np));

# cover all lags
for i in [0, 2, 4, n-1]
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

for nw in [101, 100], nr in [1000, 1001], ns in [1500, 901], nrp in [0, 500], nsp in [0, 501], nwp in [0, 50]

	r=randn(nr)
	s=randn(ns)
	w=randn(nw)
	JuMIT.DSP.fast_filt!(s,r,w,:s, nrplags=nrp, nsplags=nsp, nwplags=nwp)

	sa=randn(ns);
	ra=similar(r)
	JuMIT.DSP.fast_filt!(sa,ra,w,:r, nrplags=nrp, nsplags=nsp, nwplags=nwp)

	# dot product test
	@test dot(s, sa) ≈ dot(ra, r) 

	r=randn(nr)
	s=randn(ns)
	w=zeros(nw)
	JuMIT.DSP.fast_filt!(s,r,w,:w, nrplags=nrp, nsplags=nsp, nwplags=nwp)


	wa=randn(nw);
	ra=similar(r);
	JuMIT.DSP.fast_filt!(s,ra,wa,:r, nrplags=nrp, nsplags=nsp, nwplags=nwp)

	# dot product test
	@test dot(r, ra) ≈ dot(wa, w) 
end
