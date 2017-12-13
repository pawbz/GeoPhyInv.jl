addprocs(4)
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
z = zeros(np));
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


# check memory for 2,D
for itimes in 1:4
	# check if fast_filt will not have any allocations
			println("====================")
	for n in [2, 2^2, 2^3]
		s=zeros(n,n)
		r=zeros(n,n)
		w=zeros(n,n)
		spow2=complex.(zeros(n,n))
		rpow2=complex.(zeros(n,n))
		wpow2=complex.(zeros(n,n))
		fftplan = plan_fft!(complex.(zeros(n,n)),[1],flags=FFTW.PATIENT)
		ifftplan = plan_ifft!(complex.(zeros(n,n)),[1],flags=FFTW.PATIENT)

		@time JuMIT.DSP.fast_filt_vec!(s,r,w,spow2,rpow2, wpow2,
			:s,0,n-1,0,n-1,0,n-1,n,fftplan, ifftplan)

	end
end

for itimes in 1:4
	# check if fast_filt will not have any allocations
			println("====================")
	for n in [16, 160, 1600]
		s=zeros(n)
		r=zeros(n)
		w=zeros(n)
		spow2=complex.(zeros(n))
		rpow2=complex.(zeros(n))
		wpow2=complex.(zeros(n))
		fftplan = plan_fft!(complex.(zeros(n)),[1],flags=FFTW.PATIENT)
		ifftplan = plan_ifft!(complex.(zeros(n)),[1],flags=FFTW.PATIENT)

		@time JuMIT.DSP.fast_filt_vec!(s,r,w,spow2,rpow2, wpow2,
			:s,0,n-1,0,n-1,0,n-1,n,fftplan, ifftplan)

	end
end
