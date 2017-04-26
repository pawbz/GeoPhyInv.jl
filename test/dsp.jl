using SIT
using Base.Test

np=16;
n=7
x = randn(n); xa = similar(x)
z = complex(zeros(np),zeros(np));

SIT.DSP.nlag_npow2_pad_truncate!(x, z, n-3,2,np,1)

SIT.DSP.nlag_npow2_pad_truncate!(xa, z, n-3,2,np,-1)

@test_approx_eq x xa




# fast_filt
r = zeros(10); w = zeros(7); s = zeros(10);
ra = similar(r); wa = similar(w); sa = similar(s);

r[5]=1.; w[3]=1.; 
SIT.DSP.fast_filt!(s,r,w,:s)

SIT.DSP.fast_filt!(s,r,wa,:w)
@test_approx_eq w wa

SIT.DSP.fast_filt!(s,ra,w,:r)
@test_approx_eq ra r

SIT.DSP.fast_filt!(sa,ra,wa,:s)
@test_approx_eq sa s


