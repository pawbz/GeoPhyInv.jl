using SIT
using Base.Test

nz=50
nx=50
npropwav=1
p = zeros(nz,nx,3,npropwav);
dpdz = zeros(nz,nx,3,npropwav);
dpdx = zeros(nz,nx,3,npropwav);

#p at p[t-1] v[t-3/2]
p[20,20,:,:] = 1.0 # point
p[3:nz-2, 3:nx-2, 1,npropwav] =  randn(nz-4, nx-4) # some random

modttI = randn(nz,nx)
modrrvz = randn(nz,nx)
modrrvx = randn(nz,nx)
deltax24I=2.34
deltaz24I=1.23
δt=2.3



# po at p[t] v[t-1/2]
po = copy(p)
SIT.Fdtd.advance!(po, dpdx, dpdz, modttI, modrrvx, modrrvz, deltax24I, deltaz24I, δt, nx, nz)


poo = copy(po)
# poo at p[t+1] v[t+1/2]
SIT.Fdtd.advance!(poo, dpdx, dpdz, modttI, modrrvx, modrrvz, deltax24I, deltaz24I, δt, nx, nz)

# pb at p[t] v[t+1/2], initial conditions for adjoint propagation
pb = zeros(p);
pb[:,:,1,1] = po[:,:,1,1]
pb[:,:,2:3,1] = -1.0.*poo[:,:,2:3,1]


# adjoint propagation  p[t-1] v[]
SIT.Fdtd.advance!(pb, dpdx, dpdz, modttI, modrrvx, modrrvz, deltax24I, deltaz24I, δt, nx, nz)


error = norm(pb[:,:,1,1]-p[:,:,1,1])./norm(p[:,:,1,1])
# desired accuracy? 
@test error<1e-6
