

nz=43; nx=35
ny=randn(nz, nx);

my1=J.Fdtd.get_rhovxI(ny)
my2=J.Fdtd.get_rhovzI(ny)


myp1 = randn(size(my1));
myp2 = randn(size(my2));

# spray
nyp=zeros(nz, nx);


J.Fdtd.grad_modrr_sprayrr!(nyp, myp1, myp2)
J.Fdtd.grad_modrr_sprayrrvx!(nyp, myp1)
J.Fdtd.grad_modrr_sprayrrvz!(nyp, myp2)

# dot product test
@test LinearAlgebra.dot(vcat(my1,my2), vcat(myp1,myp2)) ≈ LinearAlgebra.dot(ny, nyp)


