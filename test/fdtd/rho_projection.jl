

nz=43; nx=35
ny=randn(nz, nx);

my1=GeoPhyInv.get_rhovxI(ny)
my2=GeoPhyInv.get_rhovzI(ny)


myp1 = randn(size(my1));
myp2 = randn(size(my2));

# spray
nyp=zeros(nz, nx);


GeoPhyInv.grad_modrr_sprayrr!(nyp, myp1, myp2)
GeoPhyInv.grad_modrr_sprayrhovxI!(nyp, myp1)
GeoPhyInv.grad_modrr_sprayrhovzI!(nyp, myp2)

# dot product test
@test LinearAlgebra.dot(vcat(my1,my2), vcat(myp1,myp2)) ≈ LinearAlgebra.dot(ny, nyp)


