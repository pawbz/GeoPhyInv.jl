
#AGeom in mgrid

# To construct a variable to type `Medium`, the first step is to create a 2-D grid
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];


include("../../src/ageom.jl")

ageom=AGeom(mgrid)
@test (ageom ∈ mgrid)

update!(ageom[1], Recs(), [0,1], [10,20],)
update!(ageom[1], Srcs(), [0,1], [10,20],)
update!(ageom, SSrcs(), [0,1], [10,20], )
update!(ageom, SSrcs(), [0,1], [10,20], )


ageom2=vcat(ageom, ageom)


ageom2=deepcopy(ageom)
@test isequal(ageom,ageom2)


