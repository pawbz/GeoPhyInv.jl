
#Geom in mgrid

# To construct a variable to type `Medium`, the first step is to create a 2-D grid
mgrid = [range(0.0, stop=10.,step=0.1), range(0.0, stop=30.,step=0.2)];


include("../../src/geom.jl")

geom=Geom(mgrid)
@test (geom ∈ mgrid)

update!(geom[1], Recs(), [0,1], [10,20],)
update!(geom[1], Srcs(), [0,1], [10,20],)
update!(geom, SSrcs(), [0,1], [10,20], )
update!(geom, SSrcs(), [0,1], [10,20], )


geom2=vcat(geom, geom)


geom2=deepcopy(geom)
@test isequal(geom,geom2)


