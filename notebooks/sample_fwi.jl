### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 474f6970-9784-45fa-9d9b-bb1abbf71929
using Pkg; Pkg.activate("GeoPhyInv")

# ╔═╡ 7d9b2219-2074-4cab-8037-9391ee0406fd
Pkg.add("SegyIO"); using SegyIO

# ╔═╡ 75ad24f7-8b57-401b-8170-d3e35f3c960f
begin
 using PlutoLinks: @revise
 using PlutoUI, PlutoTest, Plots, CUDA
end

# ╔═╡ f86410b5-87f7-4761-91a8-856406d00234
using GeoPhyInv

# ╔═╡ bd68d9b1-633b-4c9b-b759-791581714306
using Statistics, LossFunctions, LinearAlgebra, MLUtils, FiniteDifferences, SparseArrays

# ╔═╡ 903d3dc9-b4d6-44aa-97db-5e36ccbb9f98
using DSP

# ╔═╡ 7c9cba9d-7e8f-446f-a664-2e8c988f2095
using CSV, DataFrames

# ╔═╡ bbfa8c8e-1ef5-459b-b5f5-edb15ff38fa8
TableOfContents()

# ╔═╡ 893cd92d-3d0d-4fb2-b465-f61a237fc154
# begin
#     Pkg.activate(Base.current_project())
#     Pkg.instantiate()
# end

# ╔═╡ 082b3d0f-bb73-4691-8442-d6dec1310cfd
# GeoPhyInv.set_preferences(ndims=2,use_gpu=true,datatype="Float32",order=2)

# ╔═╡ 16acf176-4ffa-40f0-a4ce-8d7d8483b5cb


# ╔═╡ c60027b8-ff8e-4671-aef3-36e136096b4d
begin
    # pa_mod = SeisForwExpt(FdtdAcoustic{FullWave}(:forward, 2), Homogeneous())
	d= 0.1
	mgrid1 = [range(000.0, stop =50.0, step=d), range(0.0, stop = 50, step=d)]
	newm=AcousticMedium(mgrid1)
	newm1=AcousticMedium(mgrid1)
    # pa_mod=
	# m1=GeoPhyInv.get_modelvector(pa_mod, [:invK]); 
    # @time update!(pa_mod, m1, [:invK])
    # m2=GeoPhyInv.get_modelvector(pa_mod, [:invK])
    # @test m1 ≈ m2
end

# ╔═╡ 46bdb8a0-350e-41e5-8637-a784b17a81e1
begin
	newm.vp=GeoPhyInv.vp(Float32.(100.0 .* ones(size(newm.vp.m))))
	# medium.vs=GeoPhyInv.vs(Float32.(vs1))
	newm.rho=GeoPhyInv.rho(Float32.(1900.0 .* ones(size(newm.vp.m))))
end

# ╔═╡ 33f69a08-639c-4d3f-930d-fe6743d01a3c
begin
	newm1.vp=GeoPhyInv.vp(Float32.(100.0 .* ones(size(newm.vp.m))))
	# medium.vs=GeoPhyInv.vs(Float32.(vs1))
	newm1.rho=GeoPhyInv.rho(Float32.(1900.0 .* ones(size(newm.vp.m))))
end

# ╔═╡ 09772196-4227-4a7b-80e2-072f3a31231c
GeoPhyInv.update!(newm, [:vp], rectangle = [[12,12], [14,48]], perc =-60)

# ╔═╡ 8d458409-0d0b-4e41-a0f2-21338377a9f3
begin
	nr=24
	function get_ageom(medium)
	    ageom=AGeom(medium.grid, :surf, SSrcs(1), Recs(nr)) # surface acquisition
	    GeoPhyInv.update!(ageom, Recs(nr), [0,10],[0,33])
	    # srcs_x =rand(Uniform(10,33))
	    srcs_x=10
	     # srcs_y = 5
	#     srcs_z = rand(Uniform(100,800))
	    srcs_z = 0
	    
	#     srcs_x=-29.709914448408938 
	#     srcs_y=-80.23644969547148
	#     println([srcs_x, srcs_y])
	# srcs_x = 0.0
	# srcs_y = 0.0
	# srcs_z=0
	GeoPhyInv.update!(ageom,SSrcs(1), [srcs_z,srcs_x], [srcs_z,srcs_x])
	end
	ageom=get_ageom(newm)
	ageom1=get_ageom(newm1)
	println(ageom)
end

# ╔═╡ 2fb4a69e-d0b9-4516-a14b-f96a9a37904d
begin
	heatmap(mgrid1[2],mgrid1[1],newm1[:vp].m,seriescolor=cgrad(:Blues),yflip=true)
	# ,clim=(200,400))
	Plots.scatter!(ageom1,SSrcs(1))
	Plots.scatter!(ageom1,Recs(24))
end

# ╔═╡ 4662c24a-24e1-4080-abab-0bc42ff7e39c


# ╔═╡ 0e932217-beef-4877-b803-70af59bac5dd
begin
	dt=5e-04
	itchunk=Int64(floor(2/dt))
	dt1=5e-4
	rf1=(1/dt)/(1/dt1)
end

# ╔═╡ dae779e0-3cce-4e72-9ddf-64fdc55858a4
begin
	tgrid=range(dt,stop=2,step=dt)
	wav = ricker(25.0, tgrid, tpeak =.5);
	wav1 = ricker(4.0, tgrid, tpeak =.5);
	# dt=step(tgrid)
	# println(dt)
	function get_srcwav(wav,i=10.0)
	    # 8 wavelengths along the diagonal of the medium, and total time of 1 diagonal
	
	    # nts_perc=i
	    # a=randn(floor(Int,length(tgrid)*nts_perc*0.01))
	    # wavv=DSP.conv(a,wav)[1:length(tgrid)]
	    wavv=wav
	    # wavv=randobs(read_data(selected_files[randobs(1:size(selected_files,1))]))[1:size(tgrid,1),1]
	    global freqmax=GeoPhyInv.Utils.findfreq(wavv,tgrid,attrib=:max)
	    global freqpeak=GeoPhyInv.Utils.findfreq(wavv,tgrid,attrib=:peak)
	    global freqmin=GeoPhyInv.Utils.findfreq(wavv,tgrid,attrib=:min)
	    tmax=maximum(tgrid); nt=length(tgrid); 
	#     println("nt:\t",length(tgrid))
	#     println("freqpeak:\t",freqpeak)
	
	
	    # intialize SrcWavelets
		# Srcs(tgrid, ageom, [:p, :vx])
	    srcwav = Srcs(tgrid, ageom, [:vz])
	    GeoPhyInv.update!(srcwav, [:vz], wavv)
	    return srcwav
	end
	srcwav=get_srcwav(wav)
	srcwav1=get_srcwav(wav1,2.0)
	plot(srcwav[1].d)
end

# ╔═╡ b0e1670a-2e83-4e3b-86b0-83425169cc3e
md"""
## Load Data
#### Scan a given folder
We provide two options, where the second option overrides the first.

We search of relevant files in a given input folder.

$(@bind mainfldrname TextField((50,2), default=".", placeholder="Enter the relative folder path to scan"))

Select the file extension to scan
$(@bind file_ext Select(["csv", "h5", "dat","sg2","sgy"], default="sgy"))

#### Forcefully load a particular file
Alternatively, we manually enter the path of the file to read. 
$(@bind forced_file FilePicker())
"""

# ╔═╡ 18e4bd74-bed3-4f58-8315-f6251b67d943
@bind selected_folder Select(readdir(mainfldrname),default="Synthetic_data")

# ╔═╡ 9a44176d-3bb7-4516-95b7-278e03fa3d4a
  function filesUI()
	if(forced_file===nothing)
	filtered_files=filter(readdir(joinpath(mainfldrname,selected_folder), join=true)) do f
	            (lowercase(splitext(f)[end]) == string(".",file_ext))
	end
	return md"""
	Depending on the choices made above, select files to preprocess and analyze. By default, 15 files will be randomly selected.
	$(@bind selected_files MultiCheckBox(filtered_files, select_all=true, default=sort(first(filtered_files,15)))))
	"""
	else
		return [forced_file]
	end
end

# ╔═╡ 074aad69-bca9-4866-9579-ddec53053eb5
filesUI()

# ╔═╡ 2e22d0a7-b0ae-4a22-8aea-439a102b473f
function read_file1(fname, file_ext)
if(file_ext == "sgy")
		return Float32.(segy_read(fname).data[:,1:24])
		end
if(file_ext == "h5")
		A=Array(first(first(h5open(fname))))[:,:,1]
		return A
end
	if(file_ext=="csv")
		b=CSV.read(fname,DataFrame)
		a=[]
		for i in 1:24
			if i==1
            a=b[:,i]
        else
            a=cat(a,b[:,i],dims=2)
        end
		end
		return a
	end
	# Similarly, add hdf5 and csv support later...
end

# ╔═╡ a192f25d-0305-4f82-809f-13f16c301084
Dobs=data=broadcast(selected_files) do fname
	read_file1(fname, file_ext)
end

# ╔═╡ 6bc5b99b-0f30-46c4-94b1-ce7c63cabdab
heatmap(Dobs[1][:, 2:24], clim=(-0.1, 0.1), yflip=true, ylim=(0, 800), c=:seismic)

# ╔═╡ cc05fa9a-c2d0-445a-ba4a-11530f07f5b3
heatmap(Dobs[5])

# ╔═╡ de9a4882-cba9-4707-a900-3d82b9a2faa1
begin
# 	# typically you will run this cell once
# 	# typically you will run this cell once
	tsnaps=tgrid[1:div(length(tgrid),20):end] # store 20 snapshots
	    rfields=[:vz]  # add as many fields are possible here to record
# 	pa_true1 = SeisForwExpt(
# 	     #FdtdAcou(),
# 	   FdtdAcoustic{FullWave}(:forward, 2),
# 	     snaps_field=:vz, # comment these if not testing
# 	     tsnaps=tsnaps, # comment these if not testing
# 	    # # pml_edges=[:xmin, :xmax, :zmax, :ymin, :ymax], # if you want surface waves to be modelled
# 	    pml_faces=[:zmax, :xmax, :xmin], # absorbing BC on edges
# 	    ageom = ageom,
	    
# 	    srcwav = srcwav,
# 	    medium = newm,
# 	    # rfields=rfields,
# 	    rfields=rfields,
	    
# 	#     sflags = 1,
# 	#     rflags = 1,
# 	    tgrid = tgrid,
# 	    stressfree_faces=[:zmin], # location of the stree-free boundary
# 	    verbose = true,
# 	);
end

# ╔═╡ 1f0eea14-c416-4bdf-83d8-908fc0d8649d
# begin
#     # pa_true = SeisForwExpt(FdtdAcoustic{FullWave}(:forward, 1), RandScatterer())
#     update!(pa_true1)
#     dobs = deepcopy(pa_true1.c.data[1])
# end

# ╔═╡ 2a737f8b-92d0-4042-96eb-6b85ea4030c7
# @time update!(m1, pa_mod, [:invK])

# ╔═╡ 000a9016-2b39-428b-a635-3accd26101a1
begin
# heatmap(pa_true1[:data][:vz],yflip=true, title="vz", c=:seismic)
# p2=heatmap(pa[:data][:vy], title="vy", c=:seismic)
# p3=heatmap(pa_mod[:data][:vx],yflip=true, title="vx", c=:seismic)
# plot(p1,p3, size=(1000,400), layout=(1,2))
end

# ╔═╡ a715feb4-8270-40a4-8796-da1feb6e972e
	pa_mod = SeisForwExpt(
	     #FdtdAcou(),
	   FdtdAcoustic{FullWave}(:forward, 2),
	     # snaps_field=:vz, # comment these if not testing
	     # tsnaps=tsnaps, # comment these if not testing
	    # # pml_edges=[:xmin, :xmax, :zmax, :ymin, :ymax], # if you want surface waves to be modelled
	    pml_faces=[:zmax, :xmax, :xmin], # absorbing BC on edges
	    ageom = ageom1,
	    srcwav = srcwav,
	    medium = newm1,
	    # rfields=rfields,
	    rfields=rfields,
	    
	#     sflags = 1,
	#     rflags = 1,
	    tgrid = tgrid,
	    stressfree_faces=[:zmin], # location of the stree-free boundary
	    verbose = true,
	);

# ╔═╡ a35aefc6-4de0-4886-8366-acda8b92e22e
update!(pa_mod)

# ╔═╡ 839dd5fb-e555-4872-b6f4-342f501a4e78
begin
	pa_true=deepcopy(pa_mod)
	update!(pa_true, newm)
	update!(pa_true)
end

# ╔═╡ c1751f48-7c26-4bac-8a47-77bf9c4b3832
plot(pa_true.c.data[1][1], 99)

# ╔═╡ a2618c5d-5c6a-4a73-b1a6-cbafb794f224
# pa_true.c.data[1][1][:vz].=data[1]

# ╔═╡ 0cf61f2f-435b-4fea-b6f7-e7458cb06b89
# update!(pa_true)

# ╔═╡ 78cd9487-b705-4a64-b53b-2caa13946003
begin
heatmap(pa_mod[:data][:vz],yflip=true, title="vz", c=:seismic)
# p2=heatmap(pa[:data][:vy], title="vy", c=:seismic)
# p3=heatmap(pa_mod[:data][:vx],yflip=true, title="vx", c=:seismic)
# plot(p1,p3, size=(1000,400), layout=(1,2))
end

# ╔═╡ df0ff9a7-d7ab-4d9d-a043-d8359d75b1da
plot(pa_mod[:data])

# ╔═╡ f49439df-3db3-48ae-bf9d-1387e3c14c8f
N = 10

# ╔═╡ 6d47ae1d-2554-4d42-af6c-02ce952bf14e
@bind mpara Select([:invK, :rho])

# ╔═╡ 5582110f-327c-419a-8759-ca11844c0c96
mpara

# ╔═╡ 51030579-4d25-4e03-ae49-38071dd58a40
dobs = deepcopy(pa_mod.c.data[1])

# ╔═╡ 9d566f59-fba0-4bd5-9945-49566bf74645
dobs[1][:vz].=Dobs[1]

# ╔═╡ 26c9f7a1-f18b-4fe3-a607-ad1f6a904bb8
pa_inv, pa_src = SeisInvExpt(pa_mod, dobs, [range(0,50,length=N), range(0,50,length=N)], [mpara])

# ╔═╡ 35d6fbbc-95a9-443b-924a-8f1946e6bcb8
heatmap(Dobs[1])

# ╔═╡ 0e5eb59a-aad6-4c84-a0e4-c8558964719f
# update!(pa_inv, pa_src)

# ╔═╡ e57e1a51-d69b-4459-9ae9-8869782a0721
# pa_inv

# ╔═╡ 1b102fdb-3e82-4b11-b100-fbe1b0588543
m2 = GeoPhyInv.get_modelvector(pa_inv)

# ╔═╡ 17ea5b1d-4296-47e0-a361-aaffb3c152f8
iszero(m2)

# ╔═╡ 26b174f6-3fbc-4f02-9279-bd83908b1aa7
LV(m)=GeoPhyInv.lossvalue(m,pa_inv)

# ╔═╡ d015a19a-a5bd-4660-ae22-dde2ae179313
Grad(g,m)=GeoPhyInv.gradient!(g, m, pa_inv)

# ╔═╡ ccbe14ab-5bc6-43ab-8e7b-14076cf37340
gtest=similar(m2)

# ╔═╡ cea31a54-6e9c-4505-ac40-8e850f98f9f1
GeoPhyInv.gradient!(gtest,m2,pa_inv)

# ╔═╡ 2d9177b5-c10d-4027-8b88-574863609031
# Grad(gtest,m2)

# ╔═╡ f470fc4b-4506-4039-a21f-2d4e9ae67e2a
LV(m2)

# ╔═╡ 406f05a3-4691-4405-a272-579d44aa3ff7
pa_inv

# ╔═╡ 7a956543-5a19-4912-958c-521b3ecc7390
length(pa_mod.c.exmedium.grid[1])

# ╔═╡ 2d0eda58-022e-4836-bde8-e4af663d20fb
pa_inv.P
# gtest

# ╔═╡ 6ec9f2b1-be56-411a-8f3c-518e80d9c078
heatmap(reshape(Array(pa_inv.mfull),length(pa_mod.c.exmedium.grid[1]),length(pa_mod.c.exmedium.grid[2])))

# ╔═╡ 654bf0b4-21d4-419c-bf19-c7909e6b2134
begin
	heatmap(mgrid1[2],mgrid1[1],newm[:vp].m,seriescolor=cgrad(:Blues),yflip=true)
	# ,clim=(200,400))
	Plots.scatter!(ageom,SSrcs())
	Plots.scatter!(ageom,Recs(24))
end

# ╔═╡ 52c54384-458b-410f-b393-10d475d4f13d
GeoPhyInv.get_modelvector(pa_inv)

# ╔═╡ 12a6565b-3edd-4e8b-8895-eb49366942ec
heatmap(pa_mod.c.exmedium.grid[2],pa_mod.c.exmedium.grid[1],reshape(Array(pa_inv.P'*GeoPhyInv.get_modelvector(pa_inv)),length(pa_mod.c.exmedium.grid[1]),length(pa_mod.c.exmedium.grid[2])),yflip=true,c=:seismic)

# ╔═╡ 30598e33-a52e-4c14-895e-7ceb45865573
heatmap(pa_mod.c.exmedium.grid[2],pa_mod.c.exmedium.grid[1],reshape(Array(pa_inv.P'*gtest),length(pa_mod.c.exmedium.grid[1]),length(pa_mod.c.exmedium.grid[2])),yflip=true,c=:seismic)

# ╔═╡ db243c16-e57e-4ea4-a291-eab647188023
heatmap(pa_mod.c.exmedium.grid[2],pa_mod.c.exmedium.grid[1],Array(deepcopy(pa_mod.c.gradients[:invK])),yflip=true,c=:seismic)

# ╔═╡ 28d9ed78-fd78-46f6-8fd1-13ebfedf1561
heatmap(pa_true.c.exmedium.grid[2],pa_true.c.exmedium.grid[1],Array(deepcopy(pa_true.c.gradients[:invK])),clim=(-1e-12,1e-12),yflip=true,c=:seismic)

# ╔═╡ a11ea26a-ab85-468b-809e-34144eb6578f
# using Optim

# ╔═╡ 6180ba9a-5aee-4e42-918c-5b6bb7728655
# using Ipopt,JuMP

# ╔═╡ 43fb3211-9ffc-4ea5-a459-539cf55ba009
# pa_inv = SeisInvExpt(pa_mod, dobs, [N, N], [mpara])

# ╔═╡ Cell order:
# ╠═bbfa8c8e-1ef5-459b-b5f5-edb15ff38fa8
# ╠═474f6970-9784-45fa-9d9b-bb1abbf71929
# ╠═75ad24f7-8b57-401b-8170-d3e35f3c960f
# ╠═893cd92d-3d0d-4fb2-b465-f61a237fc154
# ╠═f86410b5-87f7-4761-91a8-856406d00234
# ╠═7d9b2219-2074-4cab-8037-9391ee0406fd
# ╠═082b3d0f-bb73-4691-8442-d6dec1310cfd
# ╠═bd68d9b1-633b-4c9b-b759-791581714306
# ╠═16acf176-4ffa-40f0-a4ce-8d7d8483b5cb
# ╠═c60027b8-ff8e-4671-aef3-36e136096b4d
# ╠═46bdb8a0-350e-41e5-8637-a784b17a81e1
# ╠═33f69a08-639c-4d3f-930d-fe6743d01a3c
# ╠═09772196-4227-4a7b-80e2-072f3a31231c
# ╠═8d458409-0d0b-4e41-a0f2-21338377a9f3
# ╠═2fb4a69e-d0b9-4516-a14b-f96a9a37904d
# ╠═4662c24a-24e1-4080-abab-0bc42ff7e39c
# ╠═0e932217-beef-4877-b803-70af59bac5dd
# ╠═903d3dc9-b4d6-44aa-97db-5e36ccbb9f98
# ╠═dae779e0-3cce-4e72-9ddf-64fdc55858a4
# ╠═b0e1670a-2e83-4e3b-86b0-83425169cc3e
# ╠═18e4bd74-bed3-4f58-8315-f6251b67d943
# ╠═074aad69-bca9-4866-9579-ddec53053eb5
# ╠═9a44176d-3bb7-4516-95b7-278e03fa3d4a
# ╠═2e22d0a7-b0ae-4a22-8aea-439a102b473f
# ╠═7c9cba9d-7e8f-446f-a664-2e8c988f2095
# ╠═a192f25d-0305-4f82-809f-13f16c301084
# ╠═6bc5b99b-0f30-46c4-94b1-ce7c63cabdab
# ╠═cc05fa9a-c2d0-445a-ba4a-11530f07f5b3
# ╠═de9a4882-cba9-4707-a900-3d82b9a2faa1
# ╠═1f0eea14-c416-4bdf-83d8-908fc0d8649d
# ╠═2a737f8b-92d0-4042-96eb-6b85ea4030c7
# ╠═000a9016-2b39-428b-a635-3accd26101a1
# ╠═a715feb4-8270-40a4-8796-da1feb6e972e
# ╠═a35aefc6-4de0-4886-8366-acda8b92e22e
# ╠═839dd5fb-e555-4872-b6f4-342f501a4e78
# ╠═c1751f48-7c26-4bac-8a47-77bf9c4b3832
# ╠═a2618c5d-5c6a-4a73-b1a6-cbafb794f224
# ╠═0cf61f2f-435b-4fea-b6f7-e7458cb06b89
# ╠═78cd9487-b705-4a64-b53b-2caa13946003
# ╠═df0ff9a7-d7ab-4d9d-a043-d8359d75b1da
# ╠═f49439df-3db3-48ae-bf9d-1387e3c14c8f
# ╠═6d47ae1d-2554-4d42-af6c-02ce952bf14e
# ╠═5582110f-327c-419a-8759-ca11844c0c96
# ╠═51030579-4d25-4e03-ae49-38071dd58a40
# ╠═9d566f59-fba0-4bd5-9945-49566bf74645
# ╠═26c9f7a1-f18b-4fe3-a607-ad1f6a904bb8
# ╠═35d6fbbc-95a9-443b-924a-8f1946e6bcb8
# ╠═0e5eb59a-aad6-4c84-a0e4-c8558964719f
# ╠═e57e1a51-d69b-4459-9ae9-8869782a0721
# ╠═1b102fdb-3e82-4b11-b100-fbe1b0588543
# ╠═17ea5b1d-4296-47e0-a361-aaffb3c152f8
# ╠═26b174f6-3fbc-4f02-9279-bd83908b1aa7
# ╠═d015a19a-a5bd-4660-ae22-dde2ae179313
# ╠═ccbe14ab-5bc6-43ab-8e7b-14076cf37340
# ╠═cea31a54-6e9c-4505-ac40-8e850f98f9f1
# ╠═2d9177b5-c10d-4027-8b88-574863609031
# ╠═f470fc4b-4506-4039-a21f-2d4e9ae67e2a
# ╠═406f05a3-4691-4405-a272-579d44aa3ff7
# ╠═7a956543-5a19-4912-958c-521b3ecc7390
# ╠═2d0eda58-022e-4836-bde8-e4af663d20fb
# ╠═6ec9f2b1-be56-411a-8f3c-518e80d9c078
# ╠═654bf0b4-21d4-419c-bf19-c7909e6b2134
# ╠═52c54384-458b-410f-b393-10d475d4f13d
# ╠═12a6565b-3edd-4e8b-8895-eb49366942ec
# ╠═30598e33-a52e-4c14-895e-7ceb45865573
# ╠═db243c16-e57e-4ea4-a291-eab647188023
# ╠═28d9ed78-fd78-46f6-8fd1-13ebfedf1561
# ╠═a11ea26a-ab85-468b-809e-34144eb6578f
# ╠═6180ba9a-5aee-4e42-918c-5b6bb7728655
# ╠═43fb3211-9ffc-4ea5-a459-539cf55ba009
