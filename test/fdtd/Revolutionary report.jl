### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 75e3051c-1273-11ed-245c-07301bfa1384
import Pkg; Pkg.activate("GeoPhyInv")


# ╔═╡ eac68826-f6e8-4376-9dc5-654c87666703
begin
	using Revise
	using GeoPhyInv
	using Statistics
	using Test
	using PlutoUI
	using Plots; gr()
end

# ╔═╡ 6b38bcf3-227f-4da8-8ac5-eb582d21269f
@init_parallel_stencil(2, true, Float32, 2)

# ╔═╡ aaccf46a-37eb-4d5c-b63a-600dfdf007d3
GeoPhyInv.Data.Array{2}

# ╔═╡ fd10cc83-86ce-4ab4-b0a1-dedcbefc1ceb
c=GeoPhyInv.get_boundary_store(GeoPhyInv.vx(), GeoPhyInv.FdtdAcoustic(), 10, 10, 7)

# ╔═╡ 12bd3624-9630-4145-89e2-44e3274f0d7a
names(c)

# ╔═╡ 57a5ded5-a17a-4bf0-b488-14f1d2f3f51a
length(c[:x]), length(c[:snap])

# ╔═╡ 7fc84cea-1ee7-45dc-95c7-4ac647f52a02
begin
	a=GeoPhyInv.Fields(GeoPhyInv.FdtdAcoustic())
	b=GeoPhyInv.Fields(GeoPhyInv.FdtdAcoustic(), "d")
	filter(x->x ∉ b, a)
end

# ╔═╡ 3fe9ba0f-2cfe-4f60-9d0f-7d314c64c195
begin
	medium = Medium(:acou_homo2D, 5); # load a simple homogeneous acoustic medium from the gallery
update!(medium, [:vp, :rho], randn_perc = 5); # add some random noise to the medium
println(medium)

	# ### AGeom
ageom = AGeom(medium.grid, :xwell, SSrcs(2)); # load a simple acquisition using `mgrid` of the medium
println(ageom)


	tgrid = range(0.0, stop = 2.0, length = 2500); # generate a time grid
wav = ricker(15.0, tgrid, tpeak = 0.25); # ricker wavelet
srcwav = Srcs(tgrid, ageom, [:p]);
update!(srcwav, [:p], wav);
	
	tsnaps=tgrid[1:div(length(tgrid),20):end] 
end


# ╔═╡ 5d9025fe-5cde-455b-b377-cbf5c7f3519e
pa = SeisForwExpt(
    FdtdAcoustic(),
    medium = medium,
	backprop_flag=1,
    ageom = ageom,
    srcwav = srcwav,
    snaps_field = :p,
    tsnaps = tsnaps,
    tgrid = tgrid,
    rfields = [:p],
    verbose = true,
);


# ╔═╡ 14580561-c81f-4c23-8daf-5af9faddc766
@time update!(pa);

# ╔═╡ de900b6a-1087-46e7-b8e8-7021e5c7b7f1
begin
	
		# change source flag and update wavelets in pa
			pa.c.sflags=[-2];
			GeoPhyInv.update!(pa,[srcwav])
			pa.c.backprop_flag=-1 # do backpropagation
end

# ╔═╡ 296434f0-4edd-4331-bed4-f19beabb7d8e
@time update!(pa);

# ╔═╡ 97c871a3-75a0-47e7-9acf-22eca2a021c6
# ### Plotting pressure snapshots after acoustic simulation
function plot_snapshots(name)
  clip_perc = 90
   pmax = maximum([maximum(pa[:snaps][it]) for it = 1:length(tsnaps)])
   pmax -= pmax * 0.01 * clip_perc

   anim = @animate for it = 1:length(tsnaps)
      plot(
           [
               (
                 f = pa[:snaps, is][it];
                   heatmap(
                       f,
                        aspect_ratio = 1,
                       xlims = (1, size(f, 2)),
                       ylims = (1, size(f, 1)),
                         yflip = true,
                        seriescolor = cgrad(:seismic),
                        clims = (-pmax, pmax),
                         legend = false,
                     )
                ) for is = 1:2
          ]...,
            title = string("snapshot at: ", round(tsnaps[it], digits = 5), " s"),
        )
     end
     return gif(anim, string(name,".gif"), fps = 1)
 end
 # plot_snapshots("acoustic_snaps")

# ╔═╡ 5dedbdc0-3a5c-4835-934f-9d61c36ecbfb
plot_snapshots("acoustic_snaps")

# ╔═╡ 8ce4180d-f2df-4a1a-a3b3-ae83125b6eb1
pa[:data].d[1]

# ╔═╡ Cell order:
# ╠═75e3051c-1273-11ed-245c-07301bfa1384
# ╠═eac68826-f6e8-4376-9dc5-654c87666703
# ╠═6b38bcf3-227f-4da8-8ac5-eb582d21269f
# ╠═aaccf46a-37eb-4d5c-b63a-600dfdf007d3
# ╠═fd10cc83-86ce-4ab4-b0a1-dedcbefc1ceb
# ╠═12bd3624-9630-4145-89e2-44e3274f0d7a
# ╠═57a5ded5-a17a-4bf0-b488-14f1d2f3f51a
# ╠═7fc84cea-1ee7-45dc-95c7-4ac647f52a02
# ╠═3fe9ba0f-2cfe-4f60-9d0f-7d314c64c195
# ╠═5d9025fe-5cde-455b-b377-cbf5c7f3519e
# ╠═14580561-c81f-4c23-8daf-5af9faddc766
# ╠═de900b6a-1087-46e7-b8e8-7021e5c7b7f1
# ╠═296434f0-4edd-4331-bed4-f19beabb7d8e
# ╠═97c871a3-75a0-47e7-9acf-22eca2a021c6
# ╠═5dedbdc0-3a5c-4835-934f-9d61c36ecbfb
# ╠═8ce4180d-f2df-4a1a-a3b3-ae83125b6eb1
