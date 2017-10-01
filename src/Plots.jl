module Plots

using PyCall
using PyPlot
#import PyCall: @pyimport
#import PyPlot: pygui
#@pyimport matplotlib2tikz
import JuMIT.Interpolation
import JuMIT.Acquisition
import JuMIT.Grid
import JuMIT.Data
import JuMIT.Models

export printfig

if !PyPlot.isdisplayok()
pygui(false)
end

"""
save current fig using matlab2tikz
"""
function printfig(fig; filename::String = "FIG")
	if filename == "FIG"
		temp = 1;
		file = join([filename string(temp) ".tex"])
		while isfile(file)
			temp += 1;
			file = join([filename string(temp) ".tex"])
		end
		tikzfile = join([file[1:end-4] ".tikz"]);
	else
		file = filename;
		tikzfile = join([file[1:end-4] ".tikz"]);
	end			
	matplotlib2tikz.save(tikzfile,fig);
	fin = open(tikzfile,"r");
	fout = open(file,"w");
	lines = readlines(fin);
	close(fin);
	rm(tikzfile);
	lines_old = ["\\documentclass[tikz]{standalone}\n";
	"\\usepackage{pgfplots}\n";
	"\\pgfplotsset{compat=newest}\n";
	"\\usepackage{amsmath}\n";
	"\\usepackage{siunitx}\n";
	"\\begin{document}\n";lines;"\n\\end{document}\n"];

	lines = ["% Generated from Julia\n";"\\begin{center}\n";lines;"\\end{center}"]
	write(fout,lines)
	close(fout)
end

"""
Plot acquisition geometry `Acquisition.Geom` on 
and model grid `M2D`.

`attrib::Symbol=:unique` : default; plots unique source and receiver positions 
`ssvec::Vector{Int64}` : plot source and receivers of only these supersources
"""
function Geom(geom::Acquisition.Geom; ssvec=nothing, fields=[:s, :r])
	if(ssvec===nothing)
		if(:r ∈ fields)
			urpos = Acquisition.Geom_get([geom],:urpos)
			b = plot(urpos[2], urpos[1], "v", color="blue",ms=10)
		else
			b=nothing
		end
		if(:s ∈ fields)
			uspos = Acquisition.Geom_get([geom],:uspos)
			a = plot(uspos[2], uspos[1], "*", color="red",ms=15)
		else
			a=nothing
		end
	else
		if(:r ∈ fields)
			rxpos = [geom.rx[iss] for iss in ssvec]
			rzpos = [geom.rz[iss] for iss in ssvec]
			b = plot(vcat(rxpos...), vcat(rzpos...), "v", color="blue",ms=10)
		else
			b=nothing
		end
		if(:s ∈ fields)
			sxpos = [geom.sx[iss] for iss in ssvec]
			szpos = [geom.sz[iss] for iss in ssvec]
			a = plot(vcat(sxpos...), vcat(szpos...), "*", color="red",ms=15)
		else
			a=nothing
		end
	end
	# return handles for sources and receivers individually to modify later
	return [a,b]
end


"""
Plot the source wavelet used for acquisition.

# Arguments

* `acqsrc::Acquisition.Src` : source acquisition parameters
"""
function Src(acqsrc::Acquisition.Src; )
	wav = acqsrc.wav[1,1][:,:]
	powwav = (abs.(fft(wav)).^2)
	powwavdb = 10. * log10.(powwav./maximum(powwav)) # power in decibel after normalizing 
	tgrid= acqsrc.tgrid
	fgrid= Grid.M1D_fft(tgrid)
	subplot(2,1,1)
	plot(tgrid.x, wav)
	xlabel(L"$t$ (s)");
	subplot(2,1,2)
	plot(fgrid.x[1:div(fgrid.nx,2)], powwavdb[1:div(fgrid.nx,2)])
	ylabel("power (dB)");
	xlabel("frequency (Hz)");
	tight_layout()
end




"""
Plot time-domain data of type `Data.TD`

# Arguments
* `td::Vector{Data.TD}` : time-domain data to be compared

# Keyword Arguments
* `ssvec::Vector{Vector{Int64}}=fill([1], length(td))` : supersource vector to be plotted
* `fieldvec::Vector{Int64}=[1]` : field vector to be plotted
* `tr_flag::Bool=false` : plot time-reversed data when true
* `attrib::Symbol=:wav` : specify type of plot
"""
function TD(td::Vector{Data.TD}; ssvec::Vector{Vector{Int64}}=fill([1], length(td)), 
	    fields::Vector{Symbol}=[:P],
	    tr_flag::Bool=false, attrib::Symbol=:wav, 
	    wclip::Vector{Float64}=[maximum(broadcast(maximum, td[id].d)) for id in 1:length(td)],
	    bclip::Vector{Float64}=[minimum(broadcast(minimum, td[id].d)) for id in 1:length(td)],
	    )

	for id=1:length(td)

		fieldvec = findin(td[id].fields, fields)
		any(ssvec[id] .> td[id].acqgeom.nss) && error("invalid ssvec")
		ns = length(ssvec[id]);
		nr = maximum(td[id].acqgeom.nr);
		if(tr_flag)
			dp = hcat(td[id].d[ssvec[id],fieldvec][end:-1:1,:]...);
			extent=[1, nr*ns*length(fieldvec), td[id].tgrid.x[1],td[id].tgrid.x[end],]
		else
			dp = hcat(td[id].d[ssvec[id],fieldvec][:,:]...);
			extent=[1, nr*ns*length(fieldvec), td[id].tgrid.x[end],td[id].tgrid.x[1],]
		end
		if(attrib == :seis)
			subplot(1,length(td),id)
			imshow(dp, cmap="gray",aspect="auto", extent=extent, vmin=bclip[id], vmax=wclip[id])
			title("common shot gather")
			xlabel("receiver index");
			ylabel(L"$t$ (s)");
			colorbar();
			tight_layout()
		elseif(attrib == :wav)
			plot(td[id].tgrid.x,dp)
		else
			error("invalid attrib")
		end
	end
end


"""
Plot the velocity and density seismic models.

# Arguments

* `model::Models.Seismic` : model that should be plotted

# Keyword Arguments

* `xlim::Vector{Float64}=[model.mgrid.x[1],model.mgrid.x[end]]` : minimum and maximum limits of the second dimension while plotting
* `zlim::Vector{Float64}=[model.mgrid.z[1],model.mgrid.z[end]]` : minimum and maximum limits of the first dimension while plotting
* `fields::Vector{Symbol}=[:vp, :ρ]` : fields that are to be plotted, see Models.Seismic_get 
* `overlay_model=nothing` : use a overlay model
* `use_bounds=false` : impose bounds from the model or not?

"""
function Seismic(model::Models.Seismic; 
		 xlim::Vector{Float64}=[model.mgrid.x[1],model.mgrid.x[end]],
		 zlim::Vector{Float64}=[model.mgrid.z[1],model.mgrid.z[end]],
		 fields::Vector{Symbol}=[:vp, :ρ],
		 overlay_model=nothing, 
		 use_bounds=false,
		 )
	#indices
	ixmin = Interpolation.indminn(model.mgrid.x, xlim[1])[1]; ixmax = Interpolation.indminn(model.mgrid.x, xlim[2])[1]
	izmin = Interpolation.indminn(model.mgrid.z, zlim[1])[1]; izmax = Interpolation.indminn(model.mgrid.z, zlim[2])[1]
	nrow = (model.mgrid.nx > model.mgrid.nz) ? length(fields) : 1
	ncol = (model.mgrid.nx > model.mgrid.nz) ? 1 : length(fields)
	for i in 1:length(fields)
		(length(fields) ≠1) && subplot(nrow,ncol,i)

		f0 = Symbol((replace("$(fields[i])", "χ", "")),0)
		m = Models.Seismic_get(model, fields[i])[izmin:izmax,ixmin:ixmax]
		ext = [model.mgrid.x[ixmin], model.mgrid.x[ixmax], model.mgrid.z[izmax], model.mgrid.z[izmin],]
		vmin = use_bounds ? Models.Seismic_get(model, f0)[1] : minimum(m)
		vmax = use_bounds ? Models.Seismic_get(model, f0)[2] : maximum(m)
		ax = imshow(m,		     cmap="gray", 		     extent=ext, vmin=vmin, vmax=vmax)
			 xlabel(L"$x$ (m)");
			 ylabel(L"$z$ (m)");
			 title(string(fields[i]))
			 (length(fields) ≠1) && colorbar(ax, fraction=0.046, pad=0.04);
		
		if(!(overlay_model ===  nothing))
			# remove mean in the backgroud model
			mbg = Models.Seismic_get(overlay_model, fields[i])[izmin:izmax,ixmin:ixmax]
			mbg -= mean(mbg)
			mbgmax=maximum(abs, mbg)
			axbg = imshow(mbg,cmap="PiYG", alpha=0.3, vmin=-1.0*mbgmax, vmax=mbgmax, extent=ext)
		end
	 end
	 tight_layout()
	 subplots_adjust(top=0.88)

end

end # module
