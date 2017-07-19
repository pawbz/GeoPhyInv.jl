module Plots

using PyCall
using PyPlot
#import PyCall: @pyimport
#import PyPlot: pygui
#@pyimport matplotlib2tikz
import SIT.Acquisition
import SIT.Grid
import SIT.Data
import SIT.Models

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
while is		file(file)
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
"""
function Geom(geom::Acquisition.Geom;
	      iss::Int64=0
	     )
	if(iss==0)
		urpos = Acquisition.Geom_get([geom],:urpos)
		uspos = Acquisition.Geom_get([geom],:uspos)

		plot(urpos[2], urpos[1], "v", color="blue",ms=10)
		plot(uspos[2], uspos[1], "*", color="red",ms=15)
	else
		plot(geom.rx[iss], geom.rz[iss], "v", color="blue",ms=10)
		plot(geom.sx[iss], geom.sz[iss], "*", color="red",ms=15)
	end
end


"""
Plot the source wavelet used for acquisition.

# Arguments

* `acqsrc::Acquisition.Src` : source acquisition parameters
"""
function Src(acqsrc::Acquisition.Src)
	plot(acqsrc.tgrid.x, acqsrc.wav[1,1][:,:])
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
	    fieldvec::Vector{Int64}=[1],
	    tr_flag::Bool=false, attrib::Symbol=:wav)
	nfield = length(fieldvec);

	for id=1:length(td)
		any(ssvec[id] .> td[id].acqgeom.nss) && error("invalid ssvec")
		ns = length(ssvec[id]);
		nr = maximum(td[id].acqgeom.nr);
		if(tr_flag)
			dp = hcat(td[id].d[ssvec[id],fieldvec][end:-1:1,:]...);
			extent=[1, nr*ns*nfield, td[id].tgrid.x[1],td[id].tgrid.x[end],]
		else
			dp = hcat(td[id].d[ssvec[id],fieldvec][:,:]...);
			extent=[1, nr*ns*nfield, td[id].tgrid.x[end],td[id].tgrid.x[1],]
		end
		if(attrib == :seis)
			subplot(1,length(td),id)
			imshow(dp, cmap="gray",aspect="auto", extent=extent)
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

"""
function Seismic(model::Models.Seismic; 
		 xlim::Vector{Float64}=[model.mgrid.x[1],model.mgrid.x[end]],
		 zlim::Vector{Float64}=[model.mgrid.z[1],model.mgrid.z[end]] 
		 )
	#indices
	ixmin = findfirst(model.mgrid.x, xlim[1]); ixmax = findfirst(model.mgrid.x, xlim[2])
	izmin = findfirst(model.mgrid.z, zlim[1]); izmax = findfirst(model.mgrid.z, zlim[2])
	nrow = (model.mgrid.nx > model.mgrid.nz) ? 2 : 1
	ncol = (model.mgrid.nx > model.mgrid.nz) ? 1 : 2
	subplot(nrow,ncol,1)
	ax = imshow(Models.χ(model.χvp[izmin:izmax,ixmin:ixmax],model.vp0,-1), 
	     cmap="gray",
	         extent=[model.mgrid.x[ixmin], 
		  model.mgrid.x[ixmax], model.mgrid.z[izmax], model.mgrid.z[izmin],]);
		 xlabel(L"$x$ (m)");
		 ylabel(L"$z$ (m)");
		 title(L"$v_p$")
		 colorbar();
	subplot(nrow,ncol,2)
	ax = imshow(Models.χ(model.χρ[izmin:izmax,ixmin:ixmax],model.ρ0,-1), 
	     cmap="gray",
	         extent=[model.mgrid.x[ixmin], 
		  model.mgrid.x[ixmax], model.mgrid.z[izmax], model.mgrid.z[izmin],]);
		 xlabel(L"$x$ (m)");
		 ylabel(L"$z$ (m)");
		 title(L"$\rho$")
		 colorbar();
	 tight_layout()

end

end # module
