module Plots

using PyCall, PyPlot
import PyCall: @pyimport
import PyPlot: pygui
@pyimport matplotlib2tikz
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
		urpos = Acquisition.Geom_get(geom,:urpos)
		uspos = Acquisition.Geom_get(geom,:uspos)

		plot(urpos[2], urpos[1], "v", color="blue",ms=10)
		plot(uspos[2], uspos[1], "*", color="red",ms=15)
	else
		plot(geom.rx[iss], geom.rz[iss], "v", color="blue",ms=10)
		plot(geom.sx[iss], geom.sz[iss], "*", color="red",ms=15)
	end
end


"""
Plot acqsrc
"""
function Src(acqsrc::Acquisition.Src, attrib::Symbol)
	
	plot(acqsrc.tgrid.x, acqsrc.wav[:,1,1])

end




"""
Plot time-domain data of type `Data.TD`

# Arguments
* `td::Data.TD` : 

# Keyword Arguments
* `ssvec::Vector{Int64}=[1]` : supersource vector to be plotted
* `fieldvec::Vector{Int64}=[1]` : field vector to be plotted
* `tr_flag::Bool=false` : plot time-reversed data when true
* `attrib::Symbol=:wav` : specify type of plot
"""
function TD(td::Data.TD; ssvec::Vector{Int64}=[1], fieldvec::Vector{Int64}=[1],
	    tr_flag::Bool=false, attrib::Symbol=:wav)
	
	nfield = length(fieldvec);
	ns = length(ssvec);
	nr = maximum(td.acqgeom.nr);

	if(tr_flag)
		dp = reshape(hcat(td.d[ssvec,fieldvec][end:-1:1,:]...),td.tgrid.nx,nr*ns*nfield);
		extent=[1, nr*ns*nfield, td.tgrid.x[1],td.tgrid.x[end],]
	else
		dp = reshape(hcat(td.d[ssvec,fieldvec][:,:]...),td.tgrid.nx,nr*ns*nfield);
		extent=[1, nr*ns*nfield, td.tgrid.x[end],td.tgrid.x[1],]
	end
	if(attrib == :seis)
		imshow(dp, cmap="gray",aspect="auto", extent=extent)

		xlabel("receiver index");
		ylabel(L"$t$ (s)");
		colorbar();
		tight_layout()
	elseif(attrib == :wav)
		plot(td.tgrid.x,dp)
	else
		error("invalid attrib")
	end



end


"""
Plot seismic model
"""
function Seismic(model::Models.Seismic)
	subplot(121)
	ax = imshow(Models.χ(model.χvp,model.vp0,-1), cmap="gray",
	         extent=[model.mgrid.x[1], model.mgrid.x[end], model.mgrid.z[end], model.mgrid.z[1],]);
		 xlabel(L"$x$ (m)");
		 ylabel(L"$z$ (m)");
		 colorbar();
	subplot(122)
	ax = imshow(Models.χ(model.χρ,model.ρ0,-1), cmap="gray",
	         extent=[model.mgrid.x[1], model.mgrid.x[end], model.mgrid.z[end], model.mgrid.z[1],]);
		 xlabel(L"$x$ (m)");
		 ylabel(L"$z$ (m)");
		 colorbar();
	 tight_layout()

end

end # module
