module Plots

using StatsBase
using RecipesBase
using FFTW
using DSP
#using ImageFiltering
import GeoPhyInv.Data
import GeoPhyInv: Medium, copyto!


@userplot Seismic
"""
Plot the velocity and density seismic models.

# Arguments

* `model::Medium` : model that should be plotted

# Keyword Arguments

* `xlim::Vector{Float64}=[model.mgrid[2][1],model.mgrid[2][end]]` : minimum and maximum limits of the second dimension while plotting
* `zlim::Vector{Float64}=[model.mgrid[1][1],model.mgrid[1][end]]` : minimum and maximum limits of the first dimension while plotting
* `fields::Vector{Symbol}=[:vp, :ρ]` : fields that are to be plotted
* `contrast_flag=false` : plot only the edges of the model
* `use_bounds=false` : adjust `clim` to the bounds in the seismic model

"""
@recipe function fseismic(p::Seismic;
		   fields=[:vp, :ρ], 
		   contrast_flag=false,
		   use_bounds=false,
		  ) 
	if(contrast_flag)
		warn("ImageFiltering bug needs to be fixed")
	end
	model=p.args[1]
	nz,nx=length.(model.mgrid)

	nrow = (nx <= nz) ? 1 : length(fields)
	ncol = (nx <= nz) ? length(fields) : 1


	layout --> (nrow,ncol)
	for (i,iff) in enumerate(fields)
		name=string(iff)

		m = model[iff]
		mx = model.mgrid[2]
		mz = model.mgrid[1]
		if(contrast_flag)
			mmin=minimum(m)
			mmax=maximum(m)
			#m=imfilter(m, Kernel.Laplacian())
			mmmin=minimum(m)
			mmmax=maximum(m)
			for j in eachindex(m)
				m[j]=mmin+((m[j]-mmmin)*inv(mmmax-mmmin))*(mmax-mmin)
			end
		end
		mmin = use_bounds ? model.bounds[iff][1] : minimum(m)
		mmax = use_bounds ? model.bounds[iff][2] : maximum(m)
		@series begin        
			subplot := i
			aspect_ratio := :equal
			seriestype := :heatmap
			legend --> true
			xlabel --> "x [m]"
			ylabel --> "z [m]"
			color --> :grays
			#xlim --> (mx[1], mx[end])
			#zlim --> (mz[1], mz[end])
			title --> name
			clim --> (mmin, mmax)
			yflip := true
			mx, mz, m
		end
	end

end
#

#"""
#save current fig using matlab2tikz
#"""
#function printfig(fig; filename::String = "FIG")
#	if filename == "FIG"
#		temp = 1;
#		file = join([filename string(temp) ".tex"])
#		while isfile(file)
#			temp += 1;
#			file = join([filename string(temp) ".tex"])
#		end
#		tikzfile = join([file[1:end-4] ".tikz"]);
#	else
#		file = filename;
#		tikzfile = join([file[1:end-4] ".tikz"]);
#	end
#	#matplotlib2tikz.save(tikzfile,fig);
#	fin = open(tikzfile,"r");
#	fout = open(file,"w");
#	lines = readlines(fin);
#	close(fin);
#	rm(tikzfile);
#	lines_old = ["\\documentclass[tikz]{standalone}\n";
#	"\\usepackage{pgfplots}\n";
#	"\\pgfplotsset{compat=newest}\n";
#	"\\usepackage{amsmath}\n";
#	"\\usepackage{siunitx}\n";
#	"\\begin{document}\n";lines;"\n\\end{document}\n"];
#
#	lines = ["% Generated from Julia\n";"\\begin{center}\n";lines;"\\end{center}"]
#	write(fout,lines)
#	close(fout)
#end

#=
@userplot Geom

"""
Plot acquisition geometry `Geom` on
and model grid.

`attrib::Symbol=:unique` : default; plots unique source and receiver positions
`ssvec::Vector{Int64}` : plot source and receivers of only these supersources
"""
@recipe function fgeom(p::Geom; ssvec=nothing, fields=[:s, :r])
	geom=p.args[1]
	if(ssvec===nothing)
		if(:r ∈ fields)
			urpos = Geom_get([geom],:urpos)
			rx=urpos[2]
			rz=urpos[1]
		else
			b=nothing
		end
		if(:s ∈ fields)
			uspos = Geom_get([geom],:uspos)
			sx=uspos[2]
			sz=uspos[1]
		else
			a=nothing
		end
	else
		if(:r ∈ fields)
			rxpos = [geom.rx[iss] for iss in ssvec]
			rzpos = [geom.rz[iss] for iss in ssvec]
			rx=vcat(rxpos...); rz=vcat(rzpos...)
		else
			b=nothing
		end
		if(:s ∈ fields)
			sxpos = [geom.sx[iss] for iss in ssvec]
			szpos = [geom.sz[iss] for iss in ssvec]
			sx=vcat(sxpos...); sz=vcat(szpos...)
		else
			a=nothing
		end
	end

	if(:r ∈ fields)
		@series begin        
			subplot --> 1
			markersize --> 7
			legend --> false
			seriestype := :scatter
			markercolor := :blue
			markershape := :utriangle
			rx, rz
		end
	end
	if(:s ∈ fields)
		@series begin        
			subplot --> 1
			legend --> false
			markersize --> 7
			markercolor := :red
			markershape := :xcross
			seriestype := :scatter
			sx, sz
		end
	end


end
=#
	

@userplot Spectrum


@recipe function pspectrum(p::Spectrum; tgrid=nothing)
	wav=p.args[1]
	if(tgrid===nothing)
		x=0:length(wav)
		tgrid=range(0.0, stop=Float64(length(wav)-1), step=1.0)
	end

	fgrid= DSP.rfftfreq(length(tgrid), inv(step(tgrid)))
	powwav = (abs.(FFTW.rfft(wav, [1])).^2)
	powwavdb = 10. * log10.(powwav./maximum(powwav)) # power in decibel after normalizing

	@series begin        
		subplot := 1
		legend := false
		ylabel := "power (dB)"
		xlabel := "frequency (Hz)"
		fgrid, powwavdb
	end
end

@userplot Src

"""
Plot the source wavelet used for acquisition.

# Arguments

* `srcwav::Src` : source acquisition parameters
"""
@recipe function psrc(p::Src)
	srcwav=p.args[1]
	tgrid=srcwav.tgrid

	wav=srcwav.wav[1,1]
	layout := (2,1)
	@series begin        
		subplot := 1
		legend := false
		xlabel := "time [s]"
		ylabel := "amplitude" 
		tgrid, wav 
	end

	fgrid= DSP.rfftfreq(length(tgrid), inv(step(tgrid)))
	powwav = (abs.(FFTW.rfft(wav, [1])).^2)
	powwavdb = 10. * log10.(powwav./maximum(powwav)) # power in decibel after normalizing
	@series begin        
		subplot := 2
		ylabel := "power [dB]"
		xlabel := "frequency [Hz]"
		legend := false
		fgrid, powwavdb
	end

end

@userplot TD
"""
Plot time-domain data of type `Data.TD`

# Arguments
* `td::Vector{Data.TD}` : time-domain data to be compared

# Keyword Arguments
* `ssvec::Vector{Vector{Int64}}=fill([1], length(td))` : supersource vector to be plotted
* `fields::Vector{Int64}=[1]` : field vector to be plotted
* `tr_flag::Bool=false` : plot time-reversed data when true
* `attrib::Symbol=:wav` : specify type of plot
"""
@recipe function ptd(p::TD;
	    fields=[:P],
            tr_flag=false, 
            ssvec=[1],
	    wclip_perc=0.0,
	    bclip_perc=0.0,
	    )

	    #wclip::Vector{Float64}=[maximum(broadcast(maximum, td[id].d)) for id in 1:length(td)],
	    #bclip::Vector{Float64}=[minimum(broadcast(minimum, td[id].d)) for id in 1:length(td)],
	   
	dat=p.args[1]
	any(ssvec .> dat.geom.nss) && error("invalid ssvec")
	ns = length(ssvec);
	nr = maximum(dat.geom.nr);
	dd=getfield(dat,fieldnames(typeof(dat))[1])
	fieldvec = findall(in(fields),dat.fields)
	if(tr_flag)
		dp = hcat(dd[ssvec,fieldvec][end:-1:1,:]...);
		dz = dat.tgrid
		dx = 1:size(dp,2)
	else
		dp = hcat(dd[ssvec,fieldvec][:,:]...);
		dz = dat.tgrid[end:-1:1]
		dx = 1:size(dp,2)
	end

	dmin=minimum(dp)
	dmax=maximum(dp)
	offset=abs(dmax-dmin)
	dmin=dmin+bclip_perc*inv(100)*offset
	dmax=dmax-wclip_perc*inv(100)*offset
	@series begin        
		subplot := 1
		aspect_ratio --> length(dz)/length(dx)
		seriestype := :heatmap
		legend --> true
		xlabel --> "receiver index"
		ylabel --> "time [s]"
		color --> :grays
		clim --> (dmin, dmax)
		yflip := true
		dx, dz, dp
	end
end



#=
function mscatter(x, y, titname="", axla=["",""])
    fact=1 
    x=x[1:fact:end]
    x[:] /= norm(x) # normalize inputs before scatter
    y=y[1:fact:end]
    y[:] /= norm(y) # normalize inputs before scatter 
    scatter(x,y,marker="x",s=0.5,color="gray")
    ax=gca()
    ax[:set_aspect]("equal")
    ax[:spines]["top"][:set_visible](false)
    ax[:spines]["right"][:set_visible](false)
    ax[:spines]["left"][:set_position]("center")
    ax[:spines]["bottom"][:set_position]("center")
    setp(ax[:get_yticklabels](false),color="none")
    setp(ax[:get_xticklabels](false),color="none")
    xmax = maximum(abs, x)
    xmax = iszero(xmax) ? 1.0 : xmax
    ymax = maximum(abs, y)
    ymax = iszero(ymax) ? 1.0 : ymax
    li = 1.2*max(xmax, ymax)
    #xlim(-li, li)
    #ylim(-li, li)
    xlabel(axla[1])
    ylabel(axla[2])
    title(titname)
    return nothing
end

function myhist(s, dist)
    fact=1
    x=s[1,1:fact:end]
    y=s[2,1:fact:end]
    nbins=div(length(x), 10)
    h1 = normalize(fit(Histogram, x, nbins=200, closed=:right))
    h2 = normalize(fit(Histogram, y, nbins=200, closed=:right))

#    h = PyPlot.plt[:hist](x, 200, color="red", histtype="barstacked", normed=true)
 #   h = PyPlot.plt[:hist](y, 200, color="blue",histtype="barstacked", normed=true)
    plot(h1.edges[1][2:end], h1.weights,color="red", linewidth=0.7)
    plot(h2.edges[1][2:end], h2.weights,color="blue", linewidth=0.7)
    ax=gca()
    xmax = 1.2*maximum(abs, x)
    xx=linspace(-xmax, xmax, 300)
    xlim(-xmax, xmax)
    yymax=0.
    for idist in 1:length(dist)
        yy = pdf(dist[idist], xx)
        yymax=maximum([maximum(yy), yymax])
        ax[:plot](xx, yy, color="black",linestyle="--", linewidth=1)
    end
    (yymax ≠ 0.) &&    ylim([0, 1.2*yymax])
    setp(ax[:get_yticklabels](),color="none")
    setp(ax[:get_xticklabels](),color="none")
    #setp(ax[:tick_params](false), labelbottom="off")
    title("Marginal Densities")
    return nothing
end
=#

end # module
