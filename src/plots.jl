@recipe function plot(medium::Medium, ageom=nothing; fields=MediumParameters(medium))
    mx = medium.grid[2]
    mz = medium.grid[1]
    as = length(mz) / length(mx)
    layout := ((as < 1) ? (length(fields), 1) : (1, length(fields)))
    margin --> 5Measures.mm
    framestyle := [:grid :grid :grid]
    if (as < 0.5)
        as1 = 0.5
        lp = :outertopleft
    elseif (as > 2)
        as1 = 2
        lp = :outertop
    else
        as1 = 1
        lp = :outertopleft
    end
    aspect_ratio --> (inv(as) * as1)
    size --> ((as < 1) ? (800, length(fields) * 300) : (400 * length(fields), 400))
    yflip := true
    legend --> lp
    xguide --> "x"
    yguide --> "z"
    xlims --> (mx[1], mx[end])
    zlims --> (mz[1], mz[end])
    for (j, field) in enumerate(fields)
        units = field in [:vp, :vs] ? " [m/s]" : " [kg/m3]"
        title1 = string(field, units) # add using here later
        @series begin
            seriestype := :heatmap
            subplot := j
            # seriescolor --> colorschemes[:roma]
            title --> title1
            clims --> Tuple(getfield(medium, field).bounds)
            seriescolor --> :grays
            mx, mz, getfield(medium, field).m
        end
        if (!(ageom === nothing))
            @assert ndims(medium) == length(dim_names(ageom))
            @series begin
                subplot := j
                title --> title1
                label --> "R"
                seriestype := :scatter
                clims --> Tuple(getfield(medium, field).bounds)
                ageom, Recs()
            end
            @series begin
                subplot := j
                title --> title1
                label --> "S"
                seriestype := :scatter
                clims --> Tuple(getfield(medium, field).bounds)
                ageom, SSrcs()
            end
        end
    end
end

@recipe function plot(ageom::Union{AGeom,AGeomss}, ::SSrcs)
    seriestype := :scatter
    legend --> false
    markersize --> 7
    markercolor := :red
    markershape := :xcross
    labels = (length(dim_names(ageom)) == 3) ? [:x, :y, :z] : [:x, :z, nothing]
    xlabel --> labels[1]
    ylabel --> labels[2]
    zlabel --> labels[3]
    dnames = (length(dim_names(ageom)) == 3) ? [:x, :y, :z] : [:x, :z]
    isa(ageom, AGeom) ? tuple([vcat([ag.s[d] for ag in ageom]...) for d in dnames]...) :
    tuple([ageom.s[d] for d in dnames]...)
end

@recipe function plot(ageom::Union{AGeom,AGeomss}, ::Recs)
    seriestype := :scatter
    markersize --> 7
    legend --> false
    markercolor := :yellow
    markershape := :utriangle
    labels = (length(dim_names(ageom)) == 3) ? [:x, :y, :z] : [:x, :z, nothing]
    xlabel --> labels[1]
    ylabel --> labels[2]
    zlabel --> labels[3]
    dnames = (length(dim_names(ageom)) == 3) ? [:x, :y, :z] : [:x, :z]
    isa(ageom, AGeom) ? tuple([vcat([ag.r[d] for ag in ageom]...) for d in dnames]...) :
    tuple([ageom.r[d] for d in dnames]...)
end

@recipe function plot(dat::NamedD{Recs}, clip_perc=0.0)
    layout := (1, length(dat.d))
    size --> (length(dat.d) * 300, 300)
    margin --> 5Measures.mm
    for (j, field) in enumerate(names(dat.d)[1])
        dp = dat[field]
        dmax = maximum(abs.(dp))
        if (!iszero(dmax))
            @series begin
                dz = dat.grid
                dx = 1:size(dp, 2)
                dmax_clip = dmax - clip_perc * inv(100) * abs(dmax)
                seriestype := :heatmap
                framestyle := :grid
                title --> field
                legend --> true
                xguide --> "receiver channel"
                yguide --> "time"
                colorbar --> false
                seriescolor --> :seismic
                clims --> (-dmax_clip, dmax_clip)
                yflip := true
                dx, dz, dp
            end
        else
            @info "iszero(data)==true; cannot plot"
        end
    end
end


@recipe function plot(dat::NamedD{Srcs})
    fgrid = FFTW.rfftfreq(length(dat.grid), inv(step(dat.grid)))
    layout := (1, 2)
    size --> (800, length(dat.d) * 200)
    margin --> 5Measures.mm
    framestyle := :grid
    for (j, field) in enumerate(names(dat.d)[1])
        dp = dat[field]
        dmax = maximum(abs.(dp))
        if (!iszero(dmax))
            D = (abs.(FFTW.rfft(dp, [1])))
            D ./= maximum(D) # normalize spectrum
            freqmax = maximum(findall(x -> x > 1e-6, D))[1]
            @series begin
                subplot := 1
                seriestype := :line
                w --> 2
                legend --> true
                label --> string(field)
                xguide --> "time"
                yguide --> "amplitude"
                dat.grid, dp
            end
            @series begin
                subplot := 2
                w --> 2
                legend --> true
                label --> string(field)
                seriestype := :line
                xguide --> "frequency (Hz)"
                yaxis := :log
                ylims --> (1e-6, 1)
                xlims --> (0, freqmax)
                fgrid, D
            end
        else
            @info "iszero(data)==true; cannot plot"
        end
    end
end


#=
Plot time-domain data of type `Records.TD`

# Arguments
* `td::Vector{Records.TD}` : time-domain data to be compared

# Keyword Arguments
* `ssvec::Vector{Vector{Int64}}=fill([1], length(td))` : supersource vector to be plotted
* `fields::Vector{Int64}=[1]` : field vector to be plotted
* `tr_flag::Bool=false` : plot time-reversed data when true
* `attrib::Symbol=:wav` : specify type of plot
"""

	ff(ssvec===nothing)
		if(:r ∈ fields)
			urpos = AGeom_get([ageom],:urpos)
			rx=urpos[2]
			rz=urpos[1]
		else
			b=nothing
		end
		if(:s ∈ fields)
			uspos = AGeom_get([ageom],:uspos)
			sx=uspos[2]
			sz=uspos[1]
		else
			a=nothing
		end
	else
		if(:r ∈ fields)
			rxpos = [ageom.r[:x][iss] for iss in ssvec]
			rzpos = [ageom.r[:z][iss] for iss in ssvec]
			rx=vcat(rxpos...); rz=vcat(rzpos...)
		else
			b=nothing
		end
		if(:s ∈ fields)
			sxpos = [ageom.s[:x][iss] for iss in ssvec]
			szpos = [ageom.s[:z][iss] for iss in ssvec]
			sx=vcat(sxpos...); sz=vcat(szpos...)
		else
			a=nothing
		end
	end

	if(:r ∈ fields)
		@series begin        
		end
	end
	if(:s ∈ fields)
		@series begin        
		end
	end


end
=#


#=
@userplot Seismic
"""
Plot the velocity and density seismic models.

# Arguments

* `model::Medium` : model that should be plotted

# Keyword Arguments

* `xlim::Vector{Float64}=[model.grid[2][1],model.grid[2][end]]` : minimum and maximum limits of the second dimension while plotting
* `zlim::Vector{Float64}=[model.grid[1][1],model.grid[1][end]]` : minimum and maximum limits of the first dimension while plotting
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
	nz,nx=length.(model.grid)

	nrow = (nx <= nz) ? 1 : length(fields)
	ncol = (nx <= nz) ? length(fields) : 1


	layout --> (nrow,ncol)
	for (i,iff) in enumerate(fields)
		name=string(iff)

		m = model[iff]
		mx = model.grid[2]
		mz = model.grid[1]
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
			xguide --> "x [m]"
			yguide --> "z [m]"
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
		xguide := "time [s]"
		yguide := "amplitude" 
		tgrid, wav 
	end

	fgrid= DSP.rfftfreq(length(tgrid), inv(step(tgrid)))
	powwav = (abs.(FFTW.rfft(wav, [1])).^2)
	powwavdb = 10. * log10.(powwav./maximum(powwav)) # power in decibel after normalizing
	@series begin        
		subplot := 2
		yguide := "power [dB]"
		xguide := "frequency [Hz]"
		legend := false
		fgrid, powwavdb
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
    xguide(axla[1])
    yguide(axla[2])
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

=#
