module Plots

using PyCall, PyPlot
import PyCall: @pyimport
import PyPlot: pygui
@pyimport matplotlib2tikz
import SIT.Acquisition
import SIT.Grid
import SIT.Data

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
"""
function Geom(geom::Acquisition.Geom,
	      mgrid::Grid.M2D,
	      attrib::Symbol
	     )

plot(geom.rx, geom.rz, "v", color="blue",ms=10)
plot(geom.sx, geom.sz, "*", color="red",ms=15)


end


"""
Plot acqsrc
"""
function Src(acqsrc::Acquisition.Src, attrib::Symbol)
	
	plot(acqsrc.tgrid.x, acqsrc.wav[:,1,1])

end




"""
Plot time-domain data of type `Data.TD`
"""
function TD(td::Data.TD, attrib::Symbol)
	
	if(attrib == :seis)
		imshow(td.d[:,:,1],cmap="gray",aspect="auto")
	else
		plot(td.tgrid.x,td.d[:,:,1])
	end


end



end

