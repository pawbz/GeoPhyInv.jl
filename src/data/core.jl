
Data=Array{NamedD{Recs},1}
function Array{NamedD{Recs},1}(grid::StepRangeLen, geom::Geom, fields::Vector{Symbol}) 
	return [NamedD(grid,Recs(geom[i].nr),fields) for i in 1:length(geom)]
end
