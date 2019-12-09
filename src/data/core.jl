
Data=Array{NamedD{Recs},1}
function Array{NamedD{Recs},1}(grid::StepRangeLen, ageom::AGeom, fields::Vector{Symbol}) 
	return [NamedD(grid,Recs(ageom[i].nr),fields) for i in 1:length(ageom)]
end
