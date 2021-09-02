
"""
A mutable type that bundles multi-component records at multiple receivers.
```julia
data=Records(tgrid, ageom, [:p, :vx])

```
Here, we initialized records, of `:p` and `:vx` fields, in time domain for receivers and supersources in `ageom`.

## Indexing
* `data[i]` : the records due to ith supersource of all fields
* `data[i][:p]` : extract field `:p` 
* `data[i][:nr]` : number of receivers (same as in `ageom`)
* `data.grid` : returns `tgrid` 
As mutable objects in Julia are like containers that might hold different values over time, `data` can be modified.
"""
Records=Array{NamedD{Recs},1}


function Array{NamedD{Recs},1}(grid::StepRangeLen, ageom::AGeom, fields::Vector{Symbol}) 
	return [NamedD(grid,Recs(ageom[i].nr),fields) for i in 1:length(ageom)]
end
