
"""
A mutable type that bundles multi-component records at receiver array.
```julia
records=Records(tgrid, ageom, [:p, :vx])

```
Here, we initialized records, of `:p` and `:vx` fields, in time domain for receivers and supersources in `ageom`.

## Indexing
* `records[i]` : the records due to ith supersource of all fields
* `records[i][:p]` : extract field `:p` 
* `records[i][:nr]` : number of receivers (same as in `ageom`)
* `records.grid` : returns `tgrid` 
As mutable objects in Julia are like containers that might hold different values over time, `records` can be modified.
"""
Records = Array{NamedD{Recs},1}


function Array{NamedD{Recs},1}(grid::StepRangeLen, ageom::AGeom, fields::Vector{Symbol})
    return [NamedD(grid, Recs(ageom[i].nr), fields) for i = 1:length(ageom)]
end
