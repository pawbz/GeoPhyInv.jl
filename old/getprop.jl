
"""
## Indexing
`pa` is an instance of `SeisInvExpt`. 
* `pa.mediumm` : inverted medium on the modeling grid (after `update!`)
* `pa.mediumi` : inverted medium on the inversion grid (after `update!`)
* `pa.ageom` : acquisition 
* `pa.srcwav`: `srcwav` used to create `pa`
* `pa[:dobs]` : observed data
* `pa[:dcal]` : modeled synthetic data
"""
function Base.getindex(pa::PFWI, s::Symbol)
	@assert s in [:dobs, :dcal]
	if(s==:dobs)
		return pa.paTD.y
	elseif(s==:dcal)
		return pa.paTD.x
	end

	#=
	if(s==:snaps)
		nzd,nxd=length.(pa.c.model.grid)
		snaps_all=zeros(nzd,nxd,length(pa.c.itsnaps))
		ip=findall(in.(iss,pa.sschunks))[1]
		issp=findall(x->x==iss,pa.sschunks[ip])[]
		@sync begin
			p=procs(pa.p)[ip]
			@sync remotecall_wait(p) do 
				pap=localpart(pa.p)
				for i in eachindex(snaps_all)
					snaps_all[i]=pap.ss[issp].snaps[i]
				end
			end
		end
		return snaps_all
	elseif(s==:data)
		(iss==0) ?  (return pa.c.data[1]) :  (return pa.c.data[1][iss])
	elseif(s==:ageom)
		return pa.c.ageom
	elseif(s==:srcwav)
		return pa.c.srcwav
	end
	=#
end

