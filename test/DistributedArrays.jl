addprocs(3)
@everywhere using DistributedArrays

@everywhere function filll!(dat,p)
	svec = localindexes(dat)[end]
	println(svec)
	for (is, s) in enumerate(svec)
		a = localpart(dat)
		println([p,s])
		a[is][:].= [p,s]
		println(localpart(dat)[is][:])
	end
end
dat=DArray((7,), workers(), [nworkers()]) do I
	[zeros(2) for i in 1:nworkers()]
end

@sync begin 
	for (ip, p) in enumerate(workers())
		@async remotecall_wait(filll!, p ,dat, p)
	end
end
println(dat)

