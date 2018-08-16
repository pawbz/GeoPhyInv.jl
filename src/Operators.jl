module Operators




"""
* `n` : input vector length

# Keyword Arguments

* `np` : pad length +ve
* `nn` : pad length -ve
* `ntimes` : number of times `n+np+nn` block is repeated
"""
function pad(T::Type, n; np=0, nn=0, ntimes=1 )
	nout=(n+np+nn)*ntimes; # total length
	neach=n+np+nn # each block length 
	nout=neach*ntimes # total length
	A=spzeros(T, nout, n)
	# need first nin elements
	for it in 1:ntimes
		for i in nn+1:n+nn
			A[i+(it-1)*neach,i-nn]=one(T)
		end
	end
	return A
end


function splitpad(T::Type, n; nfrac=1, np=0, nn=0)
	(mod(n, nfrac) ≠ 0) && error("invalid nfrac")
	nb=div(n, nfrac)
	neach=nb+np+nn
	nout=neach*nfrac
	A=spzeros(T, nout, n)
	for ii in 1:nfrac
		iii=(ii-1)*neach
		for i in 1:nb
			A[i+nn+iii,i+(ii-1)*nb]=one(T)
		end
	end
	return A
end

function splittrun(T::Type, n; np=0, nn=0, nfrac=1)
	return transpose(splitpad(T,n,np=np, nn=nn, nfrac=nfrac))
end

function blkdiag(T::Type, Ain, ntimes)
	n,m=size(Ain)
	A=spzeros(T, n*ntimes, m*ntimes)

	for i in 1:ntimes
		Av=view(A,1+n*(i-1):n*i,1+m*(i-1):m*i)
		for ii in eachindex(Ain)
			Av[ii]=Ain[ii]
		end
	end
	return A
end

end
