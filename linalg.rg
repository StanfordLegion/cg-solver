import "regent"

local c = regentlib.c
require("quiniela_io")


fspace MatEl {
  row : int32,
  col : int32,
  val : double,
}

fspace CrsEl {
	col : int32,
	val : double,
}

fspace LAConfig {
	p : int8,
	p2 : int8,
	dim : int32,
	m : int32,
	tol : double,
}

fspace MatPtr(r : region(ispace(int1d),CrsEl)) {
	el : ptr(CrsEl, r),
}

task NormInner(vec : region(ispace(int1d),double)) where reads(vec) do
	var r : double = 0
	for i in vec.ispace do
		r += vec[i]*vec[i]
	end
	return r
end

task VecAddScaleSelfInner(vec1 : region(ispace(int1d),double),
				 vec2 : region(ispace(int1d),double),
				 a : double)
where reads(vec1,vec2), writes(vec2) do
	__demand(__vectorize)
	for i in vec2.ispace do
		vec2[i] = vec1[i] + a*vec2[i]
	end
end


task VecScaleAddSelfInner(vec1 : region(ispace(int1d),double),
				 vec2 : region(ispace(int1d),double),
				 a : double)
where reads(vec1,vec2), writes(vec1) do
	__demand(__vectorize)
	for i in vec2.ispace do
		vec1[i] += a*vec2[i]
	end
end

task VecScaleAddInner(vec1 : region(ispace(int1d),double),
				 vec2 : region(ispace(int1d),double),
				 a : double,
				 result : region(ispace(int1d),double))
where reads(vec1,vec2), writes(result) do
	__demand(__vectorize)
	for i in result.ispace do
		result[i] = vec1[i] + a*vec2[i]
	end
end

task VecAddInner(vec1 : region(ispace(int1d),double),
				 vec2 : region(ispace(int1d),double),
				 result : region(ispace(int1d),double))
where reads(vec1,vec2), writes(result) do
	__demand(__vectorize)
	for i in result.ispace do
		result[i] = vec1[i] + vec2[i]
	end
end

task VecSubtInner(vec1 : region(ispace(int1d),double),
				 vec2 : region(ispace(int1d),double),
				 result : region(ispace(int1d),double))
where reads(vec1,vec2), writes(result) do
	__demand(__vectorize)
	for i in result.ispace do
		result[i] = vec1[i] - vec2[i]
	end
end

task DotInner(vec1 : region(ispace(int1d),double),
		      vec2 : region(ispace(int1d),double))
where reads(vec1,vec2) do
	var result : double = 0
	for i in vec1.ispace do
		result += vec1[i] * vec2[i]
	end
	return result
end

task DotInner2(vec1 : region(ispace(int1d),double),
		      vec2 : region(ispace(int1d),double))
where reads(vec1,vec2) do
	var result : double = 0
	for i in vec1.ispace do
		result += vec1[i] * vec2[i]
	end
	return result
end


task VecAdd(vec1 : region(ispace(int1d),double),
			vec2 : region(ispace(int1d),double),
			result : region(ispace(int1d),double),
			config : LAConfig)
where reads(vec1,vec2), writes(result) do
	var res_p_ind = ispace(int1d,config.p)
	var res_p = partition(equal,result,res_p_ind)
	for i in res_p_ind do
		VecAddInner(vec1,vec2,res_p[i])
	end
end


task MatVecMulInner(mat : region(ispace(int1d),MatEl),
					vec : region(ispace(int1d),double),
					result : region(ispace(int1d),double))
where reads(mat,vec), reduces +(result) do	
	for el in mat do
	    result[el.row] += el.val * vec[el.col]
	end
end

task MatVecCRSInner(mat : region(ispace(int1d),CrsEl),
					vec : region(ispace(int1d),double),
					result : region(ispace(int1d),double),
					rowst : int32)
where reads(mat,vec), reduces+(result) do	
	var cts : int32 = 0
	var col_last : int32 = 0
	for el in mat do
	    if el.col < col_last then
	    	cts+=1
	    end
	    result[rowst + cts] += el.val * vec[el.col]
	    col_last = el.col
	    --c.printf("rst: %d, count: %d, row: %d, col: %d\n",rowst,cts,el.row,el.col)
	end
end

task MatVecCRSInner2(mat : region(ispace(int1d),CrsEl),
					 vec : region(ispace(int1d),double),
					 result : region(ispace(int1d),double),
					 rowptrs : region(ispace(int1d),int32))
where reads relaxed(mat,vec,rowptrs), writes(result) do	
  --var total : int64 = 0
	for i in rowptrs.ispace do
    var res : double = 0
		for j = rowptrs[i],rowptrs[i+1],1 do
	    	res += mat[j].val * vec[mat[j].col]
    --    total += 1
	  end
    result[i] = res;
	end
  --c.printf("did %d\n", total)
end

task pv(v : region(ispace(int1d),double)) where reads(v) do
	c.printf("[ ")
	for i in v.ispace do
		c.printf("%f, ",v[i],i)
	end	
	c.printf("]\n")
end

task pm(A : region(ispace(int1d),MatEl)) where reads(A) do
	c.printf("======= MX =======\n")
	for e in A do
		c.printf("(%d, %d) : %f  \n", e.row, e.col, e.val)
	end
	c.printf("==================\n")
end
task pmc(A : region(ispace(int1d),CrsEl)) where reads(A) do
	c.printf("======= MX =======\n")
	for e in A do
		c.printf("(%d) : %f  \n", e.col, e.val)
	end
	c.printf("==================\n")
end
task countWraps(mat : region(ispace(int1d),MatEl)) where reads(mat) do
	var col_last : int32 = 0
	var cts : int32 = 0
	for el in mat do
		if el.col < col_last then 
			cts += 1
		end
		col_last = el.col
	end
	return cts
end

task MatVecMul(mat : region(ispace(int1d),MatEl),
               vec : region(ispace(int1d),double),
               result : region(ispace(int1d),double),
               config : LAConfig)
where reads(mat,vec,result), writes(result) do
	var mat_p_ind = ispace(int1d, config.p)
	var mat_p = partition(equal,mat,mat_p_ind)
	var rowst = region(mat_p_ind,int32)

	for i = 1,config.p,1 do
		rowst[i] = countWraps(mat_p[i-1])
		var p2 = mat_p[i]
		if mat[mat_p[i].ispace.bounds.lo].col < mat[mat_p[i].ispace.bounds.lo-1].col then
			rowst[i] += 1
		end 
	end
	for i = 1,config.p,1 do
		rowst[i] = rowst[i]+rowst[i-1]
		--c.printf("%d,%d\n",i,rowst[i])
	end
	var crsmat = region(ispace(int1d,config.m),CrsEl)
	for i in crsmat.ispace do
		crsmat[i].val = mat[i].val
		crsmat[i].col = mat[i].col
	end

	var rowptr = region(ispace(int1d,config.dim+1),int32)
	var j : int32 = 1
	for i = 1,config.m,1 do
		if mat[i].col < mat[i-1].col then
			rowptr[j] = i
			j += 1
		end
	end
	rowptr[config.dim] = config.m
	var rowptr_p = partition(equal,rowptr,mat_p_ind)
	var ts_start = c.legion_get_current_time_in_micros()
	__fence(__execution, __block)
	for i in mat_p_ind do
		MatVecCRSInner2(crsmat,vec,result,rowptr_p[i])
	end
	__fence(__execution, __block)
	var ts_stop = c.legion_get_current_time_in_micros()
  	c.printf("CRS2: %.4f sec\n",(ts_stop - ts_start) * 1e-6)
	__fence(__execution, __block)
	--pv(result)
	__fence(__execution, __block)
	fill(result,0)
	__fence(__execution, __block)
	var crs_p = partition(equal,crsmat,mat_p_ind)
	__fence(__execution, __block)
  	ts_start = c.legion_get_current_time_in_micros()
	for i in mat_p_ind do
		--c.printf("%d ======\n",i)
	    MatVecCRSInner(crs_p[i],vec,result,rowst[i])
	end
	__fence(__execution, __block)
	ts_stop = c.legion_get_current_time_in_micros()
  	c.printf("CRS: %.4f sec\n",(ts_stop - ts_start) * 1e-6)
  	--pv(result)
  	__fence(__execution, __block)
	fill(result,0)
	ts_start = c.legion_get_current_time_in_micros()
	__fence(__execution, __block)
	for i in mat_p_ind do
	    MatVecMulInner(mat_p[i],vec,result)
	end
	__fence(__execution, __block)
	ts_stop = c.legion_get_current_time_in_micros()
  	c.printf("COO: %.4f sec\n",(ts_stop - ts_start) * 1e-6)
  	--pv(result)
end

task vanillaCG(	A : region(ispace(int1d),CrsEl),
			   	rptr : region(ispace(int1d),int32),
			   	b : region(ispace(int1d),double),
			   	x : region(ispace(int1d),double),
				conf : LAConfig)
where reads(A,b,rptr,x), writes(x) do
	fill(x,0)
	-- init x,d,r
	var d = region(ispace(int1d,conf.dim),double)
	var r = region(ispace(int1d,conf.dim),double)
	var r2 = region(ispace(int1d,conf.dim),double)
	--ispaces and paritions
	var p_ind = ispace(int1d, conf.p)
	var A_p = partition(equal,A,p_ind)
	var rptr_p = partition(equal,rptr,p_ind)
	var b_p = partition(equal,b,p_ind)
	var x_p = partition(equal,x,p_ind)
	var d_p = partition(equal,d,p_ind)
	var r_p = partition(equal,r,p_ind)
	var r2_p = partition(equal,r2,p_ind)

	--e45
	var i1 = region(ispace(int1d,conf.dim),double)
	var i1_p = partition(equal,i1,p_ind)
	__demand(__parallel)
	for i in p_ind do
		MatVecCRSInner2(A,x,i1_p[i],rptr_p[i])
	end
	for i in p_ind do
		VecSubtInner(b_p[i],i1_p[i],d_p[i])
	end 
	copy(d,r)
	var i3 : double
	i3 = 0
	for i in p_ind do
		i3 += DotInner(r_p[i],r_p[i])
	end
	var a : double
	var beta : double
	var i2 : double
	var i4 : double
	var rn : double
	var tolsq : double = conf.tol*conf.tol
	for i=0,1000,1 do 
		--e46
		fill(i1,0)
		__demand(__parallel)
		for i in p_ind do
			MatVecCRSInner2(A,d,i1_p[i],rptr_p[i])
		end
		i2 = 0
		__demand(__parallel)
		for i in p_ind do
			i2 += DotInner(d_p[i],i1_p[i])
		end
		a = i3/i2
		--e46.5
		__demand(__parallel)
		for i in p_ind do
			VecScaleAddSelfInner(x_p[i],d_p[i],a)
		end
		--e47
		__demand(__parallel)		
		for i in p_ind do
			VecScaleAddInner(r_p[i],i1_p[i],(-a),r2_p[i])
		end 
		--e48
		beta = 0
		i4 = 0
		__demand(__parallel)
		for i in p_ind do
			i4 += DotInner(r2_p[i],r2_p[i])
		end
		beta = i4/i3
		--e49
		__demand(__parallel)
		for i in p_ind do
			VecAddScaleSelfInner(r2_p[i],d_p[i],beta)
		end
		copy(r2,r)
		rn = 0
		__demand(__parallel)
		for i in p_ind do
			rn += NormInner(r)
		end
		c.printf("error: %f\n",rn)
		__fence(__execution, __block)
		i3 = i4
		--__fence(__execution, __block)
		--pv(r)
		--__fence(__execution, __block)
		--pv(d)
		--__fence(__execution, __block)
	end
end

task PreConCG(	A : region(ispace(int1d),CrsEl),
			   	rptr : region(ispace(int1d),int32),
			   	b : region(ispace(int1d),double),
			   	x : region(ispace(int1d),double),
				  conf : LAConfig)
where reads(A,b,rptr,x), writes(x) do
  __fence(__execution, __block) -- This blocks to make sure we only time the pagerank computation
    var ts_start = c.legion_get_current_time_in_micros()
  __fence(__execution, __block)
	var cts : int32 = 0
	var col_last : int32 = 0
	var MI = region(ispace(int1d,conf.dim),CrsEl)
	var MIrptr = region(ispace(int1d,conf.dim+1),int32)
	for i = 0,conf.dim+1,1 do
		MIrptr[i] = i
	end
	for i in rptr.ispace do
		for j = rptr[i],rptr[i+1],1 do
	    	if A[j].col == MIrptr[i] then 
	    		MI[i].val = 1./A[j].val
	    		MI[i].col = i
	    		--c.printf("%d,%d,%f\n",MIrptr[i],A[j].col,MI[i].val)
	    	end
	    end 
	end
	fill(x,0)
	-- init x,d,r
	var d = region(ispace(int1d,conf.dim),double)
	var r = region(ispace(int1d,conf.dim),double)
	var r2 = region(ispace(int1d,conf.dim),double)
	-- ispaces and paritions
	var p_ind = ispace(int1d, conf.p)
	var A_p = partition(equal,A,p_ind)
	var rptr_p = partition(equal,rptr,p_ind)
	var MI_p = partition(equal,MI,p_ind)
	var MIrptr_p = partition(equal,MIrptr,p_ind)

	--var mxptrs = region(ispace(int1d,conf.m),MatPtr(A))
	--var mxptrs_p = partition(equal,mxptrs,p_ind)
	--unsafe_cast(ptr(MatEl, A), A[0])

	--e45
	var i1 = region(ispace(int1d,conf.dim),double)
	var i5 = region(ispace(int1d,conf.dim),double)
	var i1_p = partition(equal,i1,p_ind)
	var i5_p = partition(equal,i5,p_ind)
	__demand(__parallel)
	for i in p_ind do
		MatVecCRSInner2(A,x,i1_p[i],rptr_p[i])
	end

  VecSubtInner(b,i1,r)
	MatVecCRSInner2(MI,r,i5,MIrptr)

	copy(i5,d)
	var i3 : double
	i3 = DotInner(r, i5)

	var a : double
	var beta : double
	var i2 : double
	var i4 : double
	var rn : double
	var tolsq : double = conf.tol*conf.tol

  __fence(__execution, __block) -- This blocks to make sure we only time the pagerank computation
    var ts_end = c.legion_get_current_time_in_micros()
    c.printf("setup time %f ms\n", (ts_end-ts_start)*1e-3);
  __fence(__execution, __block)

	for it=0,1000,1 do 

		__demand(__parallel)
		for i in p_ind do
			MatVecCRSInner2(A,d,i1_p[i],rptr_p[i])
		end

		i2 = DotInner(d, i1);
		a = i3/i2
		VecScaleAddSelfInner(x,d,a)
		VecScaleAddSelfInner(r,i1,(-a))
		MatVecCRSInner2(MI,r,i5,MIrptr)
		beta = DotInner2(r,i5)/i3
		VecAddScaleSelfInner(i5,d,beta)
		i3 *= beta
    if it % 5 == 0 then
      rn = NormInner(r)
      if rn < tolsq then
        c.printf("iteration %d\n", it);
        break
      end
    end
	end
end

task AggroPreConCG(	A : region(ispace(int1d),CrsEl),
			   	rptr : region(ispace(int1d),int32),
			   	b : region(ispace(int1d),double),
			   	x : region(ispace(int1d),double),
				  conf : LAConfig)
where reads(A,b,rptr,x), writes(x) do
  --[[__fence(__execution, __block) -- This blocks to make sure we only time the pagerank computation
    var ts_start = c.legion_get_current_time_in_micros()
  __fence(__execution, __block)
  --]]
	var MI = region(ispace(int1d,conf.dim),CrsEl)
	var MIrptr = region(ispace(int1d,conf.dim+1),int32)

	var d = region(ispace(int1d,conf.dim),double)
	var r = region(ispace(int1d,conf.dim),double)
	var r2 = region(ispace(int1d,conf.dim),double)
	var p_ind = ispace(int1d, conf.p)
	var A_p = partition(equal,A,p_ind)
	var rptr_p = partition(equal,rptr,p_ind)
	var MI_p = partition(equal,MI,p_ind)
	var MIrptr_p = partition(equal,MIrptr,p_ind)

	var i1 = region(ispace(int1d,conf.dim),double)
	var i5 = region(ispace(int1d,conf.dim),double)
	var i1_p = partition(equal,i1,p_ind)
	var i5_p = partition(equal,i5,p_ind)

	var i3 : double
	var a : double
	var beta : double
	var i2 : double
	var i4 : double
	var rn : double
	var tolsq : double = conf.tol*conf.tol

	for i = 0,conf.dim+1,1 do
		MIrptr[i] = i
	end
	for i in rptr.ispace do
		for j = rptr[i],rptr[i+1],1 do
	    	if A[j].col == MIrptr[i] then 
	    		MI[i].val = 1./A[j].val
	    		MI[i].col = i
	    	end
	    end 
	end
	fill(x,0)

  --[[__fence(__execution, __block)
  var ts_end = c.legion_get_current_time_in_micros()
  c.printf("setup time %f ms\n", (ts_end-ts_start)*1e-3);
  __fence(__execution, __block)
  --]]

	__demand(__parallel)
	for i in p_ind do
		MatVecCRSInner2(A_p[i],x,i1_p[i],rptr_p[i])
	end
	VecSubtInner(b,i1,r)
	for i in p_ind do
		MatVecCRSInner2(MI,r,i5_p[i],MIrptr_p[i])
	end
  --__fence(__execution, __block)
  --c.exit(0)
	copy(i5,d)
	i3 = DotInner(r,i5)

	for it=0,1000,1 do 

		__demand(__parallel)
		for i in p_ind do
			MatVecCRSInner2(A,d,i1_p[i],rptr_p[i])
		end

		i2 = DotInner(d, i1);
		a = i3/i2
		VecScaleAddSelfInner(x,d,a)
		VecScaleAddSelfInner(r,i1,(-a))
		MatVecCRSInner2(MI,r,i5,MIrptr)
		beta = DotInner2(r,i5)/i3
		VecAddScaleSelfInner(i5,d,beta)
		i3 *= beta
    if it % 15 == 0 then
      rn = NormInner(r)
      if rn < tolsq then
        --c.printf("iteration %d\n", it);
        break
      end
    end
	end
end


task AggroNPPreConCG(	A : region(ispace(int1d),CrsEl),
			   	rptr : region(ispace(int1d),int32),
			   	b : region(ispace(int1d),double),
			   	x : region(ispace(int1d),double),
				  conf : LAConfig)
where reads(A,b,rptr,x), writes(x) do

	var d = region(ispace(int1d,conf.dim),double)
	var r = region(ispace(int1d,conf.dim),double)
	var r2 = region(ispace(int1d,conf.dim),double)
	var p_ind = ispace(int1d, conf.p)
	var A_p = partition(equal,A,p_ind)
	var rptr_p = partition(equal,rptr,p_ind)

	var i1 = region(ispace(int1d,conf.dim),double)
	var i5 = region(ispace(int1d,conf.dim),double)
	var i1_p = partition(equal,i1,p_ind)
	var i5_p = partition(equal,i5,p_ind)

	var i3 : double
	var a : double
	var beta : double
	var i2 : double
	var i4 : double
	var rn : double
	var tolsq : double = conf.tol*conf.tol

	fill(x,0)

  --[[__fence(__execution, __block)
  var ts_end = c.legion_get_current_time_in_micros()
  c.printf("setup time %f ms\n", (ts_end-ts_start)*1e-3);
  __fence(__execution, __block)
  --]]

	__demand(__parallel)
	for i in p_ind do
		MatVecCRSInner2(A,x,i1_p[i],rptr_p[i])
	end
	VecSubtInner(b,i1,r)
	copy(r,d)
	i3 = NormInner(r)

	for it=0,1000,1 do 

		__demand(__parallel)
		for i in p_ind do
			MatVecCRSInner2(A,d,i1_p[i],rptr_p[i])
		end

		i2 = DotInner(d, i1);
		a = i3/i2
		VecScaleAddSelfInner(x,d,a)
		VecScaleAddSelfInner(r,i1,(-a))
		beta = NormInner(r)/i3
		VecAddScaleSelfInner(r,d,beta)
		i3 *= beta
    if it % 2 == 0 then
      rn = NormInner(r)
      if rn < tolsq then
        --c.printf("iteration %d\n", it);
        break
      end
    end
	end
end


--[[
PreConCG:printpretty()
DotInner:printpretty()
VecScaleAddSelfInner:printpretty()
--]]
