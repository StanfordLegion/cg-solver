import "regent"

local c = regentlib.c

require("linalg")

task main()
	var conf : LAConfig
	conf.p = 15
	conf.p2 = 2
	conf.tol = 0.1
	--read A
	var matf = c.fopen("ignored/mat.bin", "rb")
	conf.dim = read_int(matf, 4)
	var Am = read_int(matf, 4)
	conf.m = Am
	c.printf("%dx%d A matrix, %d nonzero values\n",conf.dim ,conf.dim ,Am)
	var A = region(ispace(int1d,Am),MatEl)
	for el in A do
		el.row = read_int(matf, 4)
		el.col = read_int(matf, 4)
	end
	for el in A do
		el.val = read_double(matf)
	end
	--pm(A)
	var rowptr = region(ispace(int1d,conf.dim+1),int32)
	var j : int32 = 1
	for i = 1,conf.m,1 do
		if A[i].col < A[i-1].col then
			rowptr[j] = i
			j += 1
		end
	end
	rowptr[conf.dim] = conf.m
	var Acrs = region(ispace(int1d,conf.m),CrsEl)
	for i in Acrs.ispace do
		Acrs[i].val = A[i].val
		Acrs[i].col = A[i].col
	end
	__fence(__execution, __block)
	--read b
	var bf = c.fopen("ignored/b.bin", "rb")
	var bn = read_int(bf, 4)
	c.printf("%dx1 b\n",bn)
	var b = region(ispace(int1d,bn),double)
	for i in ispace(int1d,bn) do
		b[i] = read_double(bf)
	end
	var x = region(ispace(int1d,bn),double)
    conf.p = 15
    conf.p2 = 1
	__fence(__execution, __block)
  	var ts_start = c.legion_get_current_time_in_micros()
	__fence(__execution, __block)
	  AggroPreConCG(Acrs,rowptr,b,x,conf)
    __fence(__execution, __block)
	var ts_stop = c.legion_get_current_time_in_micros()
  	c.printf("%.4f sec\n",(ts_stop - ts_start) * 1e-6)

	var i1 = region(ispace(int1d,bn),double)
	MatVecMulInner(A,x,i1)
	__fence(__execution, __block)
	var s : double = 0
	for i in i1.ispace do
		s += (i1[i]-b[i])*(i1[i]-b[i]);
	end
	__fence(__execution, __block)
	c.printf("||Ax-b||^2 : %f\n",s)
	__fence(__execution, __block)
	c.printf("conf.p: %d\n",conf.p)
	c.printf("conf.p2: %d\n",conf.p2)
	__fence(__execution, __block)
end


regentlib.start(main)

