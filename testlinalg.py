import numpy as np 
import scipy.sparse as sp 
import scipy.sparse.linalg 
import time
import os
#os.environ["MKL_NUM_THREADS"] = "1"
#os.environ["NUMEXPR_NUM_THREADS"] = "1"
#os.environ["OMP_NUM_THREADS"] = "1"

np.random.seed(2)
#n = int(50e3)
#m = int(18e6)
n = int(50e3)
m = int(90e4)
#n = int(5)
#m = int(21)
inds = np.random.randint(0,n,(m,2),dtype=np.int32)
#inds = np.vstack((inds,np.arange(n)[np.newaxis].T*np.ones((1,2))))
#inds = inds[np.argsort(inds[:,0])]

vals = np.random.randn(m)
x = np.random.randn(n)


mat = sp.coo_matrix((vals,(inds[:,0],inds[:,1])),shape=(n,n))
A = mat.T.dot(mat)
inds = np.zeros((A.count_nonzero(),2))
inds[:,0] = sp.find(A)[1]
inds[:,1] = sp.find(A)[0]
vals = sp.find(A)[2]
m = A.count_nonzero()
b = A.T.dot(x)
#print(inds)
print(inds.shape,vals.shape)
f = open("ignored/mat.bin","wb")
np.array([n,m],dtype=np.int32).tofile(f)
inds.astype(np.int32).tofile(f)
vals.tofile(f)
bf = open("ignored/b.bin","wb")
np.array([n],dtype=np.int32).tofile(bf)
b.tofile(bf)
xf = open("ignored/x.bin","wb")
np.array([n],dtype=np.int32).tofile(xf)
x.tofile(xf)

A = A.tocsr()
print(np.linalg.norm(b))
print("solving:")
start = time.time()
x = sp.linalg.cg(A,b)[0]
#r = A.dot(b)
end = time.time()
print("time: " + str(end - start) + "s")
#print(r)
print(np.sum(A.dot(x)),np.sum(b))#,callback=lambda x : print(x)))
exit()
mat = sp.coo_matrix((vals,(inds[:,0],inds[:,1])),shape=(n,n))

start = time.time()
r = mat.T.dot(mat)
end = time.time()
print("time: " + str(end - start) + "s")
print("nonzero: %f" % r.count_nonzero())
print(r[:,0])