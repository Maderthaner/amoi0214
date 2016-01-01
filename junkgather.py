import numpy as np
import psana

# parallelization
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

x = []
y = []

for i in range(1,3):
	x.append([i,rank])
	y.append([i+3,rank])

xall = comm.allgather(x)
yall = comm.allgather(y)

if rank==1:
	print 'xall:',xall
	print'yall:',yall