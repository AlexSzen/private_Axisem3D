import numpy as np 
from mpi4py import MPI 

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


counts = [10, 15, 8, 12]
displs = [0, 10, 25, 33]

if rank == 0:
	data = np.ones(45)
	print('About to scatter that baby', data)
else:
	data=None

x = np.zeros(counts[rank])
data = comm.Scatterv([data, counts, displs, MPI.DOUBLE], x, root=0)

print('rank ',rank, ' has data ', x)



