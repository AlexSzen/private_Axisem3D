from netCDF4 import Dataset 

kernel_dir = '/home/alex/Desktop/phd/private_Axisem3D/AxiSEM3D/build/output/kernels/'
num_kers = 10


totNus = 0
totElems = 0
for i in range(num_kers):
	
	ker_name = "kernels_db_%s.nc4" % i
	f = Dataset(kernel_dir+ker_name, 'r')
	totNus += len(f.variables["Kernels"][:])
	totElems += len(f.variables["Nus"][:])
	f.close()
	
f_merged = Dataset(kernel_dir + "kernels_db.nc4", "w", "NETCDF4")

f_merged.createDimension("ncdim_1", 1)
f_merged.createDimension("ncdim_%s"%totElems, totElems)
f_merged.createDimension("ncdim_%s"%totNus, totNus)
f_merged.createDimension("ncdim_5", 5)
f_merged.createDimension("ncdim_6", 6)

merged_ker = f_merged.createVariable("Kernels", "f8", ("ncdim_1","ncdim_%s"%totNus, "ncdim_6", "ncdim_5", "ncdim_5") )
merged_nus = f_merged.createVariable("Nus", "i", ("ncdim_%s"%totElems) )
merged_nrs = f_merged.createVariable("Nrs", "i", ("ncdim_%s"%totElems) )

pos_nu = 0
pos_elem = 0
for i in range(num_kers):
	
	ker_name = "kernels_db_%s.nc4" % i
	f = Dataset(kernel_dir+ker_name, 'r')
	merged_ker[0, pos_nu:pos_nu + len(f.variables["Kernels"][:]), :, :, :] = f.variables["Kernels"][:]
	merged_nus[pos_elem: pos_elem + len(f.variables["Nus"][:])] = f.variables["Nus"][:]
	merged_nrs[pos_elem: pos_elem + len(f.variables["Nrs"][:])] = f.variables["Nrs"][:]

	pos_nu += len(f.variables["Kernels"][:])
	pos_elem += len(f.variables["Nus"][:])
	f.close()
	
f_merged.close()
