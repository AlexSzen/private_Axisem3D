from netCDF4 import Dataset 
import numpy as np

wvf_dir = '/home/alex/Desktop/phd/private_Axisem3D/AxiSEM3D/build/output/wavefields/'
wvf_name_merged = 'wavefield_db_fwd.nc4'
num_procs = 10


totNus = 0
totElems = 0
for i in range(num_procs):
	
	wvf_name = "wavefield_db_fwd_%s.nc4" % i
	f = Dataset(wvf_dir+wvf_name, 'r')
	totNus += np.sum(f.variables["Nus"][:])
	totElems += len(f.variables["Nus"][:])

	f.close()
	
f_merged = Dataset(wvf_dir + wvf_name_merged, "w", "NETCDF4")

f_merged.createDimension("ncdim_1", 4)
f_merged.createDimension("ncdim_%s"%totElems, totElems)
f_merged.createDimension("ncdim_%s"%totNus, totNus)
f_merged.createDimension("ncdim_5", 5)
f_merged.createDimension("ncdim_6", 6)

merged_wvf = f_merged.createVariable("displacement_wavefield", "f8", ("ncdim_1","ncdim_%s"%totNus, "ncdim_6", "ncdim_5", "ncdim_5") )
merged_nus = f_merged.createVariable("Nus", "i", ("ncdim_%s"%totElems) )
merged_nrs = f_merged.createVariable("Nrs", "i", ("ncdim_%s"%totElems) )
merged_dd = f_merged.createVariable("domain_decomposition", "i", ("ncdim_%s"%totElems) )
merged_s = f_merged.createVariable("mesh_S", "f8", ("ncdim_%s"%totElems, "ncdim_5", "ncdim_5") )
merged_z = f_merged.createVariable("mesh_Z", "f8", ("ncdim_%s"%totElems, "ncdim_5", "ncdim_5") )
merged_sem = f_merged.createVariable("sem_mesh", "i", ("ncdim_%s"%totElems, "ncdim_5", "ncdim_5") )

pos_nu = 0
pos_elem = 0
for i in range(num_procs):
	
	wvf_name = "wavefield_db_fwd_%s.nc4" % i
	f = Dataset(wvf_dir+wvf_name, 'r')
	wvf = f.variables["displacement_wavefield"][:]
	nuproc =  np.sum(f.variables["Nus"][:])
	elemproc = len(f.variables["Nus"][:])
	
	if i == 0:
		merged_dd[:] = f.variables["domain_decomposition"][:]
#	merged_wvf[:, pos_nu:pos_nu + nuproc, :, :, :] = wvf
	merged_nus[pos_elem: pos_elem + elemproc] = f.variables["Nus"][:]
	merged_nrs[pos_elem: pos_elem + elemproc] = f.variables["Nrs"][:]
	merged_s[pos_elem: pos_elem + elemproc] = f.variables["mesh_S"][:]
	merged_z[pos_elem: pos_elem + elemproc] = f.variables["mesh_Z"][:]
	merged_sem[pos_elem: pos_elem + elemproc] = f.variables["sem_mesh"][:]

	pos_nu += nuproc
	pos_elem += elemproc
	f.close()
	
f_merged.close()
