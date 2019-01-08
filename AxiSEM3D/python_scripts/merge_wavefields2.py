from netCDF4 import Dataset 
import numpy as np

wvf_dir = '/home/alex/Desktop/phd/private_Axisem3D/AxiSEM3D/build/output/instaseis_axisem3D/'
wvf_name_merged = 'wavefield_db_fwd.nc4'
num_procs = 1200


totNus = 0
totElems = 0 ### this is the dumped number of elems
len_time = 0
grandTotElems = 0 ### this is the total number of elems in the mesh

for i in range(num_procs):
    
    wvf_name = "wavefield_db_fwd_%s.nc4" % i
    f = Dataset(wvf_dir+wvf_name, 'r')
    totNus += np.sum(f.variables["Nus"][:])
    totElems += len(f.variables["Nus"][:])
    
    if i == 0:
        len_time = len(f.variables["displacement_wavefield"][:])
        grandTotElems = len(f.variables["domain_decomposition"][:])
        print(len_time)
    f.close()
print(totNus)
print(totElems)    
print(grandTotElems)
print(len_time)
f_merged = Dataset(wvf_dir + wvf_name_merged, "w", "NETCDF4")

'''
f_merged.createDimension("ncdim_%s"%len_time, len_time)
f_merged.createDimension("ncdim_%s"%totElems, totElems)
f_merged.createDimension("ncdim_%s"%totNus, totNus)
f_merged.createDimension("ncdim_5", 5)
f_merged.createDimension("ncdim_6", 6)
f_merged.createDimension("ncdim_%s"%grandTotElems, grandTotElems)
'''
f_merged.createDimension('time', len_time)
f_merged.createDimension('tot_elems', totElems)
f_merged.createDimension('tot_nus', totNus)
f_merged.createDimension('gll_points', 5)
f_merged.createDimension('comp', 6)
f_merged.createDimension('grand_tot_elems', grandTotElems)



'''
merged_wvf = f_merged.createVariable("displacement_wavefield", "f8", ("ncdim_%s"%len_time,"ncdim_%s"%totNus, "ncdim_6", "ncdim_5", "ncdim_5") )
merged_nus = f_merged.createVariable("Nus", "i", ("ncdim_%s"%totElems) )
merged_nrs = f_merged.createVariable("Nrs", "i", ("ncdim_%s"%totElems) )
merged_dd = f_merged.createVariable("domain_decomposition", "i", ("ncdim_%s"%grandTotElems) )
merged_s = f_merged.createVariable("mesh_S", "f8", ("ncdim_%s"%totElems, "ncdim_5", "ncdim_5") )
merged_z = f_merged.createVariable("mesh_Z", "f8", ("ncdim_%s"%totElems, "ncdim_5", "ncdim_5") )
merged_sem = f_merged.createVariable("sem_mesh", "i", ("ncdim_%s"%totElems, "ncdim_5", "ncdim_5") )
'''
merged_wvf = f_merged.createVariable("displacement_wavefield", "f8", ('time','tot_nus', 'comp', 'gll_points', 'gll_points') )
merged_nus = f_merged.createVariable("Nus", "i", ('tot_elems') )
merged_nrs = f_merged.createVariable("Nrs", "i", ('tot_elems') )
merged_dd = f_merged.createVariable("domain_decomposition", "i", ('grand_tot_elems') )
merged_s = f_merged.createVariable("mesh_S", "f8", ('tot_elems', 'gll_points', 'gll_points') )
merged_z = f_merged.createVariable("mesh_Z", "f8", ('tot_elems', 'gll_points', 'gll_points') )
merged_sem = f_merged.createVariable("sem_mesh", "i", ('tot_elems', 'gll_points', 'gll_points') )



pos_nu = 0
pos_elem = 0
for i in range(num_procs):
    
 
    wvf_name = "wavefield_db_fwd_%s.nc4" % i
    f = Dataset(wvf_dir+wvf_name, 'r')
    if i == 0:
        merged_dd[:] = f.variables["domain_decomposition"][:]   
    wvf = f.variables["displacement_wavefield"][:]
    nuproc =  np.sum(f.variables["Nus"][:])
    elemproc = len(f.variables["Nus"][:])
    print(wvf.shape)
    print(nuproc)
    if (elemproc != 0):

        merged_wvf[:, pos_nu:pos_nu + nuproc, :, :, :] = wvf
        merged_nus[pos_elem: pos_elem + elemproc] = f.variables["Nus"][:]
        merged_nrs[pos_elem: pos_elem + elemproc] = f.variables["Nrs"][:]
        merged_s[pos_elem: pos_elem + elemproc] = f.variables["mesh_S"][:]
        merged_z[pos_elem: pos_elem + elemproc] = f.variables["mesh_Z"][:]
        merged_sem[pos_elem: pos_elem + elemproc] = f.variables["sem_mesh"][:]

        pos_nu += nuproc
        pos_elem += elemproc
    f.close()
    
f_merged.close()
