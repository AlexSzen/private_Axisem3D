
import numpy as np

def spz2xyz(s, z, phi) :
    r = np.sqrt(np.power(s,2.)+np.power(z,2.))
    if (r<1e-6) :
         theta = 0.
    else :
        theta = np.arccos(z/r)

    x = r * np.sin((theta)) * np.cos(np.deg2rad(phi))
    y = r * np.sin((theta)) * np.sin(np.deg2rad(phi))
    z = r * np.cos((theta))

    return x, y, z


class WavefieldComputer:

    def __init__(self, wavefield, nu, nu2, s, z):
        self.wvf = wavefield
        self.nu = nu
        self.nu2 = nu2 #for nu forward
        self.s = s
        self.z = z
        self.nu_sum = np.cumsum(np.concatenate(([0],nu)))[:-1]
        self.r = np.sqrt(np.power(s,2.)+np.power(z,2.))
        self.num_steps = len(wavefield)
    
    def compute_non_fourier_slice(self, int_factor): #computes slices for real quantity defined on elements and points

        #just plot the whole slice at phi = 0
        r_inner = 0.
        r_outer = 6.371e6
        phi = 0.
        r_outer_closest = self.r[np.unravel_index(np.abs(self.r-r_outer).argmin(),self.r.shape)]
        r_inner_closest = self.r[np.unravel_index(np.abs(self.r-r_inner).argmin(),self.r.shape)]

        intfact_slice = []
        x_slice = []
        y_slice = []
        z_slice = []

        num_elems = len(self.nu)
        nPntEdge = 5

        for ielem in range(num_elems):
            for ipol in range(nPntEdge):
                for jpol in range (nPntEdge):

                    ### If within range, create the point
                    if self.r[ielem, ipol, jpol] > r_inner_closest and self.r[ielem,ipol,jpol] < r_outer_closest:
                        x_p, y_p, z_p = spz2xyz(self.s[ielem, ipol, jpol],self.z[ielem, ipol, jpol], phi)
                        x_slice.append(x_p)
                        y_slice.append(y_p)
                        z_slice.append(z_p)
                        intfact_slice.append(int_factor[ielem, ipol, jpol])

        return np.asarray(x_slice), np.asarray(y_slice), np.asarray(z_slice), np.asarray(intfact_slice)



    def compute_nu_slice(self):

        #just plot the whole slice at phi = 0
        r_inner = 0.
        r_outer = 6.371e6
        phi = 0.
        #print(np.abs(self.r-r_outer).argmin())
        r_outer_closest = self.r[np.unravel_index(np.abs(self.r-r_outer).argmin(),self.r.shape)]
        r_inner_closest = self.r[np.unravel_index(np.abs(self.r-r_inner).argmin(),self.r.shape)]

        nu_slice = []
        x_slice = []
        y_slice = []
        z_slice = []

        num_elems = len(self.nu)
        nPntEdge = 5

        for ielem in range(num_elems):
            for ipol in range(nPntEdge):
                for jpol in range (nPntEdge):

                    ### If within range, create the point
                    if self.r[ielem, ipol, jpol] > r_inner_closest and self.r[ielem,ipol,jpol] < r_outer_closest:
                        x_p, y_p, z_p = spz2xyz(self.s[ielem, ipol, jpol],self.z[ielem, ipol, jpol], phi)
                        x_slice.append(x_p)
                        y_slice.append(y_p)
                        z_slice.append(z_p)
                        nu_slice.append(self.nu[ielem])

        return np.asarray(x_slice), np.asarray(y_slice), np.asarray(z_slice), np.asarray(nu_slice)



    def compute_slice(self, r_inner, r_outer, phi, comp):

        r_outer_closest = self.r[np.unravel_index(np.abs(self.r-r_outer).argmin(),self.r.shape)]
        r_inner_closest = self.r[np.unravel_index(np.abs(self.r-r_inner).argmin(),self.r.shape)]

        ### hacky numpy solution for wvf
        wvf_slice = np.zeros((self.num_steps,1), dtype = np.float32)
        x_slice = []
        y_slice = []
        z_slice = []
        count = 0

        num_elems = len(self.nu)
        nPntEdge = 5

        for ielem in range(num_elems):
            for ipol in range(nPntEdge):
                for jpol in range (nPntEdge):

                    ### If within range, create the point
                    if self.r[ielem, ipol, jpol] > r_inner_closest and self.r[ielem,ipol,jpol] < r_outer_closest:

                        nuelem = self.nu[ielem] #nu of this elem
                        nuelem2 = self.nu2[ielem] #nu of this elem
                        ind = self.nu_sum[ielem] #index in elem and fouriers dimension of the wavefield
                        x_p, y_p, z_p = spz2xyz(self.s[ielem, ipol, jpol],self.z[ielem, ipol, jpol], phi)
                        x_slice.append(x_p)
                        y_slice.append(y_p)
                        z_slice.append(z_p)

                        ### hacky solution for the wavefield
                        if (count == 0):
                            wvf_slice[:, count] = self.wvf[:, ind, comp, ipol, jpol]
                        else:
                            if self.wvf[:, ind, comp, ipol, jpol].shape == (1,): ###another hacky thing if only one time step
                                wvf_slice = np.append(wvf_slice, np.expand_dims(self.wvf[:, ind, comp, ipol, jpol], axis=0), axis = 1 )
                            else:
                                wvf_slice = np.append(wvf_slice, self.wvf[:, ind, comp, ipol, jpol], axis = 1 )


                        for ialpha in range(1, max(nuelem,nuelem2)):
                            expval =  2 * np.exp(np.deg2rad(phi)*ialpha*1j)
                            wvf_slice[:, count] += np.real(expval * (self.wvf[:,ind+ialpha,comp,ipol,jpol] + 1j * self.wvf[:, ind+ialpha, comp+1, ipol, jpol]))

                        count += 1

        return np.asarray(x_slice), np.asarray(y_slice), np.asarray(z_slice), wvf_slice

    def compute_shell(self, r_shell, r_tolerance, sample_density, phis, inner_outer, comp):
        

        r_inner_closest = r_shell-r_tolerance
        r_outer_closest = r_shell+r_tolerance
#        r_outer_closest = self.r[np.unravel_index(np.abs(self.r-(r_shell + r_tolerance)).argmin(),self.r.shape)]
#        r_inner_closest = self.r[np.unravel_index(np.abs(self.r-(r_shell - r_tolerance)).argmin(),self.r.shape)]

        wvf_shell = np.zeros((self.num_steps, 1), dtype = np.float32)
        x_shell = []
        y_shell = []
        z_shell = []

        count = 0

        num_elems = len(self.nu)
        nPntEdge = 5

        for ielem in range(num_elems):
            for ipol in range(nPntEdge):
                for jpol in range (nPntEdge):

                    ### If within range, create the point
                    if self.r[ielem, ipol, jpol] > r_inner_closest and self.r[ielem,ipol,jpol] < r_outer_closest:

                        nuelem = self.nu[ielem] #nu of this elem
                        nuelem2 = self.nu2[ielem] #nu of this elem
                        ind = self.nu_sum[ielem] #index in elem and fouriers dimension of the wavefield
                        perimeter = ((phis[1]-phis[0])/360.) * 2 * np.pi * self.s[ielem, ipol, jpol] * 1e-6 #in thousand km
                        n_shell_points = (perimeter*sample_density).astype(int)

                        for isamp in range(n_shell_points):
                            phi = phis[0] + isamp * (phis[1]-phis[0]) / n_shell_points

                            x_p, y_p, z_p = spz2xyz(self.s[ielem, ipol, jpol],self.z[ielem, ipol, jpol], phi)
                            x_shell.append(x_p)
                            y_shell.append(y_p)
                            z_shell.append(z_p)

                            if count == 0:
                                wvf_shell[:, count] = self.wvf[:,ind,comp, ipol, jpol]
                            else:
                                if self.wvf[:, ind, comp, ipol, jpol].shape == (1,): ###another hacky thing if only one time step
                                    wvf_shell = np.append(wvf_shell, np.expand_dims(self.wvf[:, ind, comp, ipol, jpol], axis=0), axis = 1 )
                                else:
                                    wvf_shell = np.append(wvf_shell, self.wvf[:, ind, comp, ipol, jpol], axis = 1 )

                            for ialpha in range(1,max(nuelem, nuelem2)): #try truncating at something like nuFwd

                                expval = 2 * np.exp(np.deg2rad(phi)*ialpha*1j)
                                wvf_shell[:, count] += np.real( expval * (self.wvf[:,ind+ialpha,comp, ipol, jpol] + 1j * self.wvf[:,ind+ialpha,comp+1, ipol, jpol]))

                            count += 1

        return np.asarray(x_shell), np.asarray(y_shell), np.asarray(z_shell), wvf_shell
