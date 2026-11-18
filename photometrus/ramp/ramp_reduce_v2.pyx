import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import time
import numpy as np
from astropy.io import fits
import fitsio

cimport numpy as np
from cython.parallel import prange
cimport cython

#from openmp cimport omp_get_max_threads
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)

cdef inline double cheb_scale(double x, double xmax) nogil:
    return x * (2.0 / xmax) - 1

cdef double cheb_eval(double x, double *c, int deg) nogil:
    cdef:
        double b_next2 = 0.0
        double b_next1 = 0.0 
        double b 
        int i 
    for i in range(deg-1, -1, -1):
        b = 2.0 * x * b_next1 - b_next2 + c[i]
        b_next2 = b_next1
        b_next1 = b
    return b_next1 - x * b_next2


cdef double median_qs(double* a, int n) nogil:
    cdef int l = 0
    cdef int r = n - 1
    cdef int k = n // 2
    cdef int i, j
    cdef double x, t

    if n == 0:
        return 0.0

    while l < r:
        x = a[k]
        i = l
        j = r

        while i <= j:
            while a[i] < x: i += 1
            while a[j] > x: j -= 1
            if i <= j:
                t = a[i]; a[i] = a[j]; a[j] = t
                i += 1; j -= 1

        if j < k: l = i
        if k < i: r = j

    if n & 1:
        return a[k]
    else:
        x = a[0]
        for i in range(1, k):
            if a[i] > x:
                x = a[i]
        return 0.5 * (x + a[k])

# ---------------------------------------------------------------------

def do_rpc(
        np.ndarray[np.float64_t, ndim=2] data,
        np.ndarray[np.float64_t, ndim=2] superbias,
        int nframe, int wx=4096, int wy=4096, int rpc_width=128):
    cdef:
        np.ndarray[np.float64_t, ndim=2] sub_data = np.empty((wy, wx), dtype=np.float64)
        int nout = wx // rpc_width
        int l, r, i
    sub_data = data - superbias
    for i in range(nout):
        l, r = i * rpc_width, (i + 1) * rpc_width
        if nframe == 0:
            sub_data[:, l:r] -= np.nanmedian(sub_data[4092:4096, l:r])
        else:
            sub_data[:, l:r] -= np.nanmedian([sub_data[4092:4096, l:r], sub_data[0:4, l:r]])
    return sub_data


cdef class HMDReducer:
    cdef:
        double[:,:,::1] coe_C
        double[:,:,::1] coe_Ad
        double[:,:,::1] coe_Fd
        double[:,::1] satulim
        double[:,::1] darklim
        double[:,::1] superbias
        double[:,::1] Fdark_total
        double[:,::1] Ad1
        int nx, ny
        int C_deg, Ad_deg, Fd_deg

    def __init__(self,
                 double[:,:,::1] coe_C,
                 double[:,:,::1] coe_Ad,
                 double[:,::1] satulim,
                 double[:,::1] darklim,
                 double[:,::1] superbias):

        cdef int i, j

        self.coe_C = coe_C
        self.coe_Ad = coe_Ad
        self.satulim = satulim
        self.darklim = darklim
        self.superbias = superbias

        #mask = ~np.isfinite(self.superbias)
        #self.superbias[mask] = 0.0
        #self.satulim[mask] = 0.0
        #self.darklim[mask] = 0.0

        self.ny = coe_C.shape[0]
        self.nx = coe_C.shape[1]
        self.C_deg = coe_C.shape[2]
        self.Ad_deg = coe_Ad.shape[2]
        self.Fd_deg = coe_Ad.shape[2] + 1


        for i in range(self.ny):
            for j in range(self.nx):
                if not np.isfinite(self.superbias[i, j]):
                    self.superbias[i,j] = 0.0
                    self.satulim[i,j] = 0.0
                    self.darklim[i,j] = 0.0

        self.coe_Fd = np.zeros((self.ny, self.nx, self.Fd_deg), dtype=np.float64)
        self.Fdark_total = np.zeros((self.ny, self.nx), dtype=np.float64)
        self.Ad1 = np.zeros((self.ny, self.nx), dtype=np.float64)

        self._prep_dark()
        #print (self.coe_C[1000,1000,0])

    cdef _prep_dark(self):
        cdef:
            int j, i, n, N
            double *a
            double *b
            double an, Ad1, Fdm1, Fdp1, sign
            #int N = self.D_deg
        N = self.Ad_deg

        for j in range(4, self.ny-4):
            for i in range(4, self.nx-4):
                a = &self.coe_Ad[j, i, 0]
                b = &self.coe_Fd[j, i, 0]
                # --- integral & Ad(1) ----
                Ad1 = 0.0
                for n in range(N):
                    an = a[n]
                    Ad1 += an

                    if n == 0:
                        b[1] += an

                    elif n == 1:
                        b[2] += an * 0.25

                    else:
                        b[n+1] += an / (2.0*(n+1))
                        b[n-1] -= an / (2.0*(n-1))
                # --- Fd(1) & Fd(-1) ---
                Fdm1 = 0.0
                Fdp1 = 0.0
                sign = 1.0
                for n in range(N+1):
                    Fdp1 += b[n]
                    Fdm1 += sign * b[n]
                    sign = -sign
                # ---
                b[0] -= Fdm1
                self.Fdark_total[j,i] = Fdp1 - Fdm1
                self.Ad1[j,i] = Ad1

    @property
    def coe_Fd_view(self):
        return np.asarray(self.coe_Fd)


    cdef void subtract_bias(self, unsigned short[:,::1] raw, double[:,:] data, unsigned char[:,:] mask, int nf):
        cdef:
            double[:,::1] sb = self.superbias
            double[:,::1] sat = self.satulim
            int i, x, y, n, idx, l, block
            double tmp[1024]
            double med, val
            int width = 128
            int nout = 32
            double[32] medians

            unsigned short *p_raw
            double *p_sb
            double *p_sat
        with nogil:
            for block in range(nout):
                l = block * width
                n = 0
                # btm reference (4092:4096)
                for y in range(4092, 4096):
                    p_raw = &raw[y, l]
                    p_sb = &sb[y, l]
                    for x in range(width):
                        val =  <double>p_raw[x] - p_sb[x]
                        if val == val:
                            tmp[n] = val
                            n += 1
                # top reference
                if nf != 0:
                    for y in range(0, 4):
                        p_raw = &raw[y, l]
                        p_sb = &sb[y, l]
                        for x in range(width):
                            val =  <double>p_raw[x] - p_sb[x]
                            if val == val:
                                tmp[n] = val
                                n += 1
                medians[block] = median_qs(tmp, n)
            for y in range(4,4092):
                for block in range(nout):
                    med = medians[block]
                    l = block * width
                    p_raw = &raw[y,l]
                    p_sb = &sb[y,l]
                    p_sat = &sat[y,l]
                    for i in range(width):
                        val = p_raw[i] - (p_sb[i] + med)
                        data[y,l+i] = val
                        mask[y,l+i] = val > p_sat[i]


    cdef double[:,:,::1] _nlc_cds(self, double[:,:,::1] data_cube, unsigned char[:,:,::1] mask_cube):
        cdef:
            int i, x, y, z, n
            int nz
            double d, ds, dd, dd0
            double slim, dlim, fdt, ad1
            double slim_inv, dlim_inv
            int C_deg = self.C_deg
            int Fd_deg = self.Fd_deg
            double[:,:,::1] coe_C = self.coe_C
            double[:,:,::1] coe_Fd = self.coe_Fd
            double[:,::1] satulim = self.satulim
            double[:,::1] darklim = self.darklim
            double[:,::1] Fd_total = self.Fdark_total
            double[:,::1] Ad1 = self.Ad1
            double *p_d
            unsigned char *p_m
            double *p_C
            double *p_Fd

            double sc1, sc2, sc_tmp, sc, r
            double sd1, sd2, sd
            double fd0, fd_sum
            double sl, sl_unc_inv
            double gain_inv = 1.0 / 1.8
            double var_read = 25.0

            double[:,:,::1] cds_cube
            double sc0
            double *p_cds
        nz = data_cube.shape[2]
        cds_cube = np.zeros((self.ny, self.nx, nz))
        with nogil:
            for y in prange(4, 4092, num_threads=4):
                for x in range(4, 4092):
                    p_d = &data_cube[y,x,0]
                    p_m = &mask_cube[y,x,0]
                    p_C = &coe_C[y,x,0]
                    p_Fd = &coe_Fd[y,x,0]
                    slim = satulim[y,x]
                    dlim = darklim[y,x]
                    fdt = Fd_total[y,x]
                    ad1 = Ad1[y,x]
                    p_cds = &cds_cube[y,x,0]
                    if slim == 0 or dlim == 0:
                        continue
                    slim_inv = 2.0 / slim
                    dlim_inv = 2.0 / dlim
                    dd0 = -1.0
                    fd0 = 0.0
                    fd_sum = 0.0
                    sc0 = 0.0
                    for z in range(nz):
                        d = p_d[z]
                        if p_m[z] == 1 or d <= 0:
                            break
                        ds = slim_inv * d - 1.0    ## scaling
                        dd = dlim_inv * d - 1.0    ##
                        if dd == dd0:
                            break
                        sc1 = 0.0
                        sc2 = 0.0
                        for i in range(C_deg-1, -1, -1):
                            sc = 2.0 * ds * sc1 - sc2 + p_C[i]
                            sc2 = sc1
                            sc1 = sc
                        sc = (sc1 - ds * sc2) * slim *0.5
                        if dd >= 1.0:
                            sd = fdt + ad1 * (dd-1.0)
                        else:
                            sd1 = 0.0
                            sd2 = 0.0
                            for i in range(Fd_deg-1, -1, -1):
                                sd = 2.0 * ds * sd1 - sd2 + p_Fd[i]
                                sd2 = sd1
                                sd1 = sd
                            #sd = sd1 - ds * sd2 ## miss ds -> dd
                            sd = sd1 - dd * sd2
                        p_cds[z] = sc - sc0 - (sd - fd0) / (dd- dd0)
                        fd0 = sd
                        dd0 = dd
                        sc0 = sc
        return cds_cube


    cdef double[:,::1] _nlc_ramp(self, double[:,:,::1] data_cube, unsigned char[:,:,::1] mask_cube):
        cdef:
            int i, x, y, z, n
            int nz
            double d, ds, dd, dd0
            double slim, dlim, fdt, ad1
            double slim_inv, dlim_inv
            int C_deg = self.C_deg
            int Fd_deg = self.Fd_deg
            double[:,:,::1] coe_C = self.coe_C
            double[:,:,::1] coe_Fd = self.coe_Fd
            double[:,::1] satulim = self.satulim
            double[:,::1] darklim = self.darklim
            double[:,::1] Fd_total = self.Fdark_total
            double[:,::1] Ad1 = self.Ad1
            double *p_d
            unsigned char *p_m
            double *p_C
            double *p_Fd
            double *p_slim
            double *p_dlim
            double *p_Fdt
            double *p_ad1

            double sc1, sc2, sc_tmp, sc, r
            double sd1, sd2, sd
            double fd0, fd_sum
            double sl, sl_unc_inv
            double gain_inv = 1.0 / 1.8
            double var_read = 25.0

            double w, wx, wy, wxx, wxy

            double[:,::1] slope_data
            double slope
            double *p_slope
        nz = data_cube.shape[2]
        slope_data = np.zeros((self.ny, self.nx), dtype=np.float64)
        with nogil:
            for y in prange(4, 4092, num_threads=4):
                p_slim = &satulim[y,0]
                p_dlim = &darklim[y,0]
                p_Fdt = &Fd_total[y,0]
                p_ad1 = &Ad1[y,0]
                p_slope = &slope_data[y,0]
                for x in range(4, 4092):
                    p_d = &data_cube[y,x,0]
                    p_m = &mask_cube[y,x,0]
                    p_C = &coe_C[y,x,0]
                    p_Fd = &coe_Fd[y,x,0]
                    slim = p_slim[x]
                    dlim = p_dlim[x]
                    fdt = p_Fdt[x]
                    ad1 = p_ad1[x]
                    if slim == 0 or dlim == 0:
                        continue
                    slim_inv = 2.0 / slim
                    dlim_inv = 2.0 / dlim
                    dd0 = -1.0
                    fd0 = 0.0
                    fd_sum = 0.0
                    w = 0.0
                    wx = 0.0
                    wy = 0.0
                    wxx = 0.0
                    wxy = 0.0
                    for z in range(nz):
                        d = p_d[z]
                        if p_m[z] == 1 or d <= 0:
                            break
                        ds = slim_inv * d - 1.0    ## scaling
                        dd = dlim_inv * d - 1.0    ##
                        if dd == dd0:
                            break
                        sc1 = 0.0
                        sc2 = 0.0
                        for i in range(C_deg-1, -1, -1):
                            sc = 2.0 * ds * sc1 - sc2 + p_C[i]
                            sc2 = sc1
                            sc1 = sc
                        sc_tmp = sc1 - ds * sc2
                        if sc_tmp == 0.0:
                            break
                        r = sc_tmp / (ds+1.0)
                        sc = sc_tmp * slim * 0.5
                        if dd >= 1.0:
                            sd = fdt + ad1 * (dd-1.0)
                        else:
                            sd1 = 0.0
                            sd2 = 0.0
                            for i in range(Fd_deg-1, -1, -1):
                                sd = 2.0 * ds * sd1 - sd2 + p_Fd[i]
                                sd2 = sd1
                                sd1 = sd
                            #sd = sd1 - ds * sd2 ##miss ds -> dd
                            sd = sd1 - dd * sd2
                        fd_sum = fd_sum + (sd - fd0) / (dd - dd0)
                        fd0 = sd
                        dd0 = dd
                        sl = sc - fd_sum
                        sl_unc_inv = 1.0 / (abs(sc) * gain_inv + r * r * var_read)

                        w = w + sl_unc_inv
                        wx = wx + z * sl_unc_inv
                        wy = wy + sl * sl_unc_inv
                        wxx = wxx + z * z * sl_unc_inv
                        wxy = wxy + sl * z * sl_unc_inv
                    
                    if wx == 0.0:
                        slope = 65535.0
                    else:
                        slope = (w * wxy - wx * wy) / (w * wxx - wx * wx)
                    p_slope[x] = slope
        return slope_data

    def reduce(self, list raw_paths, int hdu, int mode):
        cdef:
            double[:,:,::1] data_cube
            unsigned char[:,:,::1] mask_cube
            int nf, nf_total, i
            unsigned short[:,::1] raw
            object tmp_header, header


        nf_total = len(raw_paths)
        data_cube = np.zeros((self.ny, self.nx, nf_total), dtype=np.float64)
        mask_cube = np.zeros((self.ny, self.nx, nf_total), dtype=np.uint8)

        for i, raw_path in enumerate(raw_paths):
            #st = time.time()
            try:
                with fitsio.FITS(raw_path) as fits:
                    raw = fits[hdu][:,6:6+4096]
                    tmp_header = fits[hdu].read_header()
            except FileNotFoundError:
                print (f"{raw_path} is not found")
                break
            #print (time.time()-st)
            nf = int(tmp_header["FRAME"])
            if nf == 0:
                header = tmp_header
            self.subtract_bias(raw, data_cube[:,:,i], mask_cube[:,:,i], nf)
            #data_cube[:,:,i] = data
            #mask_cube[:,:,i] = mask
        st = time.time()
        if mode == 0: # ramp
            return self._nlc_ramp(data_cube, mask_cube), header
        elif mode == 1: # cds
            return self._nlc_cds(data_cube, mask_cube), header
        else:
            raise ValueError("Invalid mode")
        print (time.time()-st)
