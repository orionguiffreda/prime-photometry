import os

os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

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
#@cython.nonecheck(False)
#@cython.cdivision(True)

cdef void calc_ramp(
        np.ndarray[np.float64_t, ndim=2] ramp,
        np.ndarray[np.float64_t, ndim=2] sum_w,
        np.ndarray[np.float64_t, ndim=2] sum_wx,
        np.ndarray[np.float64_t, ndim=2] sum_wy,
        np.ndarray[np.float64_t, ndim=2] sum_wxx,
        np.ndarray[np.float64_t, ndim=2] sum_wxy):
    cdef:
        int wy = <int> ramp.shape[0]
        int wx = <int> ramp.shape[1]
        int i, j
        double deno
    with nogil:
        for i in prange(4, wy - 4, num_threads=8):
            for j in range(4, wx - 4):
                deno = sum_w[i, j] * sum_wxx[i, j] - sum_wx[i, j] ** 2
                if deno == 0.0:
                    ramp[i, j] = 65535.0
                else:
                    ramp[i, j] = (sum_w[i, j] * sum_wxy[i, j] - sum_wx[i, j] * sum_wy[i, j]) / deno

cdef void calc_sum(
        np.ndarray[np.float64_t, ndim=2] corr_data,
        np.ndarray[np.uint8_t, ndim=2] tmp_x,
        np.ndarray[np.float64_t, ndim=2] sum_w,
        np.ndarray[np.float64_t, ndim=2] sum_wx,
        np.ndarray[np.float64_t, ndim=2] sum_wy,
        np.ndarray[np.float64_t, ndim=2] sum_wxx,
        np.ndarray[np.float64_t, ndim=2] sum_wxy,
        np.ndarray[np.uint8_t, ndim=2] mask,
        double gain_inv,
        double var_read,
        bint weight):
    cdef:
        int wy = <int> sum_w.shape[0]
        int wx = <int> sum_w.shape[1]
        int i, j
        double w, tx, cd
    with nogil:
        for i in prange(4, wy - 4, num_threads=8):
            for j in range(4, wx - 4):
                tx = tmp_x[i, j]
                cd = corr_data[i, j]
                if weight:
                    w = mask[i, j] / (cd * gain_inv + var_read)
                else:
                    w = mask[i, j]

                sum_w[i, j] += w
                sum_wx[i, j] += w * tx
                sum_wy[i, j] += w * cd
                sum_wxx[i, j] += w * tx ** 2
                sum_wxy[i, j] += w * tx * cd

cdef void subtract_dark(
        np.ndarray[np.float64_t, ndim=2] corr_data,
        np.ndarray[np.float64_t, ndim=2] Fdark,
        np.ndarray[np.float64_t, ndim=2] tmp_Fdark,
        np.ndarray[np.float64_t, ndim=2] rpc,
        np.ndarray[np.float64_t, ndim=2] tmp_rpc,
        np.ndarray[np.float64_t, ndim=2] total_dark,
        np.ndarray[np.uint8_t, ndim=2] mask):
    cdef:
        int wy = <int> corr_data.shape[0]
        int wx = <int> corr_data.shape[1]
        int i, j
        double deno
    with nogil:
        for i in prange(4, wy - 4, num_threads=8):
            for j in range(4, wx - 4):
                deno = rpc[i, j] - tmp_rpc[i, j]
                if deno != 0.0:
                    total_dark[i, j] += (Fdark[i, j] - tmp_Fdark[i, j]) / deno
                tmp_rpc[i, j] = rpc[i, j]
                tmp_Fdark[i, j] = Fdark[i, j]
                corr_data[i, j] -= total_dark[i, j] * mask[i, j]

def calc_darklim(
        np.ndarray[np.float64_t, ndim=3] coe_D,
        np.ndarray[np.float64_t, ndim=2] darklim,
        int wx=4096, int wy=4096):
    cdef:
        int D_deg = <int> coe_D.shape[2]
        int i, j, k
        np.ndarray[np.float64_t, ndim=2] Adarklim = np.empty((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] Fdarklim = np.empty((wy, wx), dtype=np.float64)
        double dl, ad, fd, coe_di

    with nogil:
        for j in prange(4, wy - 4, num_threads=8):
            for k in range(4, wx - 4):
                dl = darklim[j, k]
                fd = 0.0
                ad = 0.0
                for i in range(D_deg - 1, -1, -1):
                    coe_di = coe_D[j, k, i]
                    ad = coe_di + ad * dl
                    fd = (fd + coe_di / (i + 1)) * dl
                Adarklim[j, k] = ad
                Fdarklim[j, k] = fd
    return Adarklim, Fdarklim

def do_nlc(
        np.ndarray[np.float64_t, ndim=2] data,
        np.ndarray[np.float64_t, ndim=3] coe_R,
        np.ndarray[np.float64_t, ndim=3] coe_D,
        np.ndarray[np.float64_t, ndim=2] darklim,
        int wx=4096, int wy=4096):
    cdef:
        int R_deg = <int> coe_R.shape[2]
        int D_deg = <int> coe_D.shape[2]
        int i, j, k
        np.ndarray[np.float64_t, ndim=2] Fdark = np.empty((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] corr_data = np.empty((wy, wx), dtype=np.float64)
        double d, dl
        double cd, fd
    with nogil:
        for j in prange(4, wy - 4, num_threads=8):
            for k in range(4, wx - 4):
                d = data[j, k]
                dl = darklim[j, k]
                cd = 0.0
                fd = 0.0
                if R_deg > D_deg:
                    for i in range(R_deg - 1, -1, -1):
                        cd = (cd + coe_R[j, k, i]) * d
                        if i < D_deg:
                            if d > dl:
                                fd = (fd + coe_D[j, k, i] / (i + 1)) * dl
                            else:
                                fd = (fd + coe_D[j, k, i] / (i + 1)) * d
                else:
                    for i in range(D_deg - 1, -1, -1):
                        if d > dl:
                            fd = (fd + coe_D[j, k, i] / (i + 1)) * dl
                        else:
                            fd = (fd + coe_D[j, k, i] / (i + 1)) * d
                        if i < R_deg:
                            cd = (cd + coe_R[j, k, i]) * d
                corr_data[j, k] = cd
                Fdark[j, k] = fd

    return corr_data, Fdark

def do_mask(
        np.ndarray[np.float64_t, ndim=2] rpc,
        np.ndarray[np.float64_t, ndim=2] satulim,
        np.ndarray[np.uint8_t, ndim=2] mask,
        np.ndarray[np.uint8_t, ndim=2] satu_mask):
    cdef:
        int wy = <int> rpc.shape[0]
        int wx = <int> rpc.shape[1]
        int i, j, m

    with nogil:
        for i in prange(4, wy - 4, num_threads=8):
            for j in range(4, wx - 4):
                if rpc[i, j] > satulim[i, j]:
                    satu_mask[i, j] = 0
        for i in prange(4, wy - 4, num_threads=8):
            for j in range(4, wx - 4):
                m = mask[i, j] & satu_mask[i, j] & satu_mask[i - 1, j] & satu_mask[i + 1, j] & satu_mask[i, j - 1] & \
                    satu_mask[i, j + 1]
                rpc[i, j] *= m

    return satu_mask, rpc

cdef void make_mask(
        np.ndarray[np.float64_t, ndim=2] rpc,
        np.ndarray[np.float64_t, ndim=2] satulim,
        np.ndarray[np.uint8_t, ndim=2] mask,
        np.ndarray[np.uint8_t, ndim=2] satu_mask,
        np.ndarray[np.uint8_t, ndim=2] master_mask,
        np.ndarray[np.uint8_t, ndim=2] tmp_x):
    cdef:
        int wy = <int> rpc.shape[0]
        int wx = <int> rpc.shape[1]
        int i, j, m

    with nogil:
        for i in prange(4, wy - 4, num_threads=8):
            for j in range(4, wx - 4):
                if rpc[i, j] > satulim[i, j]:
                    satu_mask[i, j] = 0
        for i in prange(4, wy - 4, num_threads=8):
            for j in range(4, wx - 4):
                m = satu_mask[i, j] & satu_mask[i - 1, j] & satu_mask[i + 1, j] & satu_mask[i, j - 1] & satu_mask[
                    i, j + 1]
                master_mask[i, j] = m
                rpc[i, j] *= m
                tmp_x[i, j] += m

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

def do_ramp(list raw_paths,
            np.ndarray[np.float64_t, ndim=2] superbias,
            np.ndarray[np.float64_t, ndim=2] satulim,
            np.ndarray[np.uint8_t, ndim=2] mask,
            np.ndarray[np.float64_t, ndim=3] coe_R_cube,
            np.ndarray[np.float64_t, ndim=3] coe_D_cube,
            np.ndarray[np.float64_t, ndim=2] darklim,
            int hdu_idx,
            double gain_inv = 1.0 / 1.8, double var_read=25.0,
            int wx=4096, int wy=4096,
            bint weight=True):
    cdef:
        np.ndarray[np.float64_t, ndim=2] ramp = np.zeros((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] raw_data = np.empty((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] corr_data = np.empty((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] Fdark = np.empty((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] rpc = np.empty((wy, wx), dtype=np.float64)
        np.ndarray[np.uint8_t, ndim=2] satu_mask = np.ones((wy, wx), dtype=np.uint8)
        np.ndarray[np.uint8_t, ndim=2] master_mask = np.ones((wy, wx), dtype=np.uint8)
        np.ndarray[np.float64_t, ndim=2] tmp_rpc = np.zeros((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] tmp_Fdark = np.zeros((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] total_dark = np.zeros((wy, wx), dtype=np.float64)
        np.ndarray[np.uint8_t, ndim=2] tmp_x = np.zeros((wy, wx), dtype=np.uint8)
        np.ndarray[np.float64_t, ndim=2] sum_w = np.zeros((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] sum_wx = np.zeros((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] sum_wy = np.zeros((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] sum_wxx = np.zeros((wy, wx), dtype=np.float64)
        np.ndarray[np.float64_t, ndim=2] sum_wxy = np.zeros((wy, wx), dtype=np.float64)
        int nframe
    ts = time.time()
    for raw_path in raw_paths:
        try:
            with fitsio.FITS(raw_path) as fits:
                raw_data = fits[hdu_idx].read()[:, 6:4096 + 6].astype(np.float64)
                tmp_header = fits[hdu_idx].read_header()
            #hdu = fits.open(raw_path)
        except FileNotFoundError:
            print (f"{raw_path} is not found")
            break
        #nframe = int(hdu[hdu_idx].header["FRAME"])
        nframe = int(tmp_header["FRAME"])
        if nframe == 0:
            #header = hdu[hdu_idx].header
            header = tmp_header
        #raw_data = hdu[hdu_idx].data[:,6:4096+6].astype(np.float64)
        rpc = do_rpc(raw_data, superbias, nframe)
        make_mask(rpc, satulim, mask, satu_mask, master_mask, tmp_x)
        corr_data, Fdark = do_nlc(rpc, coe_R_cube, coe_D_cube, darklim)
        #do_nlc_void(corr_data, Fdark, rpc, coe_R_cube, coe_D_cube, darklim, wx, wy)

        subtract_dark(corr_data, Fdark, tmp_Fdark, rpc, tmp_rpc, total_dark, master_mask)
        calc_sum(corr_data, tmp_x, sum_w, sum_wx, sum_wy, sum_wxx, sum_wxy, master_mask, gain_inv, var_read, weight)
    calc_ramp(ramp, sum_w, sum_wx, sum_wy, sum_wxx, sum_wxy)
    print (f"end loop: {time.time() - ts}s")
    return header, ramp
