import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve, find_peaks, find_peaks_cwt
from scipy.interpolate import interp1d, RectBivariateSpline
from numpy.fft import fft, fftfreq
from astropy.visualization import imshow_norm, ZScaleInterval

import pkg_resources
import sys
import os

locs = [9, 8, 10, 7, 11, 6, 12, 5, 13, 4, 14, 3, 15, 2, 16,
        1, 32,
        17, 31, 18, 30, 19, 29, 20, 28, 21, 27, 22, 26, 23, 25, 24]

locs = np.array(locs) - 1

calib_path = "calib/"
check_path = "check/"

def set_calib_path(path):
    global calib_path
    calib_path = path

def fetch_calib_file(config, opt):
    f_name = pkg_resources.resource_filename("spirecon", "calib/" + "%s_%s.txt" % (opt, config))
    if os.path.isfile(f_name):
        return f_name
    elif os.path.isfile(calib_path + "%s_%s.txt" % (opt, config)):
        return calib_path + "%s_%s.txt" % (opt, config)
    else:
        raise ValueError("Can't find the configuration file %s_%s.txt" % (opt, config))

def save_trace(config, coeffs):
    f = open(calib_path + "trace_%s.txt" % config, "w")
    for i in range(96):
        f.write("%02d %10.4e %10.4e %10.4e\n" % (i, *coeffs[i]))
    f.close()

def save_shfts(config, shfts):
    f = open(calib_path + "wave_%s.txt" % config, "w")
    for i in range(32):
        f.write("%02d %10.2f\n" % (i, shfts[i]))
    f.close()

def read_trace(config):
    f_name = fetch_calib_file(config, 'trace')
    data = np.loadtxt(f_name)
    return data[:,1:]

def read_shfts(config):
    f_name = fetch_calib_file(config, 'wave')
    data = np.loadtxt(f_name)
    return data[:,1:]


def trace_LS(im, config, d_plot=None, show=False, d_output=None, out_im="test.fits", bs=40, p_height=1000):

    cls = ['b', 'y', 'r']
    
    b_im = np.median(im.reshape([2040//bs, bs, 2040]), axis=1)
    peaks_2d = np.zeros([2040//bs, 96]) - 1
    imshow_norm(im, cmap="gray", interval=ZScaleInterval())
    
    for i in range(2040//bs)[::-1]:

        peaks = find_peaks(b_im[i,:], height=p_height)[0]
        diffs = peaks[1:] - peaks[:-1]
        ndxs = [0]
        cntr = 0
        for j in range(len(peaks)-1):
            cntr += np.rint(diffs[j]/21).astype(int)
            ndxs.append(cntr)
        peaks_2d[i, ndxs] = peaks
    
    coeffs = []
    
    y = np.arange(0, 2040)
    
    for i in range(96):

        yy = np.arange(2040//bs)*bs+0.5*bs
        peaks = peaks_2d[:,i]
        ndx = np.where(peaks == -1)
        peaks = np.delete(peaks, ndx)
        yy = np.delete(yy, ndx)
        plt.plot(peaks, yy, cls[i%3] + '.', markersize=0.2)

        pols = np.polyfit(yy, peaks, 2)
        coeffs.append(pols)
        fun = np.poly1d(pols)
        plt.plot(fun(y), y, cls[i%3] + '-', linewidth=0.2)
    
    plt.gca().invert_yaxis()
    if d_plot is not None:
        plt.axis("off")
        plt.savefig(d_plot, bbox_inches="tight", dpi=300)
    if show:
        plt.show()

    save_trace(config, coeffs)

def order(im, config, out_sz=64, in_sz=102, shft=2, out_im=None):

    coeffs = read_trace(config)

    r_im = np.zeros([2040, out_sz*32])

    for i in range(32):

        tmp = np.zeros([2040, in_sz])

        t_coeff = coeffs[i*3+1,:]
        polf = np.poly1d(t_coeff)
        #poldf = np.poly1d([t_coeff[0]*2, t_coeff[1]])
        Y = np.arange(2040)
        X = np.arange(2040)

        lim_l, lim_h = int(t_coeff[-1])-in_sz//2, int(t_coeff[-1])+in_sz//2

        if lim_l < 0:
            lim_l = 0
        if lim_h > 2047:
            lim_h = 2048

        fun = RectBivariateSpline(np.arange(2040), X[lim_l:lim_h], im[:, lim_l:lim_h], kx=1, ky=1)
        XX, YY = np.meshgrid(np.arange(in_sz), Y)
        tmp[:,:] =  fun(YY, polf(YY)+(XX-in_sz//2)+2, grid=False)
        r_im[:,locs[i]*out_sz:(locs[i]+1)*out_sz] = tmp[:,(in_sz-out_sz)//2:-(in_sz-out_sz)//2]

    if out_im is not None:
        fits.PrimaryHDU(r_im).writeto("ord_%s.fits" % out_im, overwrite=True)

    return r_im

def calc_shfts(im, config):

    out_sz = 64

    shfts = np.zeros(32).astype(int)

    for i in range(32):
        colr = np.median(im[:,1*out_sz:2*out_sz], axis=1)
        col = np.median(im[:,locs[i]*out_sz:(locs[i]+1)*out_sz], axis=1)

        conv = fftconvolve(col, colr[::-1], mode='same')
        shfts[i] = (np.argmax(conv)-1020)

    print(shfts)

    save_shfts(config, shfts)

def apply_shfts(im, d_output, out_im="test.fits"):

    out_sz = 64
    shfts = read_shfts(d_output).astype(int)

    r_im = np.zeros([2040+(np.max(shfts)-np.min(shfts)).astype(int), out_sz*32]) + np.NaN
    minshft = np.max(shfts).astype(int)

    for i in range(32):
        #print(-shfts[i]+minshft, 2040-shfts[i]+minshft, (locs[i]+1)*out_sz)
        r_im[int(-shfts[i]+minshft):int(2040-shfts[i]+minshft),locs[i]*out_sz:(locs[i]+1)*out_sz] = im[:,locs[i]*out_sz:(locs[i]+1)*out_sz]

    fits.PrimaryHDU(r_im).writeto("shft_%s.fits" % out_im, overwrite=True)

    return r_im

def cube(im, d_output, out_im="test.fits"):

    shfts = read_shfts(d_output)
    mn_shft, mx_shft = int(np.min(shfts)), int(np.max(shfts))

    sz = im.shape

    tmp = np.reshape(im, [sz[0], 32, 64])
    tmp = np.repeat(tmp, 2, axis=1)
    tmp = np.reshape(tmp, [sz[0], 64, 64])

    fits.PrimaryHDU(tmp).writeto("cube_%s.fits" % out_im, overwrite=True)

    tmp = np.mean(im[mx_shft-mn_shft:sz[1]-(mx_shft-mn_shft),:], axis=0)
    tmp = np.reshape(tmp, [32, 64])
    tmp = np.repeat(tmp, 2, axis=0)
    tmp = np.reshape(tmp, [64, 64])

    fits.PrimaryHDU(tmp).writeto("rec_%s.fits" % out_im, overwrite=True)

    return tmp

def gfit(r_im):

    from scipy.optimize import curve_fit

    def gauss2d(X, a, b, c, d, e):
        return d*np.exp(-((X[0]-a)**2+(X[1]-b)**2)*0.5/c**2)+e

    X = np.arange(1, 65)
    XX, YY = np.meshgrid(X, X)

    X, Y = np.unravel_index(np.argmax(r_im), (64, 64))

    coo = np.vstack([np.ravel(XX), np.ravel(YY)])

    out = curve_fit(gauss2d, coo, r_im.ravel(), p0=[Y, X, 2, np.max(r_im), 0])

    return out[0][0], out[0][1]

def arc_shft(ref_im, im):

    from scipy.signal import convolve, fftconvolve, convolve2d
    from scipy.optimize import curve_fit
    import pyds9

    def gauss2d(X, a, b, c, d, e, f):
        return e*np.exp(-((X[0]-a)**2/c**2+(X[1]-b)**2/d**2)*0.5)+f

    #im = np.roll(im, 10, axis=1)
    conv = convolve(im, ref_im[::-1,::-1], mode="same")


    X = np.arange(900, 1100)
    XX, YY = np.meshgrid(X, X)

    coo = np.vstack([np.ravel(XX), np.ravel(YY)])

    X, Y = np.unravel_index(np.argmax(conv), (2048, 2048))

    out = curve_fit(gauss2d, coo, conv[900:1100,900:1100].ravel(), p0=[X, Y, 10, 50, np.max(conv), 0])

    return(out[0][0]-1024, out[0][1]-1024)
