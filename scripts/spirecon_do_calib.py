#!/usr/bin/env python3

from spirecon import trace_LS, read_trace, order, calc_shfts
import sys
from astropy.io import fits

config = sys.argv[1]
trace_im = fits.getdata(sys.argv[2])[4:-4,4:-4]
wave_im = fits.getdata(sys.argv[3])[4:-4,4:-4]

trace_LS(trace_im, config, d_plot="trace.png")
trace = read_trace(config)
order(trace_im, config)

wave_im_ord = order(wave_im, config, out_im="test")

calc_shfts(wave_im_ord, config)
