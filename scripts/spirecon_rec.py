#!/usr/bin/env python3

from spirecon import *
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('f_name', help='SPIFFI file to reconstruct')
parser.add_argument('--config', help='Override config')
parser.add_argument('--shift', help='Shift the reconstructed image since this doesn\'t use the arc lamps', type=float, default=2)
parser.add_argument('--path', help='Directory to search for files', default="/home/ydallilar/Documents/SPIFFIER/DETDATA")
parser.add_argument('--out', help='Naming of the outputs', default="out")
parser.add_argument('--gfit', action='store_true')
args = parser.parse_args()

im_to_r = '/'.join([args.path, args.f_name])
config = args.config if args.config is not None else "%s_250mas" % "".join(fits.getheader(im_to_r)["HIERARCH ESO INS3 SPGW NAME"].split("_"))
out = args.out

d_output_tr = "trace_%s.txt" % config
d_output_sh = "wave_%s.txt" % config

im = fits.getdata("%s" % im_to_r)[4:-4,4:-4]
im = order(im, config, shft=2)
im = apply_shfts(im, config, out)
im = cube(im, config, header=fits.getheader(im_to_r), out_im=out)

if args.gfit:
    X, Y = gfit(im)
    print("GFIT -> X : %.1f, Y :%.1f" % (X, Y))
