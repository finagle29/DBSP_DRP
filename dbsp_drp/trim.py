import argparse
import os
from typing import Optional, List

import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox, Button, CheckButtons

from astropy.io import fits

from dbsp_drp.splicing import adjust_and_combine_overlap
from dbsp_drp import show_spectrum

def parse(options: Optional[List[str]] = None) -> argparse.Namespace:
    argparser = argparse.ArgumentParser(description="Trim spliced spectrum to wavelength range.",
        formatter_class=argparse.RawTextHelpFormatter)

    argparser.add_argument("fname", type=str, nargs="+", help="Final data "
        "output from DBSP_DRP to trim. Can be multiple files.")

    argparser.add_argument("--low", type=float, default=3300,
        help="Min wavelength after trimming")

    argparser.add_argument("--high", type=float, default=10500,
        help="Max wavelength after trimming")

    return argparser.parse_args() if options is None else argparser.parse_args(options)

def main(args: argparse.Namespace):
    for fname in args.fname:
        if not os.path.isfile(fname):
            raise FileNotFoundError(f"File {fname} not found!")

        hdul = fits.open(fname, mode="update")

        wave = hdul['SPLICED'].data['wave']
        flux = hdul['SPLICED'].data['flux']
        err = hdul['SPLICED'].data['sigma']

        mask = (wave > args.low) & (wave < args.high)

        hdul['SPLICED'].data['wave'] = wave[mask]
        hdul['SPLICED'].data['flux'] = flux[mask]
        hdul['SPLICED'].data['sigma'] = err[mask]
