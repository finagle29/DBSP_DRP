"""
Command-line tool for tripping spliced spectrum to a specified wavelength range.
"""
import argparse
import os
from typing import Optional, List

from astropy.io import fits

def entrypoint():
    main(parse())

def parse(options: Optional[List[str]] = None) -> argparse.Namespace:
    argparser = argparse.ArgumentParser(description="Trim spliced spectrum to wavelength range.",
        formatter_class=argparse.RawTextHelpFormatter)

    argparser.add_argument("fname", type=str, nargs="+", help="Final data "
        "output from DBSP_DRP to trim. Can be multiple files.")

    argparser.add_argument("--low", type=float, default=3300,
        help="Min wavelength after trimming. Default is 3300 Å.")

    argparser.add_argument("--high", type=float, default=10500,
        help="Max wavelength after trimming. Default is 10500 Å.")

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

        hdul.close()
