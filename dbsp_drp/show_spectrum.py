import argparse
from typing import List, Optional

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

def parser(options: Optional[List[str]] = None) -> argparse.Namespace:
    argparser = argparse.ArgumentParser(description="Script to plot DBSP spectra",
        formatter_class=argparse.RawTextHelpFormatter)

    argparser.add_argument("fname", type=str, help="path to target_a.fits file")

    argparser.add_argument("--extension", type=str, default="SPLICED",
                           help="Extension name or number")

    return argparser.parse_args() if options is None else argparser.parse_args(options)

def main(args: argparse.Namespace) -> None:
    with fits.open(args.fname) as hdul:
        exts = [hdu.name for hdu in hdul if hdu.name != "PRIMARY"]
        if args.extension.upper() in exts:
            ext = args.extension
        else:
            try:
                ext = int(args.extension)
                if ext == 0 or ext >= len(hdul):
                    raise IndexError(f"Extension index {ext} out of range: "
                        f"must be between 1 and {len(hdul) - 1} inclusive.")
            except ValueError:
                raise LookupError(f"Extension '{args.extension}' not found in "
                    f"{args.fname}, and cannot be cast to an integer.\n"
                    f"\tValid extensions present in {args.fname} are {exts}.")
        spectrum = hdul[ext].data
        plot(spectrum)

def plot(spec: fits.FITS_rec) -> None:
    """
    Plots spectrum and error with sensible y-scale limits.

    Args:
        spec (fits.FITS_rec): Spectrum to plot
    """
    plt.plot(spec['wave'], spec['flux'], c='k', label='spectrum')
    plt.plot(spec['wave'], spec['sigma'], c='gray', label='error')

    top = np.abs(np.percentile(spec['flux'], 95)) * 1.5
    bottom = -0.1 * top

    plt.ylim(bottom, top)
    plt.legend()
    plt.show()
