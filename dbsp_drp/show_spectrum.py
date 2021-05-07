import argparse
import os
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
                if ext == 0 or ext >= len(hdul) or ext <= -len(hdul):
                    # check for negative oob
                    raise IndexError(f"Extension index {ext} out of range: "
                        f"must be between 1 and {len(hdul) - 1} inclusive "
                        f"(or between -1 and {-len(hdul)+1} inclusive to "
                        "index from the end).")
            except ValueError:
                raise LookupError(f"Extension '{args.extension}' not found in "
                    f"{args.fname}, and cannot be cast to an integer.\n"
                    f"\tValid extensions present in {args.fname} are {exts}.")
        spectrum = hdul[ext].data
        basename = os.path.splitext(os.path.basename(args.fname))[0]
        plot(spectrum, f"{basename}[{hdul[ext].name}]")

def plot(spec: fits.FITS_rec, title: str) -> None:
    """
    Plots spectrum and error with sensible y-scale limits.

    Args:
        spec (fits.FITS_rec): Spectrum to plot
    """
    plt.step(spec['wave'], spec['flux'], c='k', label='spectrum')
    plt.step(spec['wave'], spec['sigma'], c='gray', label='error')

    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"Flux ($10^{-17}\mathrm{erg}/\mathrm{s}/\mathrm{cm}^2/\AA$)")

    top1 = np.abs(np.percentile(spec['flux'], 95)) * 1.5
    top2 = np.max(spec['flux'][spec['wave'] > 4000]) * 1.1
    top = max(top1, top2)
    bottom = -0.05 * top

    plt.ylim(bottom, top)
    plt.legend()
    plt.title(title)
    plt.show()
