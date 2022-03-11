"""
All-purpose interactive spectrum plotter for DBSP_DRP and intermediate PypeIt files.
"""
import argparse
import os
import sys
from typing import List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

from pypeit.specobjs import SpecObjs
from pypeit.onespec import OneSpec
from pypeit import msgs
from pypeit.pypmsgs import PypeItError

def entrypoint():
    main(parser())

def parser(options: Optional[List[str]] = None) -> argparse.Namespace:
    argparser = argparse.ArgumentParser(description="All-purpose interactive "
        "spectrum plotter for DBSP_DRP and intermediate PypeIt files.",
        formatter_class=argparse.RawTextHelpFormatter)

    argparser.add_argument("fname", type=str, help="path to FITS file")

    argparser.add_argument("--extension", type=str, default=None,
                           help="Extension name or number")

    argparser.add_argument("--BOX", action="store_false", dest="show_opt",
        help="Use boxcar extraction when plotting PypeIt spec1d files. "
        "By default the optimal extraction is plotted.")

    argparser.add_argument("--COUNTS", action="store_false", dest="show_flux",
        help="Plot counts when plotting PypeIt spec1d files. "
            "By default flux is plotted.")

    return argparser.parse_args() if options is None else argparser.parse_args(options)

# all purpose spectrum viewer
def main(args: argparse.Namespace) -> None:
    """
    Attempts to parse input FITS file as ``OneSpec``, ``SpecObjs``, and
    DBSP_DRP final data output.

    Args:
        args (argparse.Namespace): input arguments

    Raises:
        LookupError: Raised if the input extension is not found in the input
            FITS file.
        IndexError: Raised if the input integer extension is outside of the
            valid range for indexing the input FITS file's extensions.
    """
    try:
        extension = int(args.extension)
    except:
        extension = args.extension
    with fits.open(args.fname) as hdul:
        if isinstance(extension, str):
            exts = [hdu.name for hdu in hdul if hdu.name != "PRIMARY"]
            if extension.upper() not in exts:
                raise LookupError(f"Extension '{extension}' not found in "
                    f"{args.fname}.\n"
                    f"\tValid extensions present in {args.fname} are {exts}.")
        elif isinstance(extension, int):
            if (extension >= len(hdul) or extension <= -len(hdul) or extension == 0):
                raise IndexError(f"Extension index {extension} out of range: "
                            f"must be between 1 and {len(hdul) - 1} inclusive "
                            f"(or between -1 and {-len(hdul)+1} inclusive to "
                            "index from the end).")
    try:
        msgs.reset(verbosity=0)
        spectrum = OneSpec.from_file(args.fname)
        msgs.reset()
        basename = os.path.splitext(os.path.basename(args.fname))[0]
        plot(spectrum.wave, spectrum.flux, spectrum.ivar ** -0.5, basename)
    except PypeItError:
        try:
            spectrum = SpecObjs.from_fitsfile(args.fname)
            if extension is None:
                extension = 1
            msgs.reset()
            if isinstance(extension, int):
                if extension > 0:
                    extension -= 1
                sobj = spectrum[extension]
            elif isinstance(extension, str):
                sobj = spectrum[spectrum['NAME'] == extension][0]

            method = 'OPT' if args.show_opt else 'BOX'
            flux = 'FLAM' if args.show_flux else 'COUNTS'
            wave_name = f"{method}_WAVE"
            flux_name = f"{method}_{flux}"
            err_name = f"{method}_{flux}_SIG"
            name = f"{spectrum.header['TARGET']}[{sobj['NAME'].split('-')[0]}]"
            plot(sobj[wave_name], sobj[flux_name], sobj[err_name], name)
        except PypeItError:
            with fits.open(args.fname) as hdul:
                exts = [hdu.name for hdu in hdul if hdu.name != "PRIMARY"]
                if "RED" in exts and "BLUE" in exts and "SPLICED" in exts:
                    # it's one of our files
                    if extension is None:
                        extension = 'SPLICED'
                    spectrum = hdul[extension].data
                    basename = os.path.splitext(os.path.basename(args.fname))[0]
                    name = f"{basename}[{hdul[extension].name}]"
                    plot(spectrum['wave'], spectrum['flux'], spectrum['sigma'], name)
                else:
                    sys.exit("Could not automatically determine the spectrum format.")

def plot(wave: np.ndarray, flux: np.ndarray, err: np.ndarray, title: str) -> None:
    """
    Plots spectrum and error with sensible y-scale limits.

    Flux and error have units of :math:`10^{-17}\\mathrm{erg}/\\mathrm{s}/\\mathrm{cm}^2/\\overset{\\circ}{\\mathrm{A}}`.

    Args:
        wave (np.ndarray): Wavelengh array
        flux (np.ndarray): Flux array
        err (np.ndarray): Flux error array
        title (str): Plot title.
    """
    plt.step(wave, flux, c='k', label='spectrum')
    plt.step(wave, err, c='gray', label='error')

    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"Flux ($10^{-17}\mathrm{erg}/\mathrm{s}/\mathrm{cm}^2/\AA$)")

    plt.ylim(sensible_ylim(wave, flux))
    plt.legend()
    plt.title(title)
    plt.show()

def sensible_ylim(wave: np.ndarray, flux: np.ndarray) -> Tuple[float, float]:
    """
    Returns a tuple of sensible y-axis limits for plotting spectra.

    Usage:

    >>> plt.step(wave, flux)
    >>> plt.ylim(sensible_ylim(wave, flux))
    >>> plt.show()

    Args:
        wave (np.ndarray): wavelength array
        flux (np.ndarray): flux array

    Returns:
        Tuple[float, float]: bottom, top
    """
    top1 = np.abs(np.percentile(flux, 95)) * 1.5

    minwave = wave[int(len(wave) * 0.1)]
    maxwave = wave[int(len(wave) * 0.9)]

    top2 = np.max(flux[(wave > minwave) & (wave < maxwave)]) * 1.1

    top = max(top1, top2)
    bottom = -0.05 * top
    return bottom, top
