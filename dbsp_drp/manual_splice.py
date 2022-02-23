import argparse
import os
from typing import Optional, List
import re

from astropy.io import fits

# needed for version info
import numpy as np
import astropy
import pypeit
import dbsp_drp

from pypeit.history import History
from dbsp_drp import splicing, adjust_splicing

def entrypoint():
    main(parse())

def parse(options: Optional[List[str]] = None) -> argparse.Namespace:
    argparser = argparse.ArgumentParser(description="Manually splice two coadded files together for DBSP.\n"
        "After preparing the spliced spectrum, the dbsp_adjust_splicing GUI pops up to allow for the splicing "
        "to be adjusted. To save the manually spliced spectrum, you MUST press save in the GUI.",
        formatter_class=argparse.RawTextHelpFormatter)

    argparser.add_argument("raw_data_path", type=str, help="Path to raw data from this reduction")

    argparser.add_argument("outfile", type=str, help="Destination of spliced spectrum.")

    argparser.add_argument("-r", "--red_file", type=str, default=None,
        help="redNNNN-redMMMM_SPATXXXX.fits file.")

    argparser.add_argument("-b", "--blue_file", type=str, default=None,
        help="blueNNNN-redMMMM_SPATXXXX.fits file.")

    return argparser.parse_args() if options is None else argparser.parse_args(options)

def find_spec1ds_spats(history: History):
    unique_spec1ds = {}
    spec1ds_spats = []

    for line in history.history:
        match = re.search(r'File (\d+) "(.*)"', line)
        if match:
            unique_spec1ds[match.group(1)] = match.group(2)

    for line in history.history:
        match = re.search(r"Object ID (.*) from file (\d+)", line)
        if match:
            objid = match.group(1)
            spec1d = unique_spec1ds[match.group(2)]
            spat = int(re.search(r"SPAT(\d+)", objid).group(1))
            spec1ds_spats.append((spec1d, spat))

    return spec1ds_spats

def main(args: argparse.Namespace):
    if args.red_file is None and args.blue_file is None:
        raise ValueError("You must provide at least one of red_file and blue_file")

    # infer output dir
    output_path = os.path.abspath(
        os.path.join(
            os.path.dirname(args.red_file if args.red_file is not None else args.blue_file),
            "..")
        )

    # read them in
    blue_history = History()
    red_history = History()

    spec_b = None
    spec_r = None

    if args.blue_file is not None:
        hdul_b = fits.open(args.blue_file)
        spec_b = hdul_b[1].data
        blue_history = History(hdul_b[0].header)

    if args.red_file is not None:
        hdul_r = fits.open(args.red_file)
        spec_r = hdul_r[1].data
        red_history = History(hdul_r[0].header)

    # find the spec1ds and spats
    blue_spec1ds_spats = find_spec1ds_spats(blue_history)
    red_spec1ds_spats = find_spec1ds_spats(red_history)

    # initial aco
    ((final_wvs, final_flam, final_flam_sig),
        (red_wvs, red_flam, red_sig),
        (blue_wvs, blue_flam, blue_sig)) = splicing.adjust_and_combine_overlap(spec_b, spec_r, False)

    # build the final data product
    primary_header = fits.Header()
    primary_header['HIERARCH DBSP_DRP_V'] = dbsp_drp.__version__
    primary_header['PYPEIT_V'] = pypeit.__version__
    primary_header['NUMPY_V'] = np.__version__
    primary_header['HIERARCH ASTROPY_V'] = astropy.__version__
    primary_header['B_COADD'] = args.blue_file
    primary_header['R_COADD'] = args.red_file
    primary_hdu = fits.PrimaryHDU(header=primary_header)

    raw_red_hdus = splicing.get_raw_hdus_from_spec1d(red_spec1ds_spats, args.raw_data_path, output_path)
    raw_blue_hdus = splicing.get_raw_hdus_from_spec1d(blue_spec1ds_spats, args.raw_data_path, output_path)

    col_wvs = fits.Column(name='wave', array=red_wvs, unit='ANGSTROM', format='D')
    col_flux = fits.Column(name='flux', array=red_flam, unit='E-17 ERG/S/CM^2/ANG', format='D')
    col_error = fits.Column(name='sigma', array=red_sig, unit='E-17 ERG/S/CM^2/ANG', format='D')
    red_hdu = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error], name="RED")

    col_wvs = fits.Column(name='wave', array=blue_wvs, unit='ANGSTROM', format='D')
    col_flux = fits.Column(name='flux', array=blue_flam, unit='E-17 ERG/S/CM^2/ANG', format='D')
    col_error = fits.Column(name='sigma', array=blue_sig, unit='E-17 ERG/S/CM^2/ANG', format='D')
    blue_hdu = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error], name="BLUE")

    col_wvs = fits.Column(name='wave', array=final_wvs, unit='ANGSTROM', format='D')
    col_flux = fits.Column(name='flux', array=final_flam, unit='E-17 ERG/S/CM^2/ANG', format='D')
    col_error = fits.Column(name='sigma', array=final_flam_sig, unit='E-17 ERG/S/CM^2/ANG', format='D')
    table_hdu = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error], name="SPLICED")

    hdul = fits.HDUList(hdus=[primary_hdu, *raw_red_hdus, *raw_blue_hdus, red_hdu, blue_hdu, table_hdu])

    # before we save, run adjust splicing GUI
    adjust_splicing.adjust_splicing_GUI(hdul, args.outfile)
