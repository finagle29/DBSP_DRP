import argparse
import os
from typing import Optional, List

import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox, Button, CheckButtons

from astropy.io import fits

from PySide2 import QtWidgets

from dbsp_drp.splicing import adjust_and_combine_overlap
from dbsp_drp import show_spectrum

def entrypoint():
    main(parse())

def parse(options: Optional[List[str]] = None) -> argparse.Namespace:
    argparser = argparse.ArgumentParser(description="Red/Blue splicing adjustment for P200 DBSP",
        formatter_class=argparse.RawTextHelpFormatter)

    argparser.add_argument("fname", type=str, nargs="+", help="Final data "
        "output from DBSP_DRP to adjust the splicing of. Can be one or "
        "multiple files.")

    return argparser.parse_args() if options is None else argparser.parse_args(options)

def main(args: argparse.Namespace):
    for fname in args.fname:
        if not os.path.isfile(fname):
            raise FileNotFoundError(f"File {fname} not found!")

        hdul = fits.open(fname, mode="update")
        if hdul[0].header['R_COADD'] and hdul[0].header['B_COADD']:
            adjust_splicing_GUI(hdul, fname)
        else:
            missing_side = 'blue' if hdul[0].header['R_COADD'] else 'red'
            print(f"Cannot adjust splicing of {fname} because {missing_side} data is not present.")

def adjust_splicing_GUI(hdul: fits.HDUList, fname: str):
    """
    Opens a Matplotlib GUI for users to manually adjust the matching of overall
    red/blue flux levels, as well as interpolate over detector gaps.

    Args:
        hdul (fits.HDUList): HDUList read in from DBSP_DRP final data output file.
        fname (str): Output path and filename to write adjusted spectrum to.
    """
    spec_r = hdul['RED'].data
    spec_b = hdul['BLUE'].data
    spliced = hdul['SPLICED'].data

    red_mult = hdul['SPLICED'].header.get('RED_MULT', 1.0)
    gaps = hdul['SPLICED'].header.get('INTERP_GAPS', False)

    # prepare for aco
    spec_r = fits.FITS_rec.from_columns([*spec_r.columns, fits.Column(name='ivar', array=spec_r['sigma'] ** -2.0, format='D')])
    spec_b = fits.FITS_rec.from_columns([*spec_b.columns, fits.Column(name='ivar', array=spec_b['sigma'] ** -2.0, format='D')])

    #fig = plt.figure(figsize=(10,5))
    fig, ax = plt.subplots(figsize=(10,5))
    fig.subplots_adjust(bottom=0.2, right=0.8)
    ax.plot(spec_b['wave'], spec_b['flux'], c='b', alpha=0.5)
    red_l, = ax.plot(spec_r['wave'], spec_r['flux'], c='r', alpha=0.5)
    spliced_l, = ax.plot(spliced['wave'], spliced['flux'], c='k')

    bottom_b, top_b = show_spectrum.sensible_ylim(spec_b['wave'], spec_b['flux'])
    bottom_r, top_r = show_spectrum.sensible_ylim(spec_r['wave'], spec_r['flux'])
    bottom_s, top_s = show_spectrum.sensible_ylim(spliced['wave'], spliced['flux'])

    bottom = min(bottom_b, bottom_r, bottom_s)
    top = max(top_b, top_r, top_s)
    ax.set_ylim(bottom, top)

    ax.set_xlabel(r"Wavelength ($\AA$)")
    ax.set_ylabel(r"Flux ($10^{-17}\mathrm{erg}/\mathrm{s}/\mathrm{cm}^2/\AA$)")

    def update(_):
        try:
            red_mult = float(text_box.text)
            interp_gaps = check.get_status()[0]

            wave, flux, _ = adjust_and_combine_overlap(spec_b, spec_r, interp_gaps, red_mult=red_mult)[0]

            red_l.set_ydata(spec_r['flux'] * red_mult)
            spliced_l.set_data(wave, flux)

            plt.draw()
        except ValueError:
            pass

    def save(_):
        try:
            red_mult = float(text_box.text)
            interp_gaps = check.get_status()[0]
            wave, flux, sig = adjust_and_combine_overlap(spec_b, spec_r, interp_gaps, red_mult=red_mult)[0]

            col_wvs = fits.Column(name='wave', array=wave, unit='ANGSTROM', format='D')
            col_flux = fits.Column(name='flux', array=flux, unit='E-17 ERG/S/CM^2/ANG', format='D')
            col_error = fits.Column(name='sigma', array=sig, unit='E-17 ERG/S/CM^2/ANG', format='D')
            table_hdu = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error], name="SPLICED")

            table_hdu.header['RED_MULT'] = red_mult
            table_hdu.header['HIERARCH INTERP_GAPS'] = interp_gaps

            hdul.pop()
            hdul.append(table_hdu)
            hdul.writeto(fname, overwrite=True)
        except Exception as e:
            QtWidgets.QMessageBox.warning(None, "Error occured while saving changes: " + e)

    axbox = fig.add_axes([0.15, 0.05, 0.4, 0.05])
    text_box = TextBox(axbox, "Red multiplier", str(red_mult))
    text_box.on_submit(update)

    axbutton = fig.add_axes([0.6, 0.05, 0.1, 0.05])
    button = Button(axbutton, "Save")
    button.on_clicked(save)

    axcheck = fig.add_axes([0.81, 0.35, 0.15, 0.3])
    axcheck.set_frame_on(False)
    check = CheckButtons(axcheck, ['interpolate gaps'], [gaps])
    check.on_clicked(update)

    plt.show()
