"""
Telluric correction for P200 DBSP.
"""
import os
from typing import List

import numpy as np
from astropy.io import fits

from pypeit.par import pypeitpar
from pypeit.core import telluric
from pypeit.spectrographs.util import load_spectrograph

def telluric_correct(coadd: str, output_path: str, spectrograph: str,
        user_config_lines: List[str], debug: bool = False, plot: bool = False):
    """
    method to telluric correct one coadded file
    """
    ## TODO:
    # factor this to be instrument-independent
    # perhaps give Instrument a method to mutate par?
    # or give it lines
    coadd_path = os.path.join(output_path, 'Science', coadd)
    spectrograph = load_spectrograph(spectrograph)
    par = spectrograph.default_pypeit_par()

    par['tellfit']['objmodel'] = 'poly'
    par['tellfit']['fit_wv_min_max'] = [5500, 11000]
    # maybe somehow choose between poly and exp??????? look at median
    par['tellfit']['model'] = 'exp'
    par['tellfit']['polyorder'] = 8

    if par['tellfit']['tell_grid'] is None:
        if par['sensfunc']['IR']['telgridfile'] is not None:
            par['tellfit']['tell_grid'] = par['sensfunc']['IR']['telgridfile']

    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=par.to_config(),
        merge_with=user_config_lines)

    # Parse the output filename
    outfile = os.path.splitext(coadd)[0] + '_tellcorr.fits'
    modelfile = os.path.splitext(coadd)[0] + '_tellmodel.fits'

    try:
        TelPoly = telluric.poly_telluric(coadd_path, par['tellfit']['tell_grid'],
            modelfile, outfile,
            z_obj=par['tellfit']['redshift'],
            func=par['tellfit']['func'],
            model=par['tellfit']['model'],
            polyorder=par['tellfit']['polyorder'],
            fit_wv_min_max=par['tellfit']['fit_wv_min_max'],
            mask_lyman_a=par['tellfit']['mask_lyman_a'],
            delta_coeff_bounds=par['tellfit']['delta_coeff_bounds'],
            minmax_coeff_bounds=par['tellfit']['minmax_coeff_bounds'],
            only_orders=par['tellfit']['only_orders'],
            debug_init=debug, disp=debug, debug=debug, show=plot)
    except ValueError:
        print(f"ERROR!! Telluric correction of {coadd} FAILED!")

def picklable_telluric_correct(args):
    telluric_correct(*args)
