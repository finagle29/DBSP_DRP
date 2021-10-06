"""
Telluric correction for P200 DBSP.
"""
import os
from typing import List

from pypeit.par import pypeitpar
from pypeit.core import telluric
from pypeit.spectrographs.util import load_spectrograph

def telluric_correct(coadd: str, output_path: str, spectrograph: str,
        user_config_lines: List[str], debug: bool = False, plot: bool = False):
    """
    Telluric correct one coadd file.

    Args:
        coadd (str): Coadd filename.
        output_path (str): reduction output path
        spectrograph (str): PypeIt name of spectrograph.
        user_config_lines (List[str]): User-provided PypeIt configuration
        debug (bool, optional): Show debugging output? Defaults to False.
        plot (bool, optional): Show debugging plots? Defaults to False.
    """
    coadd_path = os.path.join(output_path, 'Science', coadd)
    spectrograph = load_spectrograph(spectrograph)
    par = spectrograph.default_pypeit_par()

    par['telluric']['objmodel'] = 'poly'
    ## TODO: Change fit_wv_min_max based on where the red detector defect is.
    par['telluric']['fit_wv_min_max'] = [5500, 11000]
    # maybe somehow choose between poly and exp??????? look at median
    par['telluric']['model'] = 'exp'
    par['telluric']['polyorder'] = 8

    if par['telluric']['telgridfile'] is None:
        if par['sensfunc']['IR']['telgridfile'] is not None:
            par['telluric']['telgridfile'] = par['sensfunc']['IR']['telgridfile']

    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=par.to_config(),
        merge_with=user_config_lines)

    # Parse the output filename
    outfile = os.path.join(output_path, 'Science', os.path.splitext(coadd)[0] + '_tellcorr.fits')
    modelfile = os.path.join(output_path, 'Science', os.path.splitext(coadd)[0] + '_tellmodel.fits')


    try:
        TelPoly = telluric.poly_telluric(coadd_path, par['telluric']['telgridfile'],
            modelfile, outfile,
            z_obj=par['telluric']['redshift'],
            func=par['telluric']['func'],
            model=par['telluric']['model'],
            polyorder=par['telluric']['polyorder'],
            fit_wv_min_max=par['telluric']['fit_wv_min_max'],
            mask_lyman_a=par['telluric']['mask_lyman_a'],
            delta_coeff_bounds=par['telluric']['delta_coeff_bounds'],
            minmax_coeff_bounds=par['telluric']['minmax_coeff_bounds'],
            only_orders=par['telluric']['only_orders'],
            debug_init=debug, disp=debug, debug=debug, show=plot)
    except ValueError:
        print(f"ERROR!! Telluric correction of {coadd} FAILED!")

def picklable_telluric_correct(args):
    telluric_correct(*args)
