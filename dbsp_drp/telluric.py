"""
Telluric correction for P200 DBSP.
"""
import os

import numpy as np
from astropy.io import fits

from pypeit.par import pypeitpar
from pypeit.core import telluric
from pypeit.spectrographs.util import load_spectrograph

def telluric_correct(args: dict):
    """
    method to telluric correct one coadded file
    """
    spec1dfile = os.path.join(args['output_path'], 'Science', args['spec1dfile'])
    spectrograph = load_spectrograph(args['spectrograph'])
    par = spectrograph.default_pypeit_par()

    par['tellfit']['objmodel'] = 'poly'
    par['tellfit']['fit_wv_min_max'] = [5500, 11000]
    par['tellfit']['model'] = 'exp' # maybe somehow choose between poly and exp??????? look at median
    par['tellfit']['polyorder'] = 8

    if par['tellfit']['tell_grid'] is None:
        if par['sensfunc']['IR']['telgridfile'] is not None:
            par['tellfit']['tell_grid'] = par['sensfunc']['IR']['telgridfile']

    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=par.to_config(), merge_with=args['user_config_lines'])

    # Parse the output filename
    outfile = os.path.join(args['output_path'], 'Science', args['spec1dfile'].replace('.fits','_tellcorr.fits'))
    modelfile = os.path.join(args['output_path'], 'Science', args['spec1dfile'].replace('.fits','_tellmodel.fits'))


    try:
        TelPoly = telluric.poly_telluric(spec1dfile, par['tellfit']['tell_grid'], modelfile, outfile,
                                         z_obj=par['tellfit']['redshift'],
                                         func=par['tellfit']['func'], model=par['tellfit']['model'],
                                         polyorder=par['tellfit']['polyorder'],
                                         fit_wv_min_max=par['tellfit']['fit_wv_min_max'],
                                         mask_lyman_a=par['tellfit']['mask_lyman_a'],
                                         delta_coeff_bounds=par['tellfit']['delta_coeff_bounds'],
                                         minmax_coeff_bounds=par['tellfit']['minmax_coeff_bounds'],
                                         only_orders=par['tellfit']['only_orders'],
                                         debug_init=args['debug'], disp=args['debug'], debug=args['debug'], show=args['plot'])
    except ValueError:
        print(f"ERROR!! Telluric correction of {args['spec1dfile']} FAILED!")
