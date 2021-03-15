"""
Automated fluxing for P200 DBSP.
"""

import os

import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit import pypmsgs
from pypeit.par import pypeitpar
from pypeit.specobjs import SpecObjs
from pypeit import sensfunc
from pypeit import fluxcalibrate
from pypeit.scripts.flux_calib import read_fluxfile
from pypeit.spectrographs.util import load_spectrograph

def make_sensfunc(args: dict) -> str:
    """
    Makes a sensitivity function
    """
    try:
        spec1dfile = os.path.join(args['output_path'], 'Science', args['spec1dfile'])
        par = load_spectrograph(args['spectrograph']).default_pypeit_par()
        default_cfg_lines = par.to_config()
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines = default_cfg_lines, merge_with=args['user_config_lines'])

        outfile = os.path.join(args['output_path'], args['spec1dfile'].replace('spec1d', 'sens'))

        #par['sensfunc']['UVIS']['extinct_correct'] = False

        sensobj = sensfunc.SensFunc.get_instance(spec1dfile, outfile, par=par['sensfunc'], debug=args['debug'])

        if 'red' in args['spectrograph']:
            # read in spec1dfile to get wavelengths
            sobjs = SpecObjs.from_fitsfile(spec1dfile)
            wave_star = sobjs[0].OPT_WAVE
            orig_mask = sensobj.counts_mask.copy()

            mask_tell = np.ones_like(wave_star).astype(bool)
            tell_opt = np.any([((wave_star >= 6270.00) & (wave_star <= 6290.00)), # H2O
                            ((wave_star >= 6850.00) & (wave_star <= 6960.00)), # O2 telluric band
                            ((wave_star >= 7580.00) & (wave_star <= 7750.00)), # O2 telluric band
                            ((wave_star >= 7160.00) & (wave_star <= 7340.00)), # H2O
                            ((wave_star >= 8150.00) & (wave_star <= 8250.00))], axis=0) # H2O
            mask_tell[tell_opt] = False

            sensobj.counts_mask &= mask_tell


        sensobj.run()
        if 'red' in args['spectrograph']:
            sensobj.out_table['MASK_SENS'] = orig_mask

        sensobj.save()
        return os.path.basename(outfile)
    except pypmsgs.PypeItError as err:
        print(f"ERROR creating sensitivity function using {args['spec1dfile']}")
        print("Changing its frametype to science")
        print(str(err))
        return ""

def build_fluxfile(args: dict) -> str:
    """
    Writes the fluxfile for fluxing.
    """
    # args['spec1dfiles'] is a dict mapping spec1d files to the sensitivity function file they should use
    cfg_lines = args['user_config_lines'][:]
    cfg_lines.append("\n")

    # data section
    cfg_lines.append('flux read')
    for spec1d, sensfun in args['spec1dfiles'].items():
        spec_path = os.path.join(args['output_path'], 'Science', spec1d)
        sens_path = os.path.join(args['output_path'], sensfun)
        cfg_lines.append(f'  {spec_path} {sens_path}')
    cfg_lines.append('flux end')

    ofile = os.path.join(args['output_path'], f'{args["spectrograph"]}.flux')

    with open(ofile, mode='wt') as f:
        for line in cfg_lines:
            f.write(line)
            f.write("\n")

    return ofile

def flux(args: dict) -> None:
    """
    Fluxes spectra.
    """
    # Load the file
    config_lines, spec1dfiles, sensfiles = read_fluxfile(args['flux_file'])
    # Read in spectrograph from spec1dfile header
    header = fits.getheader(spec1dfiles[0])
    spectrograph = load_spectrograph(header['PYP_SPEC'])

    # Parameters
    spectrograph_def_par = spectrograph.default_pypeit_par()
    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(), merge_with=config_lines)
    # Write the par to disk
    par_outfile = os.path.join(args['output_path'], f"{os.path.basename(args['flux_file'])}.par")
    print(f"Writing the parameters to {par_outfile}")
    par.to_config(par_outfile)

    # Instantiate
    FxCalib = fluxcalibrate.FluxCalibrate.get_instance(spec1dfiles, sensfiles, par=par['fluxcalib'], debug=args['debug'])
    msgs.info('Flux calibration complete')
