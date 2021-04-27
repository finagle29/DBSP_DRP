"""
Automated fluxing for P200 DBSP.
"""

import os
from typing import List, Dict

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

def make_sensfunc(standard_file: str, output_path: str, spectrograph: str,
        user_config_lines: List[str], debug: bool = False) -> str:
    """
    Makes a sensitivity function
    """
    try:
        spec1dfile = os.path.join(output_path, 'Science', standard_file)
        par = load_spectrograph(spectrograph).default_pypeit_par()
        default_cfg_lines = par.to_config()
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines = default_cfg_lines, merge_with=user_config_lines)

        outfile = os.path.join(output_path, standard_file.replace('spec1d', 'sens'))

        #par['sensfunc']['UVIS']['extinct_correct'] = False

        sensobj = sensfunc.SensFunc.get_instance(spec1dfile, outfile, par=par['sensfunc'], debug=debug)

        # read in spec1dfile to get wavelengths
        sobjs = SpecObjs.from_fitsfile(spec1dfile)
        wave_star = sobjs[0].OPT_WAVE
        orig_mask = sensobj.counts_mask.copy()
        mask = np.ones_like(orig_mask).astype(bool)

        # Mask out telluric lines for sensfunc generation
        if par['sensfunc']['algorithm'] == 'UVIS':
            tell_opt = np.any([((wave_star >= 6270.00) & (wave_star <= 6290.00)), # H2O
                            ((wave_star >= 6850.00) & (wave_star <= 6960.00)), # O2 telluric band
                            ((wave_star >= 7580.00) & (wave_star <= 7750.00)), # O2 telluric band
                            ((wave_star >= 7160.00) & (wave_star <= 7340.00)), # H2O
                            ((wave_star >= 8150.00) & (wave_star <= 8250.00))], axis=0) # H2O
            mask[tell_opt] = False

        # Mask blueward of approx. atmospheric cutoff
        mask[wave_star < 3000] = False

        sensobj.counts_mask &= mask

        sensobj.run()

        sensobj.out_table['SENS_ZEROPOINT_GPM'] = orig_mask.T
        sensobj.out_table['SENS_ZEROPOINT_FIT_GPM'] = orig_mask.T
        sensobj.save()
        return os.path.basename(outfile)
    except pypmsgs.PypeItError as err:
        print(f"ERROR creating sensitivity function using {standard_file}")
        print("Changing its frametype to science")
        print(str(err))
        return ""

def build_fluxfile(spec1d_to_sensfunc: Dict[str,str], archived_sensfunc: str,
        output_path: str, spectrograph: str, user_config_lines: List[str]) -> str:
    """
    Writes the fluxfile for fluxing.

    Args:
        spec1d_to_sensfunc (Dict[str,str]): maps spec1d filenames to the
            sensitivity function they should use
        output_path (str): reduction output path
        spectrograph (str): spectrograph name
        user_config_lines (List[str]): list of user-supplied PypeIt
            configuration lines

    Returns:
        str: path to created fluxfile
    """
    cfg_lines = user_config_lines[:]
    cfg_lines.append("\n")

    # data section
    cfg_lines.append('flux read')
    for spec1d, sensfun in spec1d_to_sensfunc.items():
        spec_path = os.path.join(output_path, 'Science', spec1d)
        if sensfun:
            sens_path = os.path.join(output_path, sensfun)
        else:
            sens_path = archived_sensfunc
        cfg_lines.append(f'  {spec_path} {sens_path}')
    cfg_lines.append('flux end')

    ofile = os.path.join(output_path, f'{spectrograph}.flux')

    with open(ofile, mode='wt') as f:
        for line in cfg_lines:
            f.write(line)
            f.write("\n")

    return ofile

def flux(flux_file: str, output_path: str, debug: bool = False) -> None:
    """
    Fluxes spectra.
    """
    # Load the file
    config_lines, spec1dfiles, sensfiles = read_fluxfile(flux_file)
    # Read in spectrograph from spec1dfile header
    header = fits.getheader(spec1dfiles[0])
    spectrograph = load_spectrograph(header['PYP_SPEC'])

    # Parameters
    spectrograph_def_par = spectrograph.default_pypeit_par()
    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(), merge_with=config_lines)
    # Write the par to disk
    par_outfile = os.path.join(output_path, f"{os.path.basename(flux_file)}.par")
    print(f"Writing the parameters to {par_outfile}")
    par.to_config(par_outfile)

    # Instantiate
    FxCalib = fluxcalibrate.FluxCalibrate.get_instance(spec1dfiles, sensfiles, par=par['fluxcalib'], debug=debug)
    msgs.info('Flux calibration complete')
