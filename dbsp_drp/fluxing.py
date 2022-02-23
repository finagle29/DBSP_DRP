"""
Automated fluxing for P200 DBSP.
"""

import os
from typing import List, Dict

from pkg_resources import resource_filename

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

archived_sensfuncs = [
    'blue_300_3990_d55', # current wavelength range 3000 - 7500 Å
    'blue_600_4000_d55', # current wavelength range 3000 - 6200 Å
    'blue_600_4000_d68', # current wavelength range 3730 - 6800 Å
    'red_316_7500_d55',  # current wavelength range 4600 - 11000 Å
    'red_600_10000_d55', # current wavelength range 5500 - 9000 Å
    'red_1200_7100_d68', # current wavelength range 7420 - 9060 Å
    'red_1200_9400_d55', # current wavelength range 8000 - 9600 Å
]

def make_sensfunc(standard_file: str, output_path: str, spectrograph: str,
        user_config_lines: List[str], debug: bool = False) -> str:
    """
    Makes a sensitivity function.

    Args:
        standard_file (str): Filename of standard exposure.
        output_path (str): Partial path to standard exposure.
        spectrograph (str): PypeIt name of spectrograph.
        user_config_lines (List[str]): User-provided PypeIt configuration.
        debug (bool, optional): Show debugging output/plots? Defaults to False.

    Returns:
        str: Filename of sensitivity function file, or empty string on failure.
    """
    try:
        spec1dfile = os.path.join(output_path, 'Science', standard_file)
        par = load_spectrograph(spectrograph).default_pypeit_par()
        default_cfg_lines = par.to_config()
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines = default_cfg_lines, merge_with=user_config_lines)

        outfile = os.path.join(output_path, standard_file.replace('spec1d', 'sens'))


        # read in spec1dfile to get wavelengths
        sobjs = SpecObjs.from_fitsfile(spec1dfile)
        wave_star = sobjs[0].OPT_WAVE
        orig_mask = sobjs[0].OPT_MASK.reshape((len(wave_star), 1))
        mask = np.ones_like(orig_mask).astype(bool)

        if 'red' in spectrograph:
            tell_opt = np.any([((wave_star >= 6270.00) & (wave_star <= 6290.00)), # H2O
                            ((wave_star >= 6850.00) & (wave_star <= 6960.00)), # O2 telluric band
                            ((wave_star >= 7580.00) & (wave_star <= 7750.00)), # O2 telluric band
                            ((wave_star >= 7160.00) & (wave_star <= 7340.00)), # H2O
                            ((wave_star >= 8150.00) & (wave_star <= 8250.00))], axis=0) # H2O
            mask[tell_opt] = False
            par['sensfunc']['UVIS']['nresln'] = 7.5
            par['sensfunc']['UVIS']['balm_mask_wid'] = 0
        elif 'blue' in spectrograph:
            mask[wave_star < 3000] = False
            par['sensfunc']['UVIS']['nresln'] = 4.5
            par['sensfunc']['UVIS']['balm_mask_wid'] = 0

        sensobj = sensfunc.SensFunc.get_instance(spec1dfile, outfile, par=par['sensfunc'], debug=debug)

        sensobj.counts_mask &= mask

        sensobj.run()

        #sensobj.out_table['SENS_ZEROPOINT_GPM'] = orig_mask.T
        #sensobj.out_table['SENS_ZEROPOINT_FIT_GPM'] = orig_mask.T
        sensobj.to_file(outfile, overwrite=True)
        return os.path.basename(outfile)
    except (pypmsgs.PypeItError, ValueError) as err:
        print(f"ERROR creating sensitivity function using {standard_file}")
        if isinstance(err, ValueError):
            print("This standard likely has insufficient wavelength coverage in the reference spectrum.")
            print("Next time, please don't use it, or suggest a better reference spectrum with better wavelength coverage.")
        print("Changing its frametype to science")
        print(str(err))
        return ""

def build_fluxfile(spec1d_to_sensfunc: Dict[str,str], output_path: str,
        spectrograph: str, user_config_lines: List[str]) -> str:
    """
    Writes the fluxfile for fluxing.

    Uses archived sensitivity function if no standard was reduced.

    Args:
        spec1d_to_sensfunc (Dict[str,str]): maps spec1d filenames to the
            sensitivity function they should use
        output_path (str): reduction output path
        spectrograph (str): PypeIt name of spectrograph.
        user_config_lines (List[str]): User-provided PypeIt configuration.

    Returns:
        str: path to created fluxfile
    """
    ## TODO: Refactor this into flux() below.
    cfg_lines = user_config_lines[:]
    # Minor kludge to deal with PypeIt#1230. Remove after PypeIt >= v1.5.0 is required.
    if (not any('extinct_correct' in line for line in cfg_lines) and
        not any(('algorithm' in line) and ('IR' in line) for line in cfg_lines)):
        cfg_lines.append('[fluxcalib]\n')
        cfg_lines.append('extinct_correct=True\n')
    cfg_lines.append("\n")

    # data section
    cfg_lines.append('flux read')
    for spec1d, sensfun in spec1d_to_sensfunc.items():
        spec_path = os.path.join(output_path, 'Science', spec1d)
        sens_path = os.path.join(output_path, sensfun)
        if not os.path.isfile(sens_path):
            if sensfun in archived_sensfuncs:
                sens_path = resource_filename("dbsp_drp", f"data/sens_{sensfun}.fits")
            else:
                print(f"WARNING: sensitivity function {sensfun} could not be located, and is not archived by DBSP_DRP")
                print("Please add a standard star exposure taken with the same configuration to your raw data folder "
                    "and restart data reduction")
                sens_path = ''
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
    Flux spectra.

    Args:
        flux_file (str): Path to flux file
        output_path (str): reduction output path
        debug (bool, optional): Show debugging output/plots? Defaults to False.
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
