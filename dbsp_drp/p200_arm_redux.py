"""
Automated Reduction Pipeline for one arm of P200 DBSP.
"""

import os
import glob
import shutil
import fileinput

from pkg_resources import resource_filename
from typing import Tuple, List

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.stats
import astropy.table
from astropy.visualization import ZScaleInterval, ImageNormalize

from configobj import ConfigObj
from yattag import Doc, indent

from pypeit.pypeitsetup import PypeItSetup
import pypeit
import pypeit.pypeit
from pypeit import msgs
from pypeit import pypmsgs
from pypeit import sensfunc
from pypeit import coadd1d
from pypeit.core import telluric
from pypeit.spec2dobj import Spec2DObj
from pypeit.specobjs import SpecObjs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts.flux_calib import read_fluxfile
from pypeit import fluxcalibrate
from pypeit.par import pypeitpar

import dbsp_drp
from dbsp_drp.manual_tracing import ManualTracingGUI



def parse_pypeit_parameter_file(args: dict) -> None:
    user_config_lines = []
    read_lines = False
    arm = 'red' if 'red' in args['spectrograph'] else 'blue'
    other_arm = 'red' if 'blue' in arm else 'blue'

    if os.path.isfile(args['parameter_file']):
        with open(args['parameter_file']) as f:
            for line in f.readlines():
                if f'[{other_arm}]' in line:
                    read_lines = False
                if read_lines:
                    user_config_lines.append(line)
                if f'[{arm}]' in line:
                    read_lines = True
    args['user_config_lines'] = user_config_lines

def setup(args: dict) -> Tuple[PypeItSetup, str]:
    """Does PypeIt setup, without writing the .pypeit file

    Args:
        args (dict): [description]

    Raises:
        ValueError: [description]
        ValueError: [description]
        IOError: [description]

    Returns:
        Tuple[PypeItSetup, str]: [description]
    """

    # Get the output directory
    output_path = os.getcwd() if args['output_path'] is None else args['output_path']
    sort_dir = os.path.join(output_path, 'setup_files')

    # Initialize PypeItSetup based on the arguments
    if args['root'] is not None:
        ps = PypeItSetup.from_file_root(args['root'], args['spectrograph'],
                                        extension=args['extension'], output_path=sort_dir)
    else:
        # Should never reach here
        raise IOError('Need to set -r !!')

    # Run the setup
    ps.run(setup_only=True, sort_dir=sort_dir, write_bkg_pairs=args['background'])

    return (ps, output_path)

def write_setup(args: dict, context: Tuple[PypeItSetup, str]) -> List[str]:
    """
    Writes the .pypeit file
    """
    ps, output_path = context
    # Use PypeItMetaData to write the complete PypeIt file
    config_list = [item.strip() for item in args['cfg_split'].split(',')]

    ps.user_cfg.append('[calibrations]')
    ps.user_cfg.append('master_dir = Master_' + args['spectrograph'].split('_')[-1])

    user_configobj = ConfigObj(ps.user_cfg)
    user_configobj.merge(ConfigObj(args['user_config_lines']))
    ps.user_cfg = [line + "\n" for line in user_configobj.write()]

    return ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg,
                                   write_bkg_pairs=args['background'], configs=config_list)

def redux(args: dict) -> None:
    """
    Runs the reduction
    """
    splitnm = os.path.splitext(args['pypeit_file'])
    if splitnm[1] != '.pypeit':
        msgs.error("Bad extension for PypeIt reduction file."+msgs.newline()+".pypeit is required")
    logname = splitnm[0] + ".log"

    pypeIt = pypeit.pypeit.PypeIt(args['pypeit_file'], verbosity=2,
                           reuse_masters=not args['do_not_reuse_masters'],
                           overwrite=True, #args['overwrite'],
                           redux_path=args['output_path'], #args['redux_path'], # rethink these
                           calib_only=args['calib_only'],
                           logname=logname, show=args['show'])
    pypeIt.reduce_all()
    msgs.info('Data reduction complete')

    msgs.info('Generating QA HTML')
    pypeIt.build_qa()

    args['output_spec1ds'] |= set(filter(os.path.isfile, [
            os.path.abspath(pypeIt.spec_output_file(i)) \
            for i in range(len(pypeIt.fitstbl.table)) \
            if pypeIt.fitstbl.table[i]['frametype'] in ['science', 'standard']
        ]))

    args['output_spec2ds'] |= set(filter(os.path.isfile, [
            os.path.abspath(pypeIt.spec_output_file(i, True)) \
            for i in range(len(pypeIt.fitstbl.table)) \
            if pypeIt.fitstbl.table[i]['frametype'] in ['science', 'standard']
        ]))

def make_sensfunc(args: dict) -> str:
    """
    Makes a sensitivity function
    """
    try:
        par = load_spectrograph(args['spectrograph']).default_pypeit_par()
        default_cfg_lines = par.to_config()
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines = default_cfg_lines, merge_with=args['user_config_lines'])

        outfile = os.path.join(args['output_path'],
            os.path.basename(args['spec1dfile']).replace('spec1d', 'sens'))

        #par['sensfunc']['UVIS']['extinct_correct'] = False

        sensobj = sensfunc.SensFunc.get_instance(args['spec1dfile'], outfile, par=par['sensfunc'], debug=args['debug'])

        if 'red' in args['spectrograph']:
            # read in spec1dfile to get wavelengths
            sobjs = SpecObjs.from_fitsfile(args['spec1dfile'])
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
        return outfile
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
        cfg_lines.append(f'  {spec1d} {sensfun}')
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

def get_raw_header_from_coadd(coadd: str, args: dict) -> str:
    if coadd is not None:
        return fits.open(os.path.join(os.path.dirname(args['root']), f"{os.path.basename(coadd).split('-')[0]}.fits"))[0].header
    else:
        return fits.Header()

def get_raw_hdus_from_spec1d(spec1d_list: List[str], args: dict) -> List[fits.BinTableHDU]:
    ret = []
    for (spec1d, spat) in spec1d_list:
        raw_fname = os.path.join(os.path.dirname(args['root']), f"{os.path.basename(spec1d).split('_')[1].split('-')[0]}.fits")
        # get the raw header
        with fits.open(raw_fname) as raw_hdul:
            raw_header = raw_hdul[0].header.copy()
        # get the spectrum
        with fits.open(spec1d) as spec1d_hdul:
            for hdu in spec1d_hdul:
                if f'SPAT{spat:04d}' in hdu.name:
                    raw_data = hdu.data.copy()
        wave_col = fits.Column(name='wave', array=raw_data['OPT_WAVE'], unit='ANGSTROM', format='D')
        flux_col = fits.Column(name='flux', array=raw_data['OPT_FLAM'], unit='E-17 ERG/S/CM^2/ANG', format='D')
        sigma_col = fits.Column(name='sigma', array=raw_data['OPT_FLAM_SIG'], unit='E-17 ERG/S/CM^2/ANG', format='D')
        ret.append(fits.BinTableHDU.from_columns([wave_col, flux_col, sigma_col],
            name=os.path.splitext(os.path.basename(raw_fname))[0].upper(),
            header=raw_header))
    return ret


def splice(args: dict) -> None:
    """
    Splices red and blue spectra together.
    """
    for target, targets_dict in args['splicing_dict'].items():
        label = 'a'
        for _, arm_dict in targets_dict.items():
            red_dict = arm_dict.get('red', {})
            blue_dict = arm_dict.get('blue', {})

            bluefile = blue_dict.get('coadd')
            redfile = red_dict.get('coadd')
            if bluefile is None and redfile is None:
                continue

            ((final_wvs, final_flam, final_flam_sig),
                (red_wvs, red_flam, red_sig),
                (blue_wvs, blue_flam, blue_sig)) = adjust_and_combine_overlap(bluefile, redfile, target)

            primary_header = fits.Header()
            primary_header['DBSP_DRP_V'] = dbsp_drp.__version__
            primary_header['PYPEIT_V'] = pypeit.__version__
            primary_header['NUMPY_V'] = np.__version__
            primary_header['ASTROPY_V'] = astropy.__version__
            primary_hdu = fits.PrimaryHDU(header=primary_header)

            raw_red_hdus = get_raw_hdus_from_spec1d(red_dict.get('spec1ds', []), args)
            raw_blue_hdus = get_raw_hdus_from_spec1d(blue_dict.get('spec1ds', []), args)


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

            hdul.writeto(os.path.join(args['output_path'], "Science", f'{target}_{label}.fits'), overwrite=True)
            label = chr(ord(label) + 1)



def splice_old(args: dict) -> None:
    """
    Splices red and blue spectra together.
    """
    for target, targ_dict in args['splicing_dict'].items():
        bluefile = targ_dict.get('blue')
        redfile = targ_dict.get('red')

        if bluefile is None and redfile is None:
            continue

        ((final_wvs, final_flam, final_flam_sig),
            (red_wvs, red_flam, red_sig),
            (blue_wvs, blue_flam, blue_sig)) = adjust_and_combine_overlap(bluefile, redfile, target)

        primary_header = fits.Header()
        primary_header['DBSP_DRP_V'] = dbsp_drp.__version__
        primary_header['PYPEIT_V'] = pypeit.__version__
        primary_header['NUMPY_V'] = np.__version__
        primary_header['ASTROPY_V'] = astropy.__version__
        primary_hdu = fits.PrimaryHDU(header=primary_header)

        ## Here we need something a little more involved to get FITS headers in
        ## copy source red fits header
        red_header = get_raw_header_from_coadd(redfile, args)
        col_wvs = fits.Column(name='wave', array=red_wvs, unit='ANGSTROM', format='D')
        col_flux = fits.Column(name='flux', array=red_flam, unit='E-17 ERG/S/CM^2/ANG', format='D')
        col_error = fits.Column(name='sigma', array=red_sig, unit='E-17 ERG/S/CM^2/ANG', format='D')
        red_hdu = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error], header=red_header, name="RED")

        ## copy source blue fits header
        blue_header = get_raw_header_from_coadd(bluefile, args)
        col_wvs = fits.Column(name='wave', array=blue_wvs, unit='ANGSTROM', format='D')
        col_flux = fits.Column(name='flux', array=blue_flam, unit='E-17 ERG/S/CM^2/ANG', format='D')
        col_error = fits.Column(name='sigma', array=blue_sig, unit='E-17 ERG/S/CM^2/ANG', format='D')
        blue_hdu = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error], header=blue_header, name="BLUE")

        col_wvs = fits.Column(name='wave', array=final_wvs, unit='ANGSTROM', format='D')
        col_flux = fits.Column(name='flux', array=final_flam, unit='E-17 ERG/S/CM^2/ANG', format='D')
        col_error = fits.Column(name='sigma', array=final_flam_sig, unit='E-17 ERG/S/CM^2/ANG', format='D')
        table_hdu = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error], name="SPLICED")

        hdul = fits.HDUList(hdus=[primary_hdu, red_hdu, blue_hdu, table_hdu])

        hdul.writeto(os.path.join(args['output_path'], "Science", f'{target}.fits'), overwrite=True)

def adjust_and_combine_overlap(bluefile: str, redfile: str, target: str) -> Tuple[
            Tuple[np.ndarray, np.ndarray, np.ndarray],
            Tuple[np.ndarray, np.ndarray, np.ndarray],
            Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    # TODO: propagate input masks

    if bluefile is not None:
        spec_b = fits.open(bluefile)
        if redfile is None:
            return ((spec_b[1].data['wave'], spec_b[1].data['flux'], spec_b[1].data['ivar'] ** -0.5),
                (None, None, None),
                (spec_b[1].data['wave'], spec_b[1].data['flux'], spec_b[1].data['ivar'] ** -0.5))
    if redfile is not None:
        spec_r = fits.open(redfile)
        if bluefile is None:
            return ((spec_r[1].data['wave'], spec_r[1].data['flux'], spec_r[1].data['ivar'] ** -0.5),
                (spec_r[1].data['wave'], spec_r[1].data['flux'], spec_r[1].data['ivar'] ** -0.5),
                (None, None, None))

    # combination steps
    overlap_lo = spec_r[1].data['wave'][0]
    overlap_hi = spec_b[1].data['wave'][-1]
    # maybe need to fix the overlaps?
    # need more finely spaced grid to be completely contained within coarser grid

    lw = 0.5
    plt.figure(figsize=(20,10))

    mask_b = ~((spec_b[1].data['ivar'] < 0.1) & (spec_b[1].data['wave'] > overlap_lo))
    mask_r = ~((spec_r[1].data['ivar'] < 0.1) & (spec_r[1].data['wave'] < overlap_hi))

    olap_r = (spec_r[1].data['wave'] < overlap_hi)
    olap_b = (spec_b[1].data['wave'] > overlap_lo)

    red_mult = 1
    red_mult = (astropy.stats.sigma_clipped_stats(spec_b[1].data['flux'][olap_b])[1] /
        astropy.stats.sigma_clipped_stats(spec_r[1].data['flux'][olap_r])[1])
    #red_mult = np.average(spec_aag[1].data['OPT_FLAM'][olap_b], weights=spec_aag[1].data['OPT_FLAM_SIG'][olap_b] ** -2.0)\
    #   /np.average(spec_aag_red[1].data['OPT_FLAM'][olap_r], weights=spec_aag_red[1].data['OPT_FLAM_SIG'][olap_r] ** -2.0)
    if red_mult > 3 or 1/red_mult > 3:
        msgs.warn(f"For {target}, red spectrum is {red_mult} times less flux than blue spectrum in overlap region." +
            "The red and blue traces may not correspond to the same object.")


    # different dispersion.
    wvs_b = spec_b[1].data['wave'][~olap_b]
    wvs_r = spec_r[1].data['wave'][~olap_r]
    flam_b = spec_b[1].data['flux'][~olap_b]
    flam_r = spec_r[1].data['flux'][~olap_r]
    flam_sig_b = spec_b[1].data['ivar'][~olap_b] ** -0.5
    flam_sig_r = spec_r[1].data['ivar'][~olap_r] ** -0.5


    olap_wvs_r = spec_r[1].data['wave'][olap_r]
    olap_flam_r = red_mult * spec_r[1].data['flux'][olap_r]
    olap_flam_sig_r = red_mult * spec_r[1].data['ivar'][olap_r] ** -0.5
    olap_wvs_b = spec_b[1].data['wave'][olap_b][:-1]
    olap_flam_b = spec_b[1].data['flux'][olap_b][:-1]
    olap_flam_sig_b = spec_b[1].data['ivar'][olap_b][:-1] ** -0.5

    olap_flam_r_interp = np.interp(olap_wvs_b, olap_wvs_r, olap_flam_r)
    olap_flam_r_interp, olap_flam_sig_r_interp = interp_w_error(olap_wvs_b, olap_wvs_r, olap_flam_r, olap_flam_sig_r)

    olap_flams = np.array([olap_flam_r_interp, olap_flam_b])
    sigs = np.array([olap_flam_sig_r_interp, olap_flam_sig_b])
    weights = sigs ** -2.0

    olap_flam_avgd = np.average(olap_flams, axis=0, weights=weights)
    olap_flam_sig_avgd = 1.0 / np.sqrt(np.mean(weights, axis=0))

    final_wvs = np.concatenate((wvs_b, olap_wvs_b, wvs_r))
    final_flam = np.concatenate((flam_b, olap_flam_avgd, red_mult * flam_r))
    final_flam_sig = np.concatenate((flam_sig_b, olap_flam_sig_avgd, red_mult * flam_sig_r))

    #plt.errorbar(spec_b[1].data['wave'][mask_b], spec_b[1].data['flux'][mask_b], yerr=spec_b[1].data['ivar'][mask_b] ** -0.5, lw=lw)
    #plt.errorbar(spec_r[1].data['wave'][mask_r], red_mult*spec_r[1].data['flux'][mask_r], yerr=red_mult*spec_r[1].data['ivar'][mask_r] ** -0.5, lw=lw)
    #plt.xlabel('Wavelength (Angstroms)')
    #plt.ylabel('Flux (erg/s/cm${}^2$/$\mathring{A}$)')
    #plt.title(f'Fluxed spectrum of {target}')
    #plt.savefig(f'{target}_quickplot.png')


    #plt.figure(figsize=(20,10))
    #plt.errorbar(final_wvs, final_flam, yerr=final_flam_sig)
    #plt.grid()
    #plt.xlabel('Wavelength (Angstroms)')
    #plt.ylabel('Flux (erg/s/cm${}^2$/$\mathring{A}$)')
    #plt.title(f'Fluxed spectrum of {target}')
    #plt.savefig(f'{target}_quickplot.png')

    return ((final_wvs, final_flam, final_flam_sig),
        (spec_r[1].data['wave'], spec_r[1].data['flux'], spec_r[1].data['ivar'] ** -0.5),
        (spec_b[1].data['wave'], spec_b[1].data['flux'], spec_b[1].data['ivar'] ** -0.5))

def interp_w_error(x: np.ndarray, xp: np.ndarray, yp: np.ndarray, err_yp: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if len(xp) == 1:
        return np.ones_like(x) * yp[0], np.ones_like(x) * err_yp[0]
    y = np.zeros_like(x)
    yerr = np.zeros_like(x)
    slopes = np.zeros(xp.shape[0] - 1)
    #slopes = np.zeros(xp.shape[0])
    for i in range(len(slopes)):
        slopes[i] = (yp[i+1] - yp[i])/(xp[i+1] - xp[i])
    #slopes[-1] = slopes[-2]

    for i in range(len(x)):
        # find the index j into xp such that xp[j-1] <= x[i] < xp[j]
        j = np.searchsorted(xp, x[i], side='right')
        if (x[i] == xp[j-1]):
            y[i] = yp[j-1]
            yerr[i] = err_yp[j-1]
        elif (j == len(xp)):
            # extrapolating outside domain!!!
            y[i] = yp[-1]# + slopes[j-2]*(x[i] - xp[-1])
            yerr[i] = np.sqrt((((x[i] - xp[-2])*err_yp[-1]) ** 2 + ((x[i] - xp[-1])*err_yp[-2]) ** 2) / ((xp[-2] - xp[-1]) ** 2))
        elif (j == 0):
            # extrapolating outside domain!!!
            y[i] = yp[0]# + slopes[j]*(x[i] - xp[0])
            yerr[i] = np.sqrt((((x[i] - xp[0])*err_yp[1]) ** 2 + ((x[i] - xp[1])*err_yp[0]) ** 2) / ((xp[1] - xp[0]) ** 2))
        else:
            y[i] = yp[j-1] + slopes[j-1]*(x[i] - xp[j-1])
            yerr[i] = np.sqrt((((x[i] - xp[j])*err_yp[j-1]) ** 2 + ((x[i] - xp[j-1])*err_yp[j]) ** 2) / ((xp[j-1] - xp[j]) ** 2))
    return y, yerr

def save_one2dspec(ax: plt.Axes, spec: np.ndarray, edges: Tuple[np.ndarray, np.ndarray], traces: List[np.ndarray]) -> None:
    all_left, all_right = edges

    norm = ImageNormalize(spec, interval=ZScaleInterval())
    im = ax.imshow(spec, origin='upper', norm=norm, cmap='gray')
    #im, norm = imshow_norm(spec, interval=ZScaleInterval(), cmap='gray')

    plt.axis('off')

    xs = np.arange(spec.shape[0])

    for i in range(all_left.shape[1]):
        ax.plot(all_left[:, i], xs, 'green', lw=1)
        ax.plot(all_right[:, i], xs, 'red', lw=1)

    if traces is not None:
        for trace in traces:
            ax.plot(trace, xs, 'orange', lw=1)

def save_2dspecs(args: dict) -> None:
    obj_png_dict = args['qa_dict']

    arm = 'blue' if 'blue' in args['spectrograph'] else 'red'
    paths = args['output_spec2ds']

    out_path = os.path.join(args['output_path'], 'QA', 'PNGs', 'Extraction')
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    for path in paths:
        # open fits file
        spec = Spec2DObj.from_file(path, 1)
        spec1d_file = path.replace('spec2d', 'spec1d')
        if os.path.isfile(spec1d_file):
            spec1ds = SpecObjs.from_fitsfile(spec1d_file)
        else:
            spec1ds = None
            msgs.warn('Could not find spec1d file: {:s}'.format(spec1d_file) + msgs.newline() +
                  '             No objects were extracted.')

        base_name = os.path.join(out_path, os.path.basename(path).split('_')[1])

        all_left, all_right, _ = spec.slits.select_edges()
        edges = [all_left, all_right]
        traces = [spec1ds[i]['TRACE_SPAT'] for i in range(len(spec1ds))] if spec1ds is not None else None
        mask = spec.bpmmask == 0

        out_shape = spec.sciimg.shape
        fig_width = 0.15 + 5*out_shape[1]/100.
        subfig_width = out_shape[1]/100. / fig_width
        padding_width = 0.15 / fig_width
        fig = plt.figure(figsize=(fig_width, out_shape[0]/100.), dpi=100)

        # science image
        #ax = fig.add_subplot(1, 4, 1)
        ax = fig.add_axes([0, 0, subfig_width, 1])
        save_one2dspec(ax, spec.sciimg, edges, traces)

        # sky model
        ax = fig.add_axes([subfig_width + padding_width, 0, subfig_width, 1])
        save_one2dspec(ax, spec.skymodel * mask, edges, traces)

        # sky subtracted science image
        #ax = fig.add_subplot(1, 4, 2)
        ax = fig.add_axes([2*(subfig_width + padding_width), 0, subfig_width, 1])
        save_one2dspec(ax, (spec.sciimg - spec.skymodel) * mask, edges, traces)

        # sky subtracted images divided by noise
        #ax = fig.add_subplot(1, 4, 3)
        ax = fig.add_axes([3*(subfig_width + padding_width), 0, subfig_width, 1])
        save_one2dspec(ax, (spec.sciimg - spec.skymodel) * np.sqrt(spec.ivarmodel) * mask, edges, traces)

        # total residauls
        #ax = fig.add_subplot(1, 4, 4)
        ax = fig.add_axes([4*(subfig_width + padding_width), 0, subfig_width, 1])
        save_one2dspec(ax, (spec.sciimg - spec.skymodel - spec.objmodel) * np.sqrt(spec.ivarmodel) * mask, edges, traces)

        # save figure
        outname = f'{base_name}_extraction.png'
        msgs.info(f"Saving {os.path.basename(outname)}")
        fig.savefig(outname, bbox_inches='tight', pad_inches=0)
        plt.close(fig)

        airmass = spec.head0['AIRMASS']
        ut = spec.head0['UTSHUT']
        fname = spec.head0['FILENAME']

        if obj_png_dict.get(spec.head0['OBJECT']):
            obj_png_dict[spec.head0['OBJECT']]['pngs'].append(outname)
            obj_png_dict[spec.head0['OBJECT']]['airmass'].append(airmass)
            obj_png_dict[spec.head0['OBJECT']]['time'].append(ut)
            obj_png_dict[spec.head0['OBJECT']]['fname'].append(fname)
        else:
            obj_png_dict[spec.head0['OBJECT']] = {
                'pngs': [outname],
                'airmass': [airmass],
                'time': [ut],
                'fname': [fname]
            }

    args['qa_dict'] = obj_png_dict

def manual_extraction_GUI(args):
    arm = 'blue' if 'blue' in args['spectrograph'] else 'red'
    paths = args['output_spec2ds']

    gui_dict = {}
    for path in paths:
        # open fits file
        spec = Spec2DObj.from_file(path, 1)
        spec1d_file = path.replace('spec2d', 'spec1d')

        if os.path.isfile(spec1d_file):
            spec1ds = SpecObjs.from_fitsfile(spec1d_file)
        else:
            spec1ds = None

        base_name = os.path.basename(path).split('_')[1]

        all_left, all_right, _ = spec.slits.select_edges()
        edges = [all_left, all_right]
        traces = [spec1ds[i]['TRACE_SPAT'] for i in range(len(spec1ds))] if spec1ds is not None else None
        fwhms = [spec1ds[i]['FWHM'] for i in range(len(spec1ds))] if spec1ds is not None else None

        gui_dict[base_name] = {
            'spec': spec,
            'edges': edges,
            'traces': traces,
            'fwhms': fwhms
        }

    # call GUI
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gui = ManualTracingGUI(fig, ax, gui_dict)

    manual_traces = gui.manual_dict

    return manual_traces

def write_manual_pypeit_files(args: dict) -> List[str]:
    """
    Writes pypeit files based on the default pypeit file for this reduction,
    but with filtered targets and additional manual parameter lines.

    requires
    args = {
        'pypeit_file': default pypeit file,
        'targets_list': [[target1, target2], [target3, target4]] will result in
            1 & 2 being reduced together and 3 & 4 being reduced together. The
            target names are blueNNNN-ZTF21abcd.
        'manual_lines_fn': function mapping [target1, target2] to the cfg lines
            they need
        'needs_std_fn': function target1 -> bool, True if target1 needs a
            standard star to be reduced alongside it
    }
    """
    old_pypeit_file = args['pypeit_file']

    targets_list = args['targets_list']
    manual_lines_fn = args['manual_lines_fn']
    needs_std_fn = args['needs_std_fn']

    new_pypeit_files = []
    for targets in targets_list:
        if not targets:
            continue
        target_fnames = [target.split('-')[0] for target in targets]
        new_pypeit_file = f'{os.path.splitext(old_pypeit_file)[0]}_{"_".join(targets)}.pypeit'
        shutil.copy(old_pypeit_file, new_pypeit_file)

        manual_lines = manual_lines_fn(targets)

        cfg_lines = []
        setup_lines = []
        setup = False
        with open(old_pypeit_file, 'r') as old_pypeit_fd:
            for line in old_pypeit_fd.readlines():
                if 'science' in line and '|' in line and all([targ_fname not in line for targ_fname in target_fnames]):
                    pass
                elif 'standard' in line and '|' in line and any(needs_std_fn(target) for target in targets):
                    setup_lines.append(line)
                else:
                    if '# Setup' in line:
                        setup = True
                    if setup:
                        setup_lines.append(line)
                    else:
                        cfg_lines.append(line)

        # Note: the user's parameter file is merged into these config lines
        # so any manual settings user gives might override their GUI actions.
        # Should probably caution user re: this.
        final_cfg = ConfigObj(manual_lines)
        final_cfg.merge(ConfigObj(cfg_lines))
        final_lines = [line + "\n" for line in final_cfg.write()] + setup_lines
        with open(new_pypeit_file, 'w') as new_pypeit_fd:
            new_pypeit_fd.writelines(final_lines)

        new_pypeit_files.append(new_pypeit_file)
    return new_pypeit_files


def manual_extraction(args: dict) -> list:
    manual_dict = manual_extraction_GUI(args)

    args['targets_list'] = [[key] for key in manual_dict.keys()]
    args['manual_lines_fn'] = lambda targ: ['# Added by DBSP_DRP for manual extraction\n',
            '[reduce]\n',
            '[[extraction]]\n',
            '[[[manual]]]\n',
            f"spat_spec = {str(manual_dict[targ[0]]['spat_spec']).strip('[]')}\n",
            f"det = {str([1 for trace in manual_dict[targ[0]]['spat_spec']]).strip('[]')}\n",
            f"fwhm = {str(manual_dict[targ[0]]['fwhm']).strip('[]')}\n",
            "[[skysub]]\n",
            f"user_regions = {str(manual_dict[targ[0]]['bgs']).strip('[]')}\n"
        ]
    args['needs_std_fn'] = lambda targ: manual_dict[targ]['needs_std']
    ret = write_manual_pypeit_files(args)
    # cleanup args so that multiprocessing can pickle it
    del args['targets_list']
    del args['manual_lines_fn']
    del args['needs_std_fn']
    return ret

def re_redux(args: dict, pypeit_files: list) -> None:
    for pypeit_file in pypeit_files:
        these_opts = args.copy()
        these_opts['pypeit_file'] = pypeit_file
        print(f"Using pypeit file {pypeit_file}")
        redux(these_opts)

def write_extraction_QA(args: dict) -> None:
    out_path = os.path.join(args['output_path'], 'QA')

    pngs = [os.path.basename(fullpath) for fullpath in glob.glob(os.path.join(out_path, 'PNGs', 'Extraction', '*.png'))]
    pngs = sorted(pngs)

    objnames = [png.split('-')[1].split('_')[0] for png in pngs]
    qa_dict = args['qa_dict']

    doc, tag, text = Doc().tagtext()
    doc.asis('<!DOCTYPE html>')

    with tag('html'):
        with tag('head'):
            with tag('title'):
                text("DBSP Extraction QA")
            doc.stag('link', rel='stylesheet', href='dbsp_qa.css')
            with tag('script', src='./dbsp_qa.js'):
                pass
        with tag('body'):
            with tag('div'):
                doc.attr(klass='tab')
                for target in qa_dict:
                    with tag('button'):
                        doc.attr(klass='tablinks', onclick=f"openExposure(event, '{target}')")
                        text(target)
            for target, obj_dict in qa_dict.items():
                with tag('div'):
                    doc.attr(klass='tabcontent', id=target)
                    with tag('h1'):
                        text(target)
                    for i, png in enumerate(obj_dict['pngs']):
                        with tag('h2'):
                            text(obj_dict['fname'][i])
                            text('\t')
                            text(obj_dict['time'][i])
                            text('\t')
                            text(obj_dict['airmass'][i])
                        doc.stag('img', src=os.path.join('PNGs', 'Extraction', os.path.basename(png)))

    msgs.info("Writing Extraction QA page")
    result = indent(doc.getvalue())
    extraction_page = os.path.join(out_path, 'Extraction.html')
    with open(extraction_page, mode='wt') as f:
        f.write(result)

    shutil.copy(resource_filename("dbsp_drp", "/data/dbsp_qa.js"), os.path.join(out_path, "dbsp_qa.js"))
    shutil.copy(resource_filename("dbsp_drp", "/data/dbsp_qa.css"), os.path.join(out_path, "dbsp_qa.css"))
    msgs.info(f"Extraction QA page available at {extraction_page}")

def coadd(args: dict) -> List[str]:
    """
    takes in args['spec1dfile'] and coadds each spectrum in the spec1dfile
    takes in args['grouped_spats_list'], a list of dicts mapping 'fnames' to a
        list of filenames and 'spats' to a list of integer spatial pixel
        positions.
    Returns a list of filenames of coadded spectra.
    """
    outfiles = []
    for d in args['grouped_spats_list']:
        fnames = d['fnames']
        spats = d['spats']
        basename = '_'.join([os.path.basename(fname).split("_")[1].split("-")[0]
            for fname in fnames]) + "_" + \
                os.path.basename(fnames[0]).split("_")[1].split("-")[1]

        objnames = []
        for spat, fname in zip(spats, fnames):
            hdul = fits.open(fname)
            for hdu in hdul:
                if f'SPAT{spat:04d}' in hdu.name:
                    objnames.append(hdu.name)

        outfile = os.path.join(args['output_path'], "Science", f"{basename}_{'_'.join(objnames)}.fits")
        coadd_one_object(fnames, objnames, outfile, args)
        outfiles.append(outfile)
    # delete if above works
    #hdul = fits.open(args['spec1dfile'])
    #for hdu in hdul:
    #    if 'SPAT' in hdu.name:
    #        basename = os.path.basename(args['spec1dfile']).split("_")[1]
    #        outfile = os.path.join(args['output_path'], "Science", f"{basename}_{hdu.name}.fits")
    #        coadd_one_object([args['spec1dfile']], [hdu.name], outfile, args)
    #        outfiles.append(outfile)
    return outfiles

def coadd_one_object(spec1dfiles: List[str], objids: List[str], coaddfile: str, args: dict):
    par = load_spectrograph(args['spectrograph']).default_pypeit_par()
    default_cfg_lines = par.to_config()
    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines = default_cfg_lines, merge_with=args['user_config_lines'])
    # Instantiate
    coAdd1d = coadd1d.CoAdd1D.get_instance(spec1dfiles, objids, par=par['coadd1d'], debug=args['debug'], show=args['debug'])
    # Run
    coAdd1d.run()
    # Save to file
    coAdd1d.save(coaddfile)

def telluric_correct(args: dict):
    """
    method to telluric correct one coadded file
    """
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
    outfile = os.path.join(args['output_path'], "Science", (os.path.basename(args['spec1dfile'])).replace('.fits','_tellcorr.fits'))
    modelfile = os.path.join(args['output_path'], "Science", (os.path.basename(args['spec1dfile'])).replace('.fits','_tellmodel.fits'))

    try:
        TelPoly = telluric.poly_telluric(args['spec1dfile'], par['tellfit']['tell_grid'], modelfile, outfile,
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

def group_coadds(fname_to_spats: dict):
    """
    Groups coadds. Destroys input.

    Takes in dict mapping filenames to a list of integer spatial positions
    Returns list of dicts mapping 'fnames' to a list of filenames and 'spats'
        to a list of integer spatial positions.
    """
    # input is dict mapping fname to spats
    # end result is mapping from arb. label of trace -> spats list and fnames list
    THRESHOLD = 2
    result = []
    while any(fname_to_spats.values()):
        potential_group = [(spats[0], fname) for fname, spats in fname_to_spats.items()]
        potential_group.sort(key=lambda x: x[0])
        min_spat, its_fname = potential_group.pop()
        fname_to_spats[its_fname].pop()
        # new group!
        result.append({'spats': [min_spat], 'fnames': [its_fname]})
        # see if any of the others in the potential group are in:
        for spat, fname in potential_group:
            if spat - min_spat < THRESHOLD:
                result[-1]['spats'].append(spat)
                result[-1]['fnames'].append(fname)
                fname_to_spats[fname].pop()

        # filter dict to remove fnames with no spats left
        fname_to_spats = {fname: spats for fname, spats in fname_to_spats.items() if spats}

    return result
