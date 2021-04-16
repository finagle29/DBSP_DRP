"""
Automated splicing for P200 DBSP.
"""

import os
from typing import Tuple, List, Optional
import itertools
import functools

import numpy as np

from astropy.io import fits
import astropy.stats
import astropy.table

import pypeit
import pypeit.pypeit
from pypeit import msgs

import dbsp_drp
from dbsp_drp.instruments import Instrument

def get_raw_header_from_coadd(coadd: str, root: str) -> str:
    if coadd is not None:
        return fits.open(os.path.join(os.path.dirname(root), f"{os.path.basename(coadd).split('-')[0]}.fits"))[0].header
    return fits.Header()

def get_raw_hdus_from_spec1d(spec1d_list: List[str], root: str, output_path: str) -> List[fits.BinTableHDU]:
    ret = []
    for (spec1d, spat) in spec1d_list:
        raw_fname = os.path.join(os.path.dirname(root), f"{os.path.basename(spec1d).split('_')[1].split('-')[0]}.fits")
        # get the raw header
        with fits.open(raw_fname) as raw_hdul:
            raw_header = raw_hdul[0].header.copy()
        # get the spectrum
        spec1d_path = os.path.join(output_path, 'Science', spec1d)
        with fits.open(spec1d_path) as spec1d_hdul:
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


def splice(splicing_dict: dict, root: str, output_path: str, instrument: Instrument) -> None:
    """
    Splices spectra of all arms together.
    """
    for target, targets_dict in splicing_dict.items():
        label = 'a'
        for _, arm_dict in targets_dict.items():
            arm_dicts = [arm_dict.get(arm, {}) for arm in instrument.arm_prefixes]

            arm_files = [d.get('coadd') for d in arm_dicts]
            arm_files = [
                os.path.join(output_path, 'Science', arm_file) if arm_file is not None else None for arm_file in arm_files
            ]

            if not any(arm_files):
                continue

            ((final_wvs, final_flam, final_flam_sig),
                (arm_wvs, arm_flams, arm_sigs)) = adjust_and_combine_overlap_all(arm_files, target)

            primary_header = fits.Header()
            primary_header['DBSP_DRP_V'] = dbsp_drp.__version__
            primary_header['PYPEIT_V'] = pypeit.__version__
            primary_header['NUMPY_V'] = np.__version__
            primary_header['ASTROPY_V'] = astropy.__version__
            for i, arm in enumerate(instrument.arm_prefixes):
                primary_header[f'{arm[0].upper()}_COADD'] = arm_files[i]
            primary_hdu = fits.PrimaryHDU(header=primary_header)

            raw_arm_hdus = [get_raw_hdus_from_spec1d(d.get('spec1ds', []), root, output_path) for d in arm_dicts]

            arm_hdus = [None] * instrument.arms
            for i, arm in enumerate(instrument.arm_prefixes):
                col_wvs = fits.Column(name='wave', array=arm_wvs[i], unit='ANGSTROM', format='D')
                col_flux = fits.Column(name='flux', array=arm_flams[i], unit='E-17 ERG/S/CM^2/ANG', format='D')
                col_error = fits.Column(name='sigma', array=arm_sigs[i], unit='E-17 ERG/S/CM^2/ANG', format='D')
                arm_hdus[i] = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error], name=arm.upper())

            col_wvs = fits.Column(name='wave', array=final_wvs, unit='ANGSTROM', format='D')
            col_flux = fits.Column(name='flux', array=final_flam, unit='E-17 ERG/S/CM^2/ANG', format='D')
            col_error = fits.Column(name='sigma', array=final_flam_sig, unit='E-17 ERG/S/CM^2/ANG', format='D')
            table_hdu = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error], name="SPLICED")

            hdul = fits.HDUList(hdus=[primary_hdu, *itertools.chain(*raw_arm_hdus), *arm_hdus, table_hdu])

            log_msg = f"{target}_{label}.fits contains "
            log_msg += ' and '.join([os.path.basename(arm_file) for arm_file in arm_files if arm_file is not None])
            print(log_msg)

            hdul.writeto(os.path.join(output_path, "Science", f'{target}_{label}.fits'), overwrite=True)
            label = chr(ord(label) + 1)

def adjust_and_combine_overlap_all(arm_files: List[str], target: str) -> Tuple[
            Tuple[np.ndarray, np.ndarray, np.ndarray],
            Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]]:
    # TODO: propagate input masks

    specs = [None] * len(arm_files)
    specs_present = []

    # read in all data
    for i, arm_file in enumerate(arm_files):
        if arm_file is not None:
            specs[i] = fits.open(arm_file)[1].data
            specs_present.append(i)
    # if only one of the specs is not None
    if len(specs_present) == 1:
        ix = specs_present[0]
        return ((specs[ix]['wave'], specs[ix]['flux'], specs[ix]['ivar'] ** -0.5),
            ([spec['wave'] if spec is not None else None for spec in specs],
            [spec['flux'] if spec is not None else None for spec in specs],
            [spec['ivar'] ** -0.5 if spec is not None else None for spec in specs]))

    specs_for_reduction = [
        (spec['wave'], spec['flux'], spec['ivar'] ** -0.5) if spec is not None else None for spec in specs
    ]

    adjust_and_combine_overlap_reduce = functools.partial(adjust_and_combine_overlap, target=target)
    final_wvs, final_flam, final_flam_sig = functools.reduce(adjust_and_combine_overlap_reduce, specs_for_reduction)

    return ((final_wvs, final_flam, final_flam_sig),
            ([spec['wave'] if spec is not None else None for spec in specs],
            [spec['flux'] if spec is not None else None for spec in specs],
            [spec['ivar'] ** -0.5 if spec is not None else None for spec in specs]))

def adjust_and_combine_overlap(bluer_spec: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]],
        redder_spec: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]],
        target: str = "") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if bluer_spec is None:
        return redder_spec
    elif redder_spec is None:
        return bluer_spec

    spec_r_wave, spec_r_flux, spec_r_sig = redder_spec
    spec_b_wave, spec_b_flux, spec_b_sig = bluer_spec

    # combination steps
    overlap_lo = spec_r_wave[0]
    overlap_hi = spec_b_wave[-1]
    # maybe need to fix the overlaps?
    # need more finely spaced grid to be completely contained within coarser grid

    if overlap_lo > overlap_hi:
        # there is no overlap!
        # we can't adjust the flux level
        # so we just concatenate!
        final_wvs = np.concatenate([spec_b_wave, spec_r_wave])
        final_flam = np.concatenate([spec_b_flux, spec_r_flux])
        final_flam_sig = np.concatenate([spec_b_sig, spec_r_sig])
        return (final_wvs, final_flam, final_flam_sig)


    olap_r = (spec_r_wave < overlap_hi)
    olap_b = (spec_b_wave > overlap_lo)

    red_mult = 1
    red_mult = (astropy.stats.sigma_clipped_stats(spec_b_flux[olap_b])[1] /
        astropy.stats.sigma_clipped_stats(spec_r_flux[olap_r])[1])
    #red_mult = np.average(spec_aag[1].data['OPT_FLAM'][olap_b], weights=spec_aag[1].data['OPT_FLAM_SIG'][olap_b] ** -2.0)\
    #   /np.average(spec_aag_red[1].data['OPT_FLAM'][olap_r], weights=spec_aag_red[1].data['OPT_FLAM_SIG'][olap_r] ** -2.0)
    if red_mult > 3 or 1/red_mult > 3:
        msgs.warn(f"For {target}, red spectrum is {red_mult} times less flux than blue spectrum in overlap region." +
            "The red and blue traces may not correspond to the same object.")


    # different dispersion.
    wvs_b = spec_b_wave[~olap_b]
    wvs_r = spec_r_wave[~olap_r]
    flam_b = spec_b_flux[~olap_b]
    flam_r = spec_r_flux[~olap_r]
    flam_sig_b = spec_b_sig[~olap_b]
    flam_sig_r = spec_r_sig[~olap_r]


    olap_wvs_r = spec_r_wave[olap_r]
    olap_flam_r = red_mult * spec_r_flux[olap_r]
    olap_flam_sig_r = red_mult * spec_r_sig[olap_r]
    olap_wvs_b = spec_b_wave[olap_b][:-1]
    olap_flam_b = spec_b_flux[olap_b][:-1]
    olap_flam_sig_b = spec_b_sig[olap_b][:-1]

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

    return (final_wvs, final_flam, final_flam_sig)

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
