"""
Automated splicing for P200 DBSP.
"""

import os
from typing import Tuple, List

import numpy as np

from astropy.io import fits
import astropy.stats
import astropy.table

import pypeit
import pypeit.pypeit
from pypeit import msgs

import dbsp_drp

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


def splice(splicing_dict: dict, root: str, output_path: str) -> None:
    """
    Splices red and blue spectra together.
    """
    for target, targets_dict in splicing_dict.items():
        label = 'a'
        for _, arm_dict in targets_dict.items():
            red_dict = arm_dict.get('red', {})
            blue_dict = arm_dict.get('blue', {})

            bluefile = blue_dict.get('coadd')
            redfile = red_dict.get('coadd')
            if bluefile is not None:
                bluefile = os.path.join(output_path, 'Science', bluefile)
            if redfile is not None:
                redfile = os.path.join(output_path, 'Science', redfile)
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
            primary_header['B_COADD'] = bluefile
            primary_header['R_COADD'] = redfile
            primary_hdu = fits.PrimaryHDU(header=primary_header)

            raw_red_hdus = get_raw_hdus_from_spec1d(red_dict.get('spec1ds', []), root, output_path)
            raw_blue_hdus = get_raw_hdus_from_spec1d(blue_dict.get('spec1ds', []), root, output_path)


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

            log_msg = f"{target}_{label}.fits contains "
            if redfile is None:
                log_msg += f"{os.path.basename(bluefile)}"
            elif bluefile is None:
                log_msg += f"{os.path.basename(redfile)}"
            else:
                log_msg += f"{os.path.basename(redfile)} and {os.path.basename(bluefile)}"
            print(log_msg)
            hdul.writeto(os.path.join(output_path, "spliced", f'{target}_{label}.fits'), overwrite=True)
            label = chr(ord(label) + 1)

def adjust_and_combine_overlap(bluefile: str, redfile: str, target: str) -> Tuple[
            Tuple[np.ndarray, np.ndarray, np.ndarray],
            Tuple[np.ndarray, np.ndarray, np.ndarray],
            Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    # TODO: propagate input masks

    if bluefile is not None:
        spec_b = fits.open(bluefile)[1].data
        if redfile is None:
            return ((spec_b['wave'], spec_b['flux'], spec_b['ivar'] ** -0.5),
                (None, None, None),
                (spec_b['wave'], spec_b['flux'], spec_b['ivar'] ** -0.5))
    if redfile is not None:
        spec_r = fits.open(redfile)[1].data
        if bluefile is None:
            return ((spec_r['wave'], spec_r['flux'], spec_r['ivar'] ** -0.5),
                (spec_r['wave'], spec_r['flux'], spec_r['ivar'] ** -0.5),
                (None, None, None))

    # combination steps
    overlap_lo = spec_r['wave'][0]
    overlap_hi = spec_b['wave'][-1]
    # maybe need to fix the overlaps?
    # need more finely spaced grid to be completely contained within coarser grid

    olap_r = (spec_r['wave'] < overlap_hi)
    olap_b = (spec_b['wave'] > overlap_lo)

    red_mult = 1
    red_mult = (astropy.stats.sigma_clipped_stats(spec_b['flux'][olap_b])[1] /
        astropy.stats.sigma_clipped_stats(spec_r['flux'][olap_r])[1])
    #red_mult = np.average(spec_aag[1].data['OPT_FLAM'][olap_b], weights=spec_aag[1].data['OPT_FLAM_SIG'][olap_b] ** -2.0)\
    #   /np.average(spec_aag_red[1].data['OPT_FLAM'][olap_r], weights=spec_aag_red[1].data['OPT_FLAM_SIG'][olap_r] ** -2.0)
    if red_mult > 3 or 1/red_mult > 3:
        msgs.warn(f"For {target}, red spectrum is {red_mult} times less flux than blue spectrum in overlap region." +
            "The red and blue traces may not correspond to the same object.")


    # different dispersion.
    wvs_b = spec_b['wave'][~olap_b]
    wvs_r = spec_r['wave'][~olap_r]
    flam_b = spec_b['flux'][~olap_b]
    flam_r = spec_r['flux'][~olap_r]
    flam_sig_b = spec_b['ivar'][~olap_b] ** -0.5
    flam_sig_r = spec_r['ivar'][~olap_r] ** -0.5


    olap_wvs_r = spec_r['wave'][olap_r]
    olap_flam_r = red_mult * spec_r['flux'][olap_r]
    olap_flam_sig_r = red_mult * spec_r['ivar'][olap_r] ** -0.5
    olap_wvs_b = spec_b['wave'][olap_b][:-1]
    olap_flam_b = spec_b['flux'][olap_b][:-1]
    olap_flam_sig_b = spec_b['ivar'][olap_b][:-1] ** -0.5

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

    return ((final_wvs, final_flam, final_flam_sig),
        (spec_r['wave'], spec_r['flux'], spec_r['ivar'] ** -0.5),
        (spec_b['wave'], spec_b['flux'], spec_b['ivar'] ** -0.5))

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
