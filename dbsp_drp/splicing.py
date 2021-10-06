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

def get_raw_hdus_from_spec1d(spec1d_list: List[Tuple[str, int]], root: str,
        output_path: str) -> List[fits.BinTableHDU]:
    """
    Returns list of ``fits.BinTableHDU`` s, each containing the raw header and
    the 1D spectrum, of the input spec1d files.

    Args:
        spec1d_list (List[Tuple[str, int]]): List of (spec1d filename, spatial
            pixel coordinate)
        root (str): Path to raw data files, possibly including filename prefix.
        output_path (str): reduction output path

    Returns:
        List[fits.BinTableHDU]: List of raw data headers and data from input
            spec1d files.
    """
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


def splice(splicing_dict: dict, interpolate_gaps: bool, root: str, output_path: str) -> None:
    """
    Splices red and blue spectra together.

    .. code-block::

        splicing_dict[target_name][position_along_slit][arm] = {
            'spec1ds': [(spec1d_filename_1, spatial_pixel_1), (spec1d_filename_2, spatial_pixel_2)],
            'coadd': coadd_filename
        }

    Args:
        splicing_dict (dict): Guides splicing.
        interpolate_gaps (bool): Interpolate across gaps in wavelength coverage?
        root (str): Path to raw data files, possibly including filename prefix.
        output_path (str): reduction output path
    """
    for target, targets_dict in splicing_dict.items():
        label = 'a'
        for _, arm_dict in targets_dict.items():
            red_dict = arm_dict.get('red', {})
            blue_dict = arm_dict.get('blue', {})

            bluefile = blue_dict.get('coadd')
            redfile = red_dict.get('coadd')
            spec_b = None
            spec_r = None
            if bluefile is not None:
                bluefile = os.path.join(output_path, 'Science', bluefile)
                spec_b = fits.open(bluefile)[1].data
            if redfile is not None:
                redfile = os.path.join(output_path, 'Science', redfile)
                spec_r = fits.open(redfile)[1].data
            if bluefile is None and redfile is None:
                continue

            ((final_wvs, final_flam, final_flam_sig),
                (red_wvs, red_flam, red_sig),
                (blue_wvs, blue_flam, blue_sig)) = adjust_and_combine_overlap(spec_b, spec_r, interpolate_gaps)

            primary_header = fits.Header()
            primary_header['HIERARCH DBSP_DRP_V'] = dbsp_drp.__version__
            primary_header['PYPEIT_V'] = pypeit.__version__
            primary_header['NUMPY_V'] = np.__version__
            primary_header['HIERARCH ASTROPY_V'] = astropy.__version__
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

            table_hdu.header['HIERARCH INTERP_GAPS'] = interpolate_gaps

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

def adjust_and_combine_overlap(
    spec_b: fits.FITS_rec,
    spec_r: fits.FITS_rec,
    interpolate_gaps: bool,
    red_mult: float = 1.0
) -> Tuple[
        Tuple[np.ndarray, np.ndarray, np.ndarray],
        Tuple[np.ndarray, np.ndarray, np.ndarray],
        Tuple[np.ndarray, np.ndarray, np.ndarray]
]:
    """
    Takes in red and blue spectra, adjusts overall flux level by red_mult, and
    combines overlap region.

    In the overlap region, the red spectrum is linearly interpolated to match
    the blue spectrum's wavelength spacing.

    Args:
        spec_b (fits.FITS_rec): blue spectrum
        spec_r (fits.FITS_rec): red spectrum.
        interpolate_gaps (bool): Interpolate across gaps in wavelength coverage?
        red_mult (float, optional): Factor multiplied into the red spectrum to
            match overal flux level with the blue spectrum. Defaults to 1.0.

    Raises:
        ValueError: Raised when both `spec_b` and `spec_r` are empty or None.

    Returns:
        Tuple[
            Tuple[np.ndarray, np.ndarray, np.ndarray],
            Tuple[np.ndarray, np.ndarray, np.ndarray],
            Tuple[np.ndarray, np.ndarray, np.ndarray]
        ]: (blue, red, combined) spectra, where each spectrum is a tuple of
            (wavelengths, flux, error)
    """
    if ((spec_b is None or not spec_b['wave'].shape[0]) and
        (spec_r is None or not spec_r['wave'].shape[0])):
        raise ValueError("Both arguments cannot be empty or None.")
    # TODO: propagate input masks
    if spec_r is None or not spec_r['wave'].shape[0]:
        return ((spec_b['wave'], spec_b['flux'], spec_b['ivar'] ** -0.5),
                (None, None, None),
                (spec_b['wave'], spec_b['flux'], spec_b['ivar'] ** -0.5))
    if spec_b is None or not spec_b['wave'].shape[0]:
        return ((spec_r['wave'], red_mult*spec_r['flux'], red_mult*spec_r['ivar'] ** -0.5),
                (spec_r['wave'], spec_r['flux'], spec_r['ivar'] ** -0.5),
                (None, None, None))

    # combination steps
    overlap_lo = spec_r['wave'][0]
    overlap_hi = spec_b['wave'][-1]
    # maybe need to fix the overlaps?
    # need more finely spaced grid to be completely contained within coarser grid

    if overlap_lo > overlap_hi:
        # there is no overlap!
        # we can't adjust the flux level
        # so we just concatenate!
        final_wvs = np.concatenate([spec_b['wave'], spec_r['wave']])
        final_flam = np.concatenate([spec_b['flux'], spec_r['flux']*red_mult])
        final_flam_sig = np.concatenate([spec_b['ivar'] ** -0.5, (spec_r['ivar'] ** -0.5) * red_mult])
        return ((final_wvs, final_flam, final_flam_sig),
            (spec_r['wave'], spec_r['flux'], spec_r['ivar'] ** -0.5),
            (spec_b['wave'], spec_b['flux'], spec_b['ivar'] ** -0.5))

    olap_r = (spec_r['wave'] < overlap_hi)
    olap_b = (spec_b['wave'] > overlap_lo)

    ## 05/25/2021 red_mult is not really necessary, spectra look better without it.
    ## 06/25/2021 keeping red_mult as an argument for manual_splicing
    #red_mult = (astropy.stats.sigma_clipped_stats(spec_b['flux'][olap_b])[1] /
    #    astropy.stats.sigma_clipped_stats(spec_r['flux'][olap_r])[1])


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

    olap_flam_r_interp, olap_flam_sig_r_interp = interp_w_error(olap_wvs_b, olap_wvs_r, olap_flam_r, olap_flam_sig_r, interpolate_gaps)

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

def interp_w_error(x: np.ndarray, xp: np.ndarray, yp: np.ndarray,
    err_yp: np.ndarray, interpolate_gaps: bool) -> Tuple[np.ndarray, np.ndarray]:
    """
    Linearly interpolate the data points (``xp``, ``yp``) with ``err_yp``
    uncertainty onto the grid ``x``.

    Args:
        x (np.ndarray): destination x data
        xp (np.ndarray): source x data
        yp (np.ndarray): source y data
        err_yp (np.ndarray): source y error data
        interpolate_gaps (bool): Interpolate across gaps in ``xp``?

    Returns:
        Tuple[np.ndarray, np.ndarray]: Interpolated y and error.
    """
    if len(xp) == 1:
        return np.ones_like(x) * yp[0], np.ones_like(x) * err_yp[0]

    y = np.zeros_like(x)
    yerr = np.zeros_like(x)
    slopes = np.zeros(xp.shape[0] - 1)

    dxp = np.diff(xp)
    mean_dxp, _, _ = astropy.stats.sigma_clipped_stats(dxp)

    for i in range(len(slopes)):
        slopes[i] = (yp[i+1] - yp[i])/dxp[i]
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
            # If we are interpolating a gap larger than 5 times the avg d lambda
            if (xp[j] - xp[j-1]) > mean_dxp * 5:
                if interpolate_gaps:
                    # err is max of edge points
                    yerr[i] = max(err_yp[j], err_yp[j-1])
                else:
                    # err is infinite, so this point is completely discounted
                    yerr[i] = np.inf
            else:
                yerr[i] = np.sqrt((((x[i] - xp[j])*err_yp[j-1]) ** 2 + ((x[i] - xp[j-1])*err_yp[j]) ** 2) / ((xp[j-1] - xp[j]) ** 2))
    return y, yerr
