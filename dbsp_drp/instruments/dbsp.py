import os
from typing import List, Dict

import numpy as np
from astropy.table import Table
from astropy.io import fits

from dbsp_drp.instruments import Instrument

class DBSP(Instrument):
    """
    Class implementing instrument-specific properties for DBSP.
    """

    @property
    def arms(self) -> int:
        """
        Number of arms the instrument has

        `Arm' is used here to denote a path of light from aperture to focal plane.
        """
        return 2

    @property
    def arm_prefixes(self) -> List[str]:
        """
        List of file prefixes used by each arm. The order of the arms here
        determines the order that the arms will be reduced in, as well as
        setting the order of all other per-arm instrument attributes. These
        prefixes are used as arm names.

        For example, for DBSP this is return ['blue', 'red'], since DBSP
        filenames are blue0001.fits and red0001.fits.
        """
        return ['blue', 'red']

    @property
    def arm_telluric(self) -> List[bool]:
        """
        List of whether or not each arm needs to be telluric corrected.

        For example, for DBSP this is [False, True], since only the red arm
        needs to be telluric corrected.
        """
        return [False, True]

    @property
    def arm_names_pypeit(self) -> List[str]:
        """
        List of PypeIt spectrograph names for each arm.

        For example, for DBSP this is ['p200_dbsp_blue', 'p200_dbsp_red'].
        """
        return ['p200_dbsp_blue', 'p200_dbsp_red']

    # time threshholds
    @property
    def coadd_threshholds(self) -> Dict[str, float]:
        """
        Dict mapping arm names to allowable time between consecutive exposures
        for those exposures to be coadded.

        TODO:
        Special value 0 indicates allow arbitrary time between consecutive
        exposures. Special value -1 indicates all exposures of the same object
        will be coadded, even if nonconsecutive.

        TODO:
        I think this needs a setter.
        """
        return {
            'blue': 60.,
            'red': 60.
        }

    # for coadding, IDing traces as the same
    @property
    def pixel_pos_threshholds(self) -> Dict[str, int]:
        """
        Dict mapping arm names to allowable spatial pixel motion between
        coadded exposures.
        """
        return {
            'blue': 2,
            'red': 2
        }

    @property
    def frac_pos_threshholds(self) -> Dict[str, float]:
        """
        Dict mapping arm names to allowable spatial motion (as a fraction of
        the slit width) between coadded exposures.
        """
        ## TODO: this isn't used, maybe remove.
        return {
            'blue': 0.01,
            'red': 0.01
        }

    def _get_std_trace(self, std_path: str) -> float:
        max_sn = -1
        max_fracpos = -1
        with fits.open(std_path) as hdul:
            # loop through trace hdus
            for hdu in hdul:
                if not 'SPAT' in hdu.name:
                    continue

                # look at s/n
                if 'OPT_COUNTS' in hdu.data.dtype.names:
                    this_sn = np.nanmedian(hdu.data['OPT_COUNTS']/hdu.data['OPT_COUNTS_SIG'])
                elif 'BOX_COUNTS' in hdu.data.dtype.names:
                    this_sn = np.nanmedian(hdu.data['BOX_COUNTS']/hdu.data['BOX_COUNTS_SIG'])
                else:
                    this_sn = -1

                if this_sn > max_sn:
                    max_sn = this_sn
                    max_fracpos = hdu.header['SPAT_FRACPOS']

        if max_fracpos == -1:
            raise Exception(f"Error! No HDUs in {os.path.basename(std_path)} have median S/N > 0.")
        return max_fracpos

    def calibrate_trace_matching(self, spec1d_table: Table, output_path: str) -> float:
        """
        Calibrates trace-matching across arms. Call this once before calling
        convert_fracpos.

        Args:
            spec1d_table (Table): table of spec1d reduction outputs

        Returns:
            float: matching tolerance
        """
        stds = spec1d_table['frametype'] == 'standard'
        red_mask = spec1d_table['arm'] == 'red'
        blue_mask = spec1d_table['arm'] == 'blue'

        FRACPOS_SUM = 1.0
        FRACPOS_TOL = 0.025

        if (stds & blue_mask).any() and (stds & red_mask).any():
            std_fracpos_sums = []
            for row in spec1d_table[stds]:
                # find closest mjd frame of other arm
                if not spec1d_table.loc[row['filename']]['processed']:
                    other_arm = spec1d_table['arm'] != row['arm']
                    corresponding_row = spec1d_table[other_arm & stds][np.abs(spec1d_table[other_arm & stds]['mjd'] - row['mjd']).argmin()]
                    this_path = os.path.join(output_path, 'Science', row['filename'])
                    corresponding_path = os.path.join(output_path, 'Science', corresponding_row['filename'])
                    std_fracpos_sums.append(self._get_std_trace(this_path) +
                        self._get_std_trace(corresponding_path))
                    spec1d_table.loc[row['filename']]['processed'] = True
                    spec1d_table.loc[corresponding_row['filename']]['processed'] = True
            FRACPOS_SUM = np.mean(std_fracpos_sums)
            FRACPOS_TOL = FRACPOS_SUM * .025

        self.FRACPOS_SUM = FRACPOS_SUM
        self.FRACPOS_TOL = FRACPOS_TOL
        return FRACPOS_TOL

    def convert_fracpos(self, arm: str, fracpos: float) -> float:
        """
        Maps (arm, fracpos) to a float such that the traces of the same object
        on different arms will have the same (or close) output from this
        function.

        Args:
            arm (str): spectrograph arm
            fracpos (float): fractional position across the slit

        Returns:
            float: object identifier
        """
        if arm == 'blue':
            return self.FRACPOS_SUM - fracpos
        return fracpos

    @classmethod
    def detect_instrument(cls, hdulist: fits.HDUList) -> bool:
        """
        Returns True if the input HDUList was taken by this instrument.

        Args:
            hdulist (HDUList)

        Returns:
            bool: True if hdulist was taken by this instrument
        """
        return 'dbsp' in hdulist[0].header['FPA'].lower()
