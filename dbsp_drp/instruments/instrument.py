from abc import ABC, abstractmethod
from typing import List, Dict, Tuple, Union

from astropy.table import Table
from astropy.io.fits import HDUList


class Instrument(ABC):
    """
    Abstract class for holding instrument-specific properties and code.
    """

    @property
    @abstractmethod
    def table_edit_columns(self) -> Tuple[str, ...]:
        """
        Tuple of columns names in PypeItMetaData that should be displayed in table_edit.
        """

    @property
    @abstractmethod
    def table_edit_editable_columns(self) -> List[str]:
        """
        List of column names in PypeItMetaData that table_edit should be allowed to edit.
        """

    @abstractmethod
    def column_to_header_kw(self, column: str) -> str:
        """
        Given a column name in PypeItMetaData, returns the corresponding keyword
        in the FITS header.
        """

    @abstractmethod
    def table_edit_hdu_formatter(self, column: str, value: Union[float, str]) -> str:
        """
        Given a PypeItMetaData column name and a corresponding value, formats
        that value for updating a FITS header.

        Args:
            column (str): Column name in PypeItMetaData
            value (Union[float, str]): Value in PypeItMetaData

        Returns:
            str: Formatted value for putting into FITS header.
        """

    @property
    @abstractmethod
    def arms(self) -> int:
        """
        Number of arms the instrument has

        `Arm' is used here to denote a path of light from aperture to focal plane.
        """

    @property
    @abstractmethod
    def arm_prefixes(self) -> List[str]:
        """
        List of file prefixes used by each arm. The order of the arms here
        determines the order that the arms will be reduced in, as well as
        setting the order of all other per-arm instrument attributes. These
        prefixes are used as arm names.

        For example, for DBSP this is return ['blue', 'red'], since DBSP
        filenames are blue0001.fits and red0001.fits.
        """

    @property
    @abstractmethod
    def arm_telluric(self) -> List[bool]:
        """
        List of whether or not each arm needs to be telluric corrected.

        For example, for DBSP this is [False, True], since only the red arm
        needs to be telluric corrected.
        """

    @property
    @abstractmethod
    def arm_names_pypeit(self) -> List[str]:
        """
        List of PypeIt spectrograph names for each arm.

        For example, for DBSP this is ['p200_dbsp_blue', 'p200_dbsp_red'].
        """

    @property
    def pypeit_name_to_arm(self) -> Dict[str, str]:
        """
        Dict mapping the PypeIt spectrograph names to arm names.
        """
        return dict(zip(self.arm_names_pypeit, self.arm_prefixes))

    # time threshholds
    @property
    @abstractmethod
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

    # for coadding, IDing traces as the same
    @property
    @abstractmethod
    def pixel_pos_threshholds(self) -> Dict[str, int]:
        """
        Dict mapping arm names to allowable spatial pixel motion between
        coadded exposures.
        """

    @property
    @abstractmethod
    def frac_pos_threshholds(self) -> Dict[str, float]:
        """
        Dict mapping arm names to allowable spatial motion (as a fraction of
        the slit width) between coadded exposures.
        """

    @property
    @abstractmethod
    def archived_sensfuncs(self) -> List[str]:
        """
        List of paths to archived sensfuncs.
        """

    @abstractmethod
    def calibrate_trace_matching(self, spec1d_table: Table, output_path: str) -> float:
        """
        Calibrates trace-matching across arms. Call this once before calling
        convert_fracpos.

        Args:
            spec1d_table (Table): table of spec1d reduction outputs

        Returns:
            float: matching tolerance
        """

    @abstractmethod
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

    @classmethod
    @abstractmethod
    def detect_instrument(cls, hdulist: HDUList) -> bool:
        """
        Returns True if the input HDUList was taken by this instrument.

        Args:
            hdulist (HDUList)

        Returns:
            bool: True if hdulist was taken by this instrument
        """
