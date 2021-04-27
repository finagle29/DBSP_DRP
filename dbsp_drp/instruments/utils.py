import sys

from astropy.io import fits
from dbsp_drp import instruments

def guess_instrument_from_file(path: str) -> instruments.Instrument:
    with fits.open(path) as hdul:
        for ins in instruments.instruments:
            if ins.detect_instrument(hdul):
                print(f"Automatically detected instrument: {ins.__name__}.")
                return ins()

        print("ERROR: Unknown spectrograph!")
        print(f"The instrument for reduction could not be determined from {path}.")
        print("Please use the --instrument flag to specify the instrument for reduction.")
        sys.exit(1)
