import os
import sys
import glob

from astropy.io import fits
from astropy.coordinates import Angle


def _get_ra_or_dec(is_ra, header):
    kw = 'RA' if is_ra else 'DEC'
    unit = 'hour' if is_ra else 'deg'
    try:
        coord = header[kw]
    except:
        coord_str = input(f"File {path} is missing {kw} keyword in header. Header may be severely malformed.\
             Enter the {kw} of {header['OBJECT']} or an empty string if not known.")
        if coord_str != '':
            try:
                coord = Angle(coord_str, unit=unit)
                header[{kw}] = coord.to_string()
            except:
                header[{kw}] = coord_str

def main(basepath):
    paths = glob.glob(os.path.join(basepath, '*.fits'))
    if not paths:
        paths = glob.glob(f'{basepath}*.fits')
    for path in paths:
        with fits.open(path, mode='update', ignore_missing_end=True) as hdul:
            header = hdul[0].header
            try:
                ang = Angle(header['ANGLE'], unit='deg').to_string()
                header['ANGLE'] = ang
            except ValueError:
                try:
                    ang = Angle(header['GRATING']).to_string()
                    grat = header['ANGLE']
                    header['ANGLE'] = ang
                    header['GRATING'] = grat
                except:
                    raise ValueError(f"ANGLE in header ({header['ANGLE']}) is not parseable as an angle\
                                    and neither is GRATING ({header['GRATING']}), the usual suspect")
            _get_ra_or_dec(True, header)
            _get_ra_or_dec(False, header)
