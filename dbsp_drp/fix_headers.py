import os
import io
import sys
import glob

from astropy.io import fits
from astropy.coordinates import Angle


def _get_ra_or_dec(path, is_ra, header):
    kw = 'RA' if is_ra else 'DEC'
    unit = 'hour' if is_ra else 'deg'
    try:
        coord = header[kw]
    except KeyError:
        coord_str = input(f"File {path} is missing {kw} keyword in header. Header may be severely malformed.\
             Enter the {kw} of {header['OBJECT']} or an empty string if not known: ")
        if coord_str != '':
            try:
                coord = Angle(coord_str, unit=unit)
                header[kw] = coord.to_string(unit=unit, sep=':')
            except ValueError:
                header[kw] = coord_str

def main(basepath):
    paths = glob.glob(os.path.join(basepath, '*.fits'))
    if not paths:
        paths = glob.glob(f'{basepath}*.fits')
    for path in paths:
        if (os.path.getsize(path) % 2880):
            size = os.path.getsize(path)
            print(f"WARNING: {os.path.basename(path)} is {size} bytes long, which is NOT a multiple of 2880.")
            print("Therefore it is an invalid FITS file.")
            print("It will now be truncated to the next smallest multiple of 2880 bytes long")
            with open(path, 'r+b') as f:
                f.seek(-(size % 2880), io.SEEK_END)
                f.truncate()
        with fits.open(path, mode='update') as hdul:
            header = hdul[0].header
            try:
                ang = Angle(header['ANGLE'], unit='deg').to_string(unit='deg', sep=(' deg ', ' min'), fields=2)
            except ValueError:
                try:
                    ang = Angle(header['GRATING']).to_string(unit='deg', sep=(' deg ', ' min'), fields=2)
                    grat = header['ANGLE']
                    header['ANGLE'] = ang
                    header['GRATING'] = grat
                except:
                    raise ValueError(f"ANGLE in header ({header['ANGLE']}) is not parseable as an angle\
                                    and neither is GRATING ({header['GRATING']}), the usual suspect")
            _get_ra_or_dec(path, True, header)
            _get_ra_or_dec(path, False, header)
