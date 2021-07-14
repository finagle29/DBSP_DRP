import os
import io
import glob
from datetime import datetime

from astropy.io import fits
from astropy.coordinates import Angle, SkyCoord, AltAz, EarthLocation
from astropy.time import Time
import astropy.units as u


def _get_ra_or_dec(fname, is_ra, header, modified_str, prompt_user) -> bool:
    kw = 'RA' if is_ra else 'DEC'
    unit = 'hour' if is_ra else 'deg'
    coord_str = None
    try:
        coord_str = header[kw]
        try:
            coord = Angle(coord_str, unit=unit)
            return True
        except ValueError as ve:
            print(f"File {fname} has {kw} that could not be parsed: {coord_str}.")
            raise KeyError from ve
    except KeyError:
        if prompt_user:
            if coord_str is not None:
                prompt_str = f"File {fname} has bad {kw} keyword in header. Bad value is {coord_str}.\n"
            else:
                prompt_str = f"File {fname} is missing {kw} keyword in header. "
            prompt_str += ("Header may be severely malformed.\n"
                    f"Enter the {kw} of `{header['OBJECT']}` or an empty string if not known: ")
            coord_str = input(prompt_str)
            if coord_str != '':
                try:
                    coord = Angle(coord_str, unit=unit)

                    header[kw] = coord.to_string(unit=unit, sep=':')
                    header.comments[kw] += modified_str + ", manually entered."
                    return True
                except ValueError:
                    header[kw] = coord_str
                    header.comments[kw] += modified_str + ", manually entered, but could not be parsed into an angle."
                    return False
        else:
            if coord_str is None:
                print(f"File {fname} is missing {kw} keyword in header. Be sure to correct this!")
            return False

def main(path_prefix, throw_errors, prompt_user):
    dt = datetime.now()
    modified_str = f" Modified {dt.isoformat(timespec='seconds')}"

    good_files = []

    paths = glob.glob(os.path.join(path_prefix, '*.fits'))
    if not paths:
        paths = glob.glob(f'{path_prefix}*.fits')
    for path in paths:
        fname = os.path.basename(path)
        if (os.path.getsize(path) % 2880):
            size = os.path.getsize(path)
            print(f"WARNING: {fname} is {size} bytes long, which is NOT a multiple of 2880.")
            print("Therefore it is an invalid FITS file.")
            print("It will now be truncated to the next smallest multiple of 2880 bytes long")
            with open(path, 'r+b') as f:
                f.seek(-(size % 2880), io.SEEK_END)
                f.truncate()
        if os.path.getsize(path) == 0:
            print(f"WARNING: {fname} is empty, it is 0 bytes.")
            continue
        with fits.open(path, mode='update') as hdul:
            header = hdul[0].header
            try:
                ang = Angle(header['ANGLE'].lower(), unit='deg').to_string(unit='deg', sep=(' deg ', ' min'), fields=2)
            except ValueError as exc:
                try:
                    ang = Angle(header['GRATING'].lower()).to_string(unit='deg', sep=(' deg ', ' min'), fields=2)
                    grat = header['ANGLE']

                    print(f"In {fname} ANGLE and GRATING were swapped, correcting them.")

                    header['ANGLE'] = ang
                    header.comments['ANGLE'] += modified_str + ", was swapped with GRATING."

                    header['GRATING'] = grat
                    header.comments['GRATING'] += modified_str + ", was swapped with ANGLE."
                except:
                    msg = (f"In {fname} ANGLE in header ({header['ANGLE']}) is not parseable as an angle" +
                        f"and neither is GRATING ({header['GRATING']}), the usual suspect.")
                    if throw_errors:
                        raise ValueError(msg) from exc
                    else:
                        print(msg)
            if header.get('IMGTYPE', '') == 'object':
                # we NEED to have parseable RA/DECs
                ra_good = _get_ra_or_dec(fname, True, header, modified_str, True)
                while not ra_good:
                    print(f"File {fname} has IMGTYPE object and therefore must have valid RA.")
                    print("Though RA can be corrected later if you ran dbsp_reduce without -i")
                    ra_good = _get_ra_or_dec(fname, True, header, modified_str, True)

                dec_good = _get_ra_or_dec(fname, False, header, modified_str, True)
                while not dec_good:
                    print(f"File {fname} has IMGTYPE object and therefore must have valid DEC.")
                    print("Though DEC can be corrected later if you ran dbsp_reduce without -i")
                    dec_good = _get_ra_or_dec(fname, False, header, modified_str, True)
            else:
                ra_good = _get_ra_or_dec(fname, True, header, modified_str, prompt_user)
                dec_good = _get_ra_or_dec(fname, False, header, modified_str, prompt_user)

            if (ra_good and dec_good and {'RA', 'DEC', 'UTSHUT'}.issubset(header)
                and (header.get('AIRMASS', None) is None)):
                print(f'In {fname} RA and DEC are present, but AIRMASS is absent. Calculating airmass from RA/DEC.')
                skycoord = SkyCoord(header['RA'], header['DEC'], unit=(u.hour, u.deg))
                time = Time(header['UTSHUT'])
                altaz = AltAz(obstime=time, location=EarthLocation.of_site('Palomar'))

                header['AIRMASS'] = f'{skycoord.transform_to(altaz).secz:.3f}'
                header.comments['AIRMASS'] += modified_str + ", calculated from present RA/DEC."
        good_files.append(path)

    return good_files
