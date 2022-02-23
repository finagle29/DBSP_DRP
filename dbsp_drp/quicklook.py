import argparse
import os
import time
import sys
from typing import List, Optional
import glob
from multiprocessing import Process

import numpy as np
from astropy.io import fits

from pkg_resources import resource_filename

from pypeit.pypeitsetup import PypeItSetup
from pypeit.core import framematch
from pypeit import pypeit
from pypeit import fluxcalibrate
from pypeit.scripts import show_2dspec

from dbsp_drp import show_spectrum

def entrypoint():
    main(parse())

def get_cfg_lines(spectrograph: str) -> List[str]:
    """
    Get standard quicklook PypeIt configuration for ``spectrograph``.

    Args:
        spectrograph (str): PypeIt name of spectrograph.

    Returns:
        List[str]: Standard PypeIt quicklook configuration for ``spectrograph``.
    """
    cfg_lines = [
        "[rdx]",
        f"spectrograph = {spectrograph}",
        "[calibrations]",
        f"master_dir = Master_{spectrograph.split('_')[-1]}",
        "raise_chk_error = False",
        "[scienceframe]",
        "[[process]]",
        "mask_cr = False",
        "[baseprocess]",
        "use_biasimage = False",
        "[reduce]",
        "[[extraction]]",
        "skip_optimal = True",
        "[[findobj]]",
        "skip_second_find = True"
    ]
    return cfg_lines

def parse(options: Optional[List[str]] = None) -> argparse.Namespace:
    argparser = argparse.ArgumentParser(description="Quicklook for P200 DBSP",
        formatter_class=argparse.RawTextHelpFormatter)

    argparser.add_argument("fname", type=str, nargs="+",
        help="file to take a quick look at, or else red/blue\n"
        "to just perform rough calibrations")

    argparser.add_argument("--no-show", default=False, action="store_true",
        help="Set this flag to suppress opening of plots")

    return argparser.parse_args() if options is None else argparser.parse_args(options)

def main(args: argparse.Namespace):
    t = time.perf_counter()
    # need an arc frame and a flat frame
    if all('red' in os.path.basename(fname) for fname in args.fname):
        spectrograph = 'p200_dbsp_red'
    elif all('blue' in os.path.basename(fname) for fname in args.fname):
        spectrograph = 'p200_dbsp_blue'
    else:
        sys.exit(f"Raw files must be from same spectrograph: {args.fname}.")

    arm = spectrograph.split('_')[-1]

    CFG_LINES = get_cfg_lines(spectrograph)

    flatimg = ""
    arcimg = ""
    sciimg = args.fname

    calib_only = any(not os.path.isfile(fname) for fname in args.fname)

    if calib_only:
        root = args.fname[0].rstrip('0123456789.fits')
        paths = glob.glob(f'{root}*.fits')

        for path in paths:
            with fits.open(path) as hdul:
                if not flatimg:
                    if hdul[0].header['OBJECT'] == 'flat' or hdul[0].header['IMGTYPE'] == 'flat':
                        flatimg = path
                if not arcimg:
                    if hdul[0].header['OBJECT'] == 'arcs' or hdul[0].header['IMGTYPE'] == 'cal':
                        arcimg = path
                if flatimg and arcimg:
                    break

        if not (flatimg and arcimg):
            raise Exception(f"Could not find a flat and an arc frame in the same directory as {root}!")
        files = [arcimg, flatimg]
    else:
        files = args.fname

    ps = PypeItSetup(files, path="./", spectrograph_name=spectrograph,
        cfg_lines = CFG_LINES)
    ps.build_fitstbl(strict=False)

    bm = framematch.FrameTypeBitMask()
    file_bits = np.zeros(len(files), dtype=bm.minimum_dtype())
    if calib_only:
        file_bits[0] = bm.turn_on(file_bits[0], ['arc', 'tilt'])
        file_bits[1] = bm.turn_on(file_bits[1], ['pixelflat', 'trace', 'illumflat'])
    else:
        for i in range(len(files)):
            file_bits[i] = bm.turn_on(file_bits[i], 'science')
        ps.fitstbl['calib_id'] = 1

    asrt = np.array([ps.fitstbl['filename'].data.tolist().index(os.path.basename(fname)) for fname in files])
    ps.fitstbl.set_frame_types(file_bits[asrt])
    ps.fitstbl.set_combination_groups()

    ps.fitstbl['setup'] = 'A'

    ## Hacky but needed to workaround PypeIt
    # pypeit.PypeIt crashes if all ra/dec data is missing, so we fake some
    if sum(ps.fitstbl['ra'] != None) == 0:
        ps.fitstbl['ra'][0] = 0
    if sum(ps.fitstbl['dec'] != None) == 0:
        ps.fitstbl['dec'][0] = 0

    ofiles = ps.fitstbl.write_pypeit(configs='A', cfg_lines=CFG_LINES)

    pypeIt = pypeit.PypeIt(ofiles[0], verbosity=0,
                           reuse_masters=True, overwrite=True,
                           logname='dbsp_ql.log', show=False, calib_only=calib_only)
    if calib_only:
        pypeIt.calib_all()
    else:
        pypeIt.reduce_all()
    pypeIt.build_qa()

    output_spec2ds = list(filter(lambda f: os.path.isfile(os.path.join('Science', f)), [
            pypeIt.spec_output_file(i, True) \
            for i in range(len(pypeIt.fitstbl.table)) \
            if pypeIt.fitstbl.table[i]['frametype'] in ['science']
        ]))

    output_spec1ds = list(filter(lambda f: os.path.isfile(os.path.join('Science', f)), [
            pypeIt.spec_output_file(i) \
            for i in range(len(pypeIt.fitstbl.table)) \
            if pypeIt.fitstbl.table[i]['frametype'] in ['science']
        ]))

    if output_spec1ds and not calib_only:
        spec = ps.spectrograph
        config = '_'.join([
            arm,
            spec.get_meta_value(sciimg[0], 'dispname').replace('/', '_'),
            spec.get_meta_value(sciimg[0], 'dichroic').lower()
        ])
        sensfiles = [resource_filename("dbsp_drp", f"data/sens_{config}.fits")]
        FxCalib = fluxcalibrate.FluxCalibrate.get_instance(output_spec1ds, sensfiles, par=ps.par['fluxcalib'])

    print(f"Time elapsed: {time.perf_counter() - t}s.")

    if not calib_only and not args.no_show:
        p1 = Process(target = show_spec2d_helper, args=(output_spec2ds[0],))
        p1.start()
        if output_spec1ds:
            with fits.open(output_spec1ds[0]) as hdul:
                specs = len(hdul) - 2
            parr = [ None ] * specs
            for i in range(specs):
                parr[i] = Process(target = show_spec1d_helper,
                    args=(str(i+1), output_spec1ds[0]))
                parr[i].start()

def show_spec2d_helper(file):
    return show_2dspec.Show2DSpec.main(show_2dspec.Show2DSpec.parse_args([file]))

def show_spec1d_helper(exten, file):
    return show_spectrum.main(
        show_spectrum.parser(['--BOX', '--exten', exten, file])
    )
