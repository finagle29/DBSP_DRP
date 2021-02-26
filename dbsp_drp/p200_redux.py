"""
Automatic Reduction Pipeline for P200 DBSP.
"""

import argparse
import os
import time
import copy
import multiprocessing
from typing import Optional, List

from astropy.io import fits
from astropy.table import Table, Column, Row
import numpy as np
from pypeit.pypeitsetup import PypeItSetup
import pypeit.display
import tqdm

import matplotlib as mpl
mpl.use("Qt5Agg")
from matplotlib import pyplot as plt

from dbsp_drp import p200_arm_redux
from dbsp_drp import table_edit
from dbsp_drp import fix_headers


def parser(options: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses command line arguments

    Args:
        options (Optional[List[str]], optional): List of command line arguments. Defaults to sys.argv[1:].

    Returns:
        argparse.Namespace: Parsed arguments
    """
    # Define command line arguments.
    argparser = argparse.ArgumentParser(description="Automatic Data Reduction Pipeline for P200 DBSP",
        formatter_class=argparse.RawTextHelpFormatter)

    # Argument for fully-automatic (i.e. nightly) or with user-checking file typing
    argparser.add_argument('-i', '--no-interactive', default=True, action='store_false',
                           help='Interactive file-checking?')

    # Argument for input file directory
    argparser.add_argument('-r', '--root', type=str, default=None,
                           help='File path+root, e.g. /data/DBSP_20200127')

    argparser.add_argument('-d', '--output_path', default=None,
                           help='Path to top-level output directory.  '
                                'Default is the current working directory.')

    # Argument for specifying only red/blue

    argparser.add_argument('-a', '--arm', default=None,
                           help='[red, blue] to only reduce one arm (null splicing)')

    argparser.add_argument('-m', '--manual-extraction', default=False, action='store_true',
                           help='manual extraction')

    argparser.add_argument('--debug', default=False, action='store_true',
                           help='debug')
    argparser.add_argument('-j', '--jobs', type=int, default=1,
                            help='Number of processes to use')

    argparser.add_argument('-p', '--parameter-file', type=str, default="",
                           help="Path to parameter file. The parameter file should be formatted as follows:\n\n"
                            "[blue]\n"
                            "** PypeIt parameters for the blue side goes here **\n"
                            "[red]\n"
                            "** PypeIt parameters for the red side goes here **\n"
                            "EOF\n\n"
                            "The [red/blue] parameter blocks are optional, and their order does not matter.")

    argparser.add_argument('-t', '--skip-telluric', default=False, action='store_true',
                           help='Skip telluric correction')

    argparser.add_argument('-c', '--null-coadd', default=False, action='store_true',
                           help="Don't coadd consecutive exposures of the same target.\n"
                                "By default consective exposures will be coadded.")

    return argparser.parse_args() if options is None else argparser.parse_args(options)

def interactive_correction(ps: PypeItSetup) -> None:
    """Allows for human correction of FITS headers and frame typing.

    Launches a GUI via dbsp_drp.table_edit, which handles saving updated FITS headers.
    table_edit depends on the current DBSP headers.

    Todo:
        Make table to FITS header mapping mutable

    :param ps: PypeIt metadata object created in p200_arm_redux.setup
    :type ps: PypeItSetup
    """
    # function for interactively correcting the fits table
    fitstbl = ps.fitstbl
    fitstbl.table.sort('filename')
    deleted_files = []
    table_edit.main(fitstbl.table, deleted_files)
    files_to_remove = []
    for rm_file in deleted_files:
        for data_file in ps.file_list:
            if rm_file in data_file:
                files_to_remove.append(data_file)
                break
    for rm_file in files_to_remove:
        ps.file_list.remove(rm_file)

def main(args):

    t = time.perf_counter()

    if os.path.isdir(args.output_path):
        os.chdir(args.output_path)
    else:
        os.makedirs(args.output_path, exist_ok=True)

    if args.arm:
        do_red = args.arm.lower() == 'red'
        do_blue = args.arm.lower() == 'blue'
    else:
        do_red = True
        do_blue = True

    # Build
    options_blue = {
        'root': os.path.join(args.root, 'blue'),
        'spectrograph': 'p200_dbsp_blue',
        'output_path': args.output_path,
        'extension': '.fits',
        'background': None,
        'cfg_split': 'all',
        'calib_only': False,
        'show': args.debug,
        'plot': args.debug,
        'do_not_reuse_masters': False,
        'debug': args.debug,
        'qa_dict': {},
        'manual_extraction': args.manual_extraction,
        'output_spec1ds': set(),
        'output_spec2ds': set(),
        'parameter_file': args.parameter_file
    }

    #options_blue['show'] = True
    #options_blue['calib_only'] = True
    #options_blue['do_not_reuse_masters'] = True

    options_red = copy.deepcopy(options_blue)
    options_red['spectrograph'] = 'p200_dbsp_red'
    options_red['root'] = os.path.join(args.root, 'red')
    options_red['qa_dict'] = options_blue['qa_dict']

    p200_arm_redux.parse_pypeit_parameter_file(options_blue)
    p200_arm_redux.parse_pypeit_parameter_file(options_red)

    if args.debug:
        pypeit.display.display.connect_to_ginga(raise_err=True, allow_new=True)

    if do_red:
        fix_headers.main(options_red['root'])
        context = p200_arm_redux.setup(options_red)
        # optionally use interactive correction
        if args.no_interactive:
            interactive_correction(context[0])
        pypeit_file_red = p200_arm_redux.write_setup(options_red, context)[0]
        # Grab pypeit file from write_setup
        options_red['pypeit_file'] = pypeit_file_red

    if do_blue:
        fix_headers.main(options_blue['root'])
        context = p200_arm_redux.setup(options_blue)
        if args.no_interactive:
            interactive_correction(context[0])
        pypeit_file_blue = p200_arm_redux.write_setup(options_blue, context)[0]
        # Grab pypeit file from write_setup
        options_blue['pypeit_file'] = pypeit_file_blue


    #options_red['calib_only'] = True
    #options_blue['calib_only'] = True
    plt.switch_backend("agg")
    if do_red:
        p200_arm_redux.redux(options_red)
        p200_arm_redux.save_2dspecs(options_red)
    if do_blue:
        p200_arm_redux.redux(options_blue)
        p200_arm_redux.save_2dspecs(options_blue)

    if do_red or do_blue:
        p200_arm_redux.write_extraction_QA(options_red)

    # TODO: use a do/while loop to iterate on the manual extraction GUI until user is satisfied
    if args.manual_extraction:
        # wait for user acknowledgement
        input("Ready for manual extraction? If using GNU screen/tmux behind ssh, make sure to check that $DISPLAY is correct.")

        if do_red:
            red_manual_pypeit_files = p200_arm_redux.manual_extraction(options_red)
        if do_blue:
            blue_manual_pypeit_files = p200_arm_redux.manual_extraction(options_blue)
        if do_red:
            p200_arm_redux.re_redux(options_red, red_manual_pypeit_files)
            p200_arm_redux.save_2dspecs(options_red)
        if do_blue:
            p200_arm_redux.re_redux(options_blue, blue_manual_pypeit_files)
            p200_arm_redux.save_2dspecs(options_blue)

    # spec1d_blueNNNN-OBJ_DBSPb_YYYYMMMDDTHHMMSS.SPAT.fits
    fname_len = 72
    # sens_blueNNNN-OBJ_DBSPb_YYYYMMMDDTHHMMSS.SPAT.fits
    sensfunc_len = 70
    # Find standards and make sensitivity functions
    spec1d_table = Table(names=('filename', 'arm', 'object', 'frametype',
                            'airmass', 'mjd', 'sensfunc', 'exptime'),
                         dtype=(f'U{fname_len}', 'U4', 'U20', 'U8',
                            float, float, f'U{sensfunc_len}', float))

    # Ingest spec_1d tables
    spec1ds = options_red['output_spec1ds'] | options_blue['output_spec1ds']
    for spec1d in spec1ds:
        path = os.path.join(args.output_path, 'Science', spec1d)
        with fits.open(path) as hdul:
            head0 = hdul[0].header
            head1 = hdul[1].header
            arm = 'red' if 'red' in head0['PYP_SPEC'] else 'blue'
            spec1d_table.add_row((spec1d, arm, head0['TARGET'],
                head1['OBJTYPE'], head0['AIRMASS'],
                head0['MJD'], '', head0['EXPTIME']))
    spec1d_table.add_index('filename')
    spec1d_table.sort(['arm', 'mjd'])

    if do_red:
        for row in spec1d_table[(spec1d_table['arm'] == 'red') & (spec1d_table['frametype'] == 'standard')]:
            options_red['spec1dfile'] = row['filename']
            sensfunc = p200_arm_redux.make_sensfunc(options_red)
            if sensfunc == "":
                spec1d_table['frametype'][spec1d_table['filename'] == row['filename']] = 'science'
            else:
                spec1d_table['sensfunc'][spec1d_table['filename'] == row['filename']] = sensfunc
    if do_blue:
        for row in spec1d_table[(spec1d_table['arm'] == 'blue') & (spec1d_table['frametype'] == 'standard')]:
            options_blue['spec1dfile'] = row['filename']
            sensfunc = p200_arm_redux.make_sensfunc(options_blue)
            if sensfunc == "":
                spec1d_table['frametype'][spec1d_table['filename'] == row['filename']] = 'science'
            else:
                spec1d_table['sensfunc'][spec1d_table['filename'] == row['filename']] = sensfunc

    stds = spec1d_table['frametype'] == 'standard'
    standards_fluxing = []
    for row in spec1d_table:
        arm = spec1d_table['arm'] == row['arm']
        if row['frametype'] == 'science':
            best_sens = spec1d_table[stds & arm]['sensfunc'][np.abs(spec1d_table[stds & arm]['airmass'] - row['airmass']).argmin()]
            standards_fluxing.append(best_sens)
        elif row['frametype'] == 'standard':
            if (stds & arm).sum() == 1:
                best_sens = spec1d_table[stds & arm]['sensfunc'][np.abs(spec1d_table[stds & arm]['airmass'] - row['airmass']).argmin()]
                standards_fluxing.append(best_sens)
            else:
                best_sens = spec1d_table[stds & arm]['sensfunc'][np.abs(spec1d_table[stds & arm]['airmass'] - row['airmass']).argsort()[1]]
                standards_fluxing.append(best_sens)

    spec1d_table['sensfunc'] = standards_fluxing

    # build fluxfile
    if do_red:
        options_red['spec1dfiles'] = {row['filename']: row['sensfunc'] for row in spec1d_table if row['arm'] == 'red'}
        red_fluxfile = p200_arm_redux.build_fluxfile(options_red)
        options_red['flux_file'] = red_fluxfile
    if do_blue:
        options_blue['spec1dfiles'] = {row['filename']: row['sensfunc'] for row in spec1d_table if row['arm'] == 'blue'}
        blue_fluxfile = p200_arm_redux.build_fluxfile(options_blue)
        options_blue['flux_file'] = blue_fluxfile

    # flux data
    if do_red:
        p200_arm_redux.flux(options_red)
    if do_blue:
        p200_arm_redux.flux(options_blue)

    # coadd - intelligent coadding of multiple files
    # first make a column "coaddID" that is the same for frames to be coadded
    coaddIDs = []
    if args.null_coadd:
        coaddIDs = range(len(spec1d_table))
    else:
        previous_row : Row = None
        S_PER_DAY = 24 * 60 * 60
        thresh = 15
        for i, row in enumerate(spec1d_table):
            if i == 0:
                coaddIDs.append(0)
            else:
                # if this is the same object as the last one
                # and they were taken consecutively
                if ((row['object'] == previous_row['object']) and
                    ((row['mjd']*S_PER_DAY - previous_row['mjd']*S_PER_DAY
                        - previous_row['exptime']) < previous_row['exptime'])):
                    coaddIDs.append(coaddIDs[-1])
                else:
                    coaddIDs.append(coaddIDs[-1] + 1)
            previous_row = row

    spec1d_table.add_column(coaddIDs, name="coadd_id")

    #options_blue['debug'] = True
    #options_red['debug'] = True
    # figure out where on detector likely target is
    all_spats = []
    all_fracpos = []
    # for each spec1d file
    for filename in spec1d_table['filename']:
        path = os.path.join(args.output_path, 'Science', filename)
        with fits.open(path) as hdul:
            spats = []
            fracpos = []
            for i in range(1, len(hdul) - 1):
                # grab all of its extensions' spatial positions
                spats.append(int(hdul[i].name.split('-')[0].lstrip('SPAT')))
                fracpos.append(hdul[i].header['SPAT_FRACPOS'])
            spats.sort()
            fracpos.sort()
            all_spats.append(spats)
            all_fracpos.append(fracpos)
    # add to table???
    spec1d_table.add_column(all_spats, name="spats")
    spec1d_table.add_column(all_fracpos, name="fracpos")
    # this needs to be dtype object to allow for variable length lists
    spec1d_table.add_column(Column(name="coadds", dtype=object, length=len(spec1d_table)))
    spec1d_table.add_column([False]*len(all_spats), name="processed")

    ## Deprecated I think
    """
    # for each arm
        # find median/mean of spatial positions
    if do_blue:
        blue_spats = spec1d_table[spec1d_table['arm'] == 'blue']['spats']
        blue_spats = [spat for spats in blue_spats for spat in spats]
        avg_blue_spat = np.mean(blue_spats)
        std_blue_spats = np.std(blue_spats)
    if do_red:
        red_spats = spec1d_table[spec1d_table['arm'] == 'red']['spats']
        red_spats = [spat for spats in red_spats for spat in spats]
        avg_red_spat = np.mean(red_spats)
        std_red_spats = np.std(red_spats)
    """

    # coadd
    #options_red['debug'] = True
    #options_red['plot'] = True
    # iterate over coadd_ids
    coadd_to_spec1d = {}
    for coadd_id in set(coaddIDs):
        subtable = spec1d_table[spec1d_table['coadd_id'] == coadd_id]
        fname_spats = {row['filename']: row['spats'].copy() for row in subtable}
        grouped_spats_list = p200_arm_redux.group_coadds(fname_spats)
        if all(subtable['arm'] == 'red'):
            options_red['grouped_spats_list'] = grouped_spats_list
            coadds = p200_arm_redux.coadd(options_red)
        if all(subtable['arm'] == 'blue'):
            options_blue['grouped_spats_list'] = grouped_spats_list
            coadds = p200_arm_redux.coadd(options_blue)
        assert all(subtable['arm'] == 'red') or all(subtable['arm'] == 'blue'),\
            "Something went wrong with coadding..."
        for row in subtable:
            spec1d_table.loc[row['filename']]['coadds'] = coadds
        for i, coadd in enumerate(coadds):
            coadd_to_spec1d[coadd] = list(zip(grouped_spats_list[i]['fnames'], grouped_spats_list[i]['spats']))

    if not args.skip_telluric:
        # telluric correct
        if do_red:
            tellcorr_inputs = []
            tell_coadd_fnames = set()
            for row in spec1d_table[spec1d_table['arm'] == 'red']:
                if isinstance(row['coadds'], list):
                    for obj in row['coadds']:
                        if not obj in tell_coadd_fnames:
                            tmp = options_red.copy()
                            tmp['spec1dfile'] = obj
                            tellcorr_inputs.append(tmp)
                            tell_coadd_fnames.add(obj)
    #                    options_red['spec1dfile'] = obj
    #                    p200_arm_redux.telluric_correct(options_red)
            if args.jobs == 1:
                # do it in series
                for tellcorr_input in tqdm.tqdm(tellcorr_inputs):
                    p200_arm_redux.telluric_correct(tellcorr_input)
            else:
                pool = multiprocessing.Pool(args.jobs)
                list(tqdm.tqdm(pool.imap(p200_arm_redux.telluric_correct, tellcorr_inputs), total=len(tellcorr_inputs)))
            #for opts in tellcorr_inputs:
            #    pool.apply_async(p200_arm_redux.telluric_correct, opts, error_callback=lambda e: print("error!"))
                pool.close()
                pool.join()

            # Maybe do something here to verify that telluric correction succeeded
            # and if so, change the coadd names
            for coadd in tell_coadd_fnames:
                tell = coadd.replace(".fits", "_tellcorr.fits")
                tellpath = os.path.join(args.output_path, 'Science', tell)
                coaddpath = os.path.join(args.output_path, 'Science', coadd)
                # check if tell exists and is newer than coadd
                if os.path.isfile(tellpath) and (os.path.getmtime(tellpath) > os.path.getmtime(coaddpath)):
                    # modify coadd
                    for row in spec1d_table:
                        if coadd in row['coadds']:
                            ix = row['coadds'].index(coadd)
                            spec1d_table.loc[row['filename']]['coadds'][ix] = tell
                    coadd_to_spec1d[tell] = coadd_to_spec1d[coadd]
                    del coadd_to_spec1d[coadd]


    # current splicing - make sure spatial fraction is similar on blue/red
    # TODO: handle multiple observations of same target throughout night with null coadding
    # splice data
    splicing_dict = {}
    blue_mask = spec1d_table['arm'] == 'blue'
    red_mask = spec1d_table['arm'] == 'red'

    ## Need to find red + blue fracpos for standards
    # hopefully standards only have one star each?
    # or should i actually try to do matching there
    if do_red or do_blue:
        FRACPOS_TOL = 0.01
        if do_red and do_blue:
            # real matching + splicing
            std_fracpos_sums = []
            for row in spec1d_table[stds]:
                # find closest mjd frame of other arm
                if not row['processed']:
                    other_arm = spec1d_table['arm'] != row['arm']
                    corresponding_row = spec1d_table[other_arm & stds][np.abs(spec1d_table[other_arm & stds]['mjd'] - row['mjd']).argmin()]
                    std_fracpos_sums.append(row['fracpos'][0] + corresponding_row['fracpos'][0])
                    spec1d_table.loc[row['filename']]['processed'] = True
                    spec1d_table.loc[corresponding_row['filename']]['processed'] = True
            FRACPOS_SUM = np.mean(std_fracpos_sums)
            FRACPOS_TOL = FRACPOS_SUM * .01

        # setup splicing dict
        splicing_dict = {}
        # for each target
        for row in spec1d_table:
            target = row['object']
            arm = row['arm']
            # for each of its fracpos
            for i, fracpos in enumerate(row['fracpos']):
                coadd = row['coadds'][i]
                targ_dict = splicing_dict.get(target)
                # normalize fracpos to red
                if do_red and do_blue and arm == 'blue':
                    fracpos = FRACPOS_SUM - fracpos
                # if it's not in the dict
                if targ_dict is None:
                    # put it in the dict
                    splicing_dict[target] = {fracpos: {
                        arm: {
                            'spec1ds': coadd_to_spec1d[coadd],
                            'coadd': coadd
                        }
                    }}
                # else
                else:
                    close_enough = False
                    # for each existing fracpos
                    for fracpos_existing in list(targ_dict):
                        # if its close enough
                        if abs(fracpos_existing - fracpos) < FRACPOS_TOL:
                            # put it in the dict
                            splicing_dict[target][fracpos_existing][arm] = {
                                'spec1ds': coadd_to_spec1d[coadd],
                                'coadd': coadd
                            }
                            close_enough = True
                            break
                    if not close_enough:
                        # If this fracpos isn't close enough to any others
                        splicing_dict[target][fracpos] = {arm: {
                            'spec1ds': coadd_to_spec1d[coadd],
                            'coadd': coadd
                        }}
        # And now, actually splice!
        options_red['splicing_dict'] = splicing_dict
        p200_arm_redux.splice(options_red)

    """
    # if do_red or do_blue ????
    if do_red or do_blue:
        # make splicing dict
        for row in spec1d_table:
            target = row['object']
            spacings = np.array(row['spats'])
            if row['arm'] == 'blue':
                target_spat = avg_blue_spat
            else:
                target_spat = avg_red_spat
            best_ix = np.abs(spacings - target_spat).argmin()

            if row['arm'] == 'blue':
                best_spec = row['coadds'][best_ix]
            else:
                best_spec = row['coadds'][best_ix].replace(".fits", "_tellcorr.fits")
                if not os.path.isfile(best_spec):
                    # Telluric correction failed.
                    print(f"ERROR!! Telluric correction of spectrum {row['coadds'][best_ix]} FAILED!")
                    print(f"The final spectrum for {target} will NOT be telluric corrected.")
                    best_spec = row['coadds'][best_ix]


            if splicing_dict.get(target):
                splicing_dict[target][row['arm']] = best_spec
            else:
                splicing_dict[target] = {row['arm']: best_spec}
        # splice
        options_red['splicing_dict'] = splicing_dict
        p200_arm_redux.splice(options_red)
        """

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))
