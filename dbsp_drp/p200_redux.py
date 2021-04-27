"""
Automatic Reduction Pipeline for P200 DBSP.
"""

import argparse
import os
import time
import multiprocessing
from typing import Optional, List
import pickle
import glob

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table, Column, Row

from pypeit.pypeitsetup import PypeItSetup
import pypeit.display

import tqdm

from dbsp_drp import reduction, qa, fluxing, coadding, telluric, splicing
from dbsp_drp import table_edit
from dbsp_drp import fix_headers
from dbsp_drp import instruments


def parser(options: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses command line arguments

    Args:
        options (Optional[List[str]], optional): List of command line arguments.
            Defaults to sys.argv[1:].

    Returns:
        argparse.Namespace: Parsed arguments
    """
    # Define command line arguments.
    argparser = argparse.ArgumentParser(description="Automatic Data Reduction Pipeline for P200 DBSP",
        formatter_class=argparse.RawTextHelpFormatter)

    # Argument for fully-automatic (i.e. nightly) or with user-checking file typing
    argparser.add_argument('-i', '--no-interactive', default=False, action='store_true',
                           help='Interactive file-checking?')

    # Argument for input file directory
    argparser.add_argument('-r', '--root', type=str, default=None,
                           help='File path+root, e.g. /data/DBSP_20200127')

    argparser.add_argument('-d', '--output_path', default=None,
                           help='Path to top-level output directory.  '
                                'Default is the current working directory.')

    # Argument for specifying only red/blue

    argparser.add_argument('-a', '--arm', default=None,
                           help='Space-separated list of arms to reduce. '
                                'Leave blank to reduce all. '
                                'For DBSP use `red` and `blue`. '
                                'For NGPS use `u` `g` `r` `i`')

    argparser.add_argument('--instrument', default=None,
                           help='Use this to manually set the instrument. '
                                'dbsp for dbsp, ngps for ngps. Usually not '
                                'needed.')

    argparser.add_argument('-m', '--manual-extraction', default=False, action='store_true',
                           help='manual extraction')

    argparser.add_argument('--debug', default=False, action='store_true',
                           help='debug')
    argparser.add_argument('-j', '--jobs', type=int, default=1, metavar='N',
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

    :param ps: PypeIt metadata object created in dbsp_drp.reduction.setup
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

    instrument: instruments.Instrument = None
    if args.instrument is None:
        # infer spectrograph from raw data path
        raw_data = glob.glob(os.path.join(args.root, "*.fits"))
        instrument = instruments.guess_instrument_from_file(raw_data[0])
    else:
        for ins in instruments.instruments:
            if ins.__name__.lower() == args.instrument.lower():
                instrument = ins()
        if instrument is None:
            raise NotImplementedError(f"Instrument {args.instrument} has not yet been implemented")


    do_arms: List[bool]
    if args.arm:
        do_arms = [False] * instrument.arms
        arms = args.arm.lower().split()
        for i, arm in enumerate(instrument.arm_prefixes):
            if arm in arms:
                do_arms[i] = True
    else:
        do_arms = [True] * instrument.arms

    roots = [os.path.join(args.root, arm) for arm in instrument.arm_prefixes]
    qa_dict = {}

    user_config_lines = [reduction.parse_pypeit_parameter_file(args.parameter_file, arm, instrument.arm_prefixes) for arm in instrument.arm_prefixes]

    if args.debug:
        pypeit.display.display.connect_to_ginga(raise_err=True, allow_new=True)

    pypeit_files = [''] * instrument.arms
    for i in range(instrument.arms):
        if do_arms[i]:
            fix_headers.main(roots[i])
            context = reduction.setup(roots[i], '.fits', args.output_path, instrument.arm_names_pypeit[i])
            # optionally use interactive correction
            if not args.no_interactive:
                interactive_correction(context[0])
            pypeit_files[i] = reduction.write_setup(context, 'all', instrument.arm_names_pypeit[i], user_config_lines[i])[0]

    plt.switch_backend("agg")
    # TODO: parallelize this
    # Would need to look like
    # Splitting up the .pypeit files into bits and pieces
    # Oooh what if I just do the calibration first
    # and then parallelize the reduction
    output_spec1ds = [set()] * instrument.arms
    output_spec2ds = [set()] * instrument.arms
    for i in range(instrument.arms):
        if do_arms[i]:
            output_spec1ds[i], output_spec2ds[i] = reduction.redux(pypeit_files[i], args.output_path)
            qa_dict = qa.save_2dspecs(qa_dict, output_spec2ds[i], args.output_path)

    qa.write_extraction_QA(qa_dict, args.output_path, type(instrument).__name__)

    for i in range(instrument.arms):
        if do_arms[i]:
            verification_counter = 0
            tmp_pypeit_files = reduction.verify_spec1ds(output_spec1ds[i], verification_counter, args.output_path)
            while tmp_pypeit_files:
                verification_counter += 1

                out_1d, out_2d = reduction.re_redux(tmp_pypeit_files, args.output_path)
                tmp_pypeit_files = reduction.verify_spec1ds(out_1d, verification_counter, args.output_path)
                qa_dict = qa.save_2dspecs(qa_dict, out_2d, args.output_path)

                output_spec1ds[i] |= out_1d
                output_spec2ds[i] |= out_2d

    # TODO: use a do/while loop to iterate on the manual extraction GUI until user is satisfied
    if args.manual_extraction:
        # wait for user acknowledgement
        input("Ready for manual extraction? If using GNU screen/tmux behind ssh, make sure to check that $DISPLAY is correct.")
        plt.switch_backend("Qt5Agg")

        manual_pypeit_files = [''] * instrument.arms
        for i in range(instrument.arms):
            if do_arms[i]:
                manual_pypeit_files[i] = reduction.manual_extraction(output_spec2ds[i], pypeit_files[i], args.output_path)
        for i in range(instrument.arms):
            if do_arms[i] and manual_pypeit_files[i]:
                out_1d, out_2d = reduction.re_redux(manual_pypeit_files[i], args.output_path)
                qa.save_2dspecs(qa_dict, out_2d, args.output_path)

                output_spec1ds[i] |= out_1d
                output_spec2ds[i] |= out_2d

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
    spec1ds = set.union(*output_spec1ds)
    for spec1d in spec1ds:
        path = os.path.join(args.output_path, 'Science', spec1d)
        with fits.open(path) as hdul:
            head0 = hdul[0].header
            head1 = hdul[1].header
            arm = instrument.pypeit_name_to_arm[head0['PYP_SPEC']]
            spec1d_table.add_row((spec1d, arm, head0['TARGET'],
                head1['OBJTYPE'], head0['AIRMASS'],
                head0['MJD'], '', head0['EXPTIME']))
    spec1d_table.add_index('filename')
    spec1d_table.sort(['arm', 'mjd'])

    for i, arm in enumerate(instrument.arm_prefixes):
        if do_arms[i]:
            for row in spec1d_table[(spec1d_table['arm'] == arm) * (spec1d_table['frametype'] == 'standard')]:
                sensfunc = fluxing.make_sensfunc(row['filename'], args.output_path, instrument.arm_names_pypeit[i], user_config_lines[i])
                if sensfunc:
                    spec1d_table.loc[row['filename']]['sensfunc'] = sensfunc
                else:
                    spec1d_table.loc[row['filename']]['frametype'] = 'science'

    for i, arm in enumerate(instrument.arm_prefixes):
        if do_arms[i]:
            arm_mask = spec1d_table['arm'] == arm
            stds = (spec1d_table['frametype'] == 'standard') & arm_mask
            if np.any(stds):
                for row in spec1d_table[arm_mask]:
                    if row['frametype'] == 'science':
                        best_sens = spec1d_table[stds]['sensfunc'][np.abs(spec1d_table[stds]['airmass'] - row['airmass']).argmin()]
                    elif row['frametype'] == 'standard':
                        if stds.sum() == 1:
                            # if there is only 1 standard, it already has itself as its own sensfunc!
                            best_sens = row['sensfunc']
                        else:
                            best_sens = spec1d_table[stds]['sensfunc'][np.abs(spec1d_table[stds]['airmass'] - row['airmass']).argsort()[1]]
                    spec1d_table.loc[row['filename']]['sensfunc'] = best_sens
            else:
                for filename in spec1d_table[arm]['filename']:
                    spec1d_table.loc[filename]['sensfunc'] = ''

    # build fluxfile
    fluxfiles = [''] * instrument.arms
    for i, arm in enumerate(instrument.arm_prefixes):
        if do_arms[i]:
            spec1d_to_sensfunc = {row['filename']: row['sensfunc'] for row in spec1d_table if row['arm'] == arm}
            fluxfiles[i] = fluxing.build_fluxfile(spec1d_to_sensfunc,
                instrument.archived_sensfuncs[i], args.output_path,
                instrument.arm_names_pypeit[i], user_config_lines[i])

    # flux data
    for i in range(instrument.arms):
        if do_arms[i]:
            fluxing.flux(fluxfiles[i], args.output_path)

    # coadd - intelligent coadding of multiple files
    # first make a column "coaddID" that is the same for frames to be coadded
    # TODO: when there are multiple exposures of an object, splice/output all of them
    coaddIDs = []
    if args.null_coadd:
        coaddIDs = range(len(spec1d_table))
    else:
        previous_row : Row = None
        S_PER_DAY = 24 * 60 * 60
        thresh = instrument.coadd_threshholds
        for i, row in enumerate(spec1d_table):
            if i == 0:
                coaddIDs.append(0)
            else:
                # if this is the same object as the last one
                # and they were taken consecutively
                if ((row['arm'] == previous_row['arm']) and
                    (row['object'] == previous_row['object']) and
                    ((row['mjd']*S_PER_DAY - previous_row['mjd']*S_PER_DAY
                        - previous_row['exptime']) < thresh[row['arm']])):
                    coaddIDs.append(coaddIDs[-1])
                else:
                    coaddIDs.append(coaddIDs[-1] + 1)
            previous_row = row

    spec1d_table.add_column(coaddIDs, name="coadd_id")

    # figure out where on detector likely target is
    spec1d_table.add_column(Column(name="spats", dtype=object, length=len(spec1d_table)))
    spec1d_table.add_column(Column(name="fracpos", dtype=object, length=len(spec1d_table)))
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
            spec1d_table.loc[filename]['spats'] = spats
            spec1d_table.loc[filename]['fracpos'] = fracpos
    # add to table???
    # this needs to be dtype object to allow for variable length lists
    spec1d_table.add_column(Column(name="coadds", dtype=object, length=len(spec1d_table)))
    spec1d_table.add_column([False]*len(all_spats), name="processed")

    # coadd
    # iterate over coadd_ids
    coadd_to_spec1d = {}
    for coadd_id in set(coaddIDs):
        subtable = spec1d_table[spec1d_table['coadd_id'] == coadd_id]
        fname_spats = {row['filename']: row['spats'].copy() for row in subtable}
        grouped_spats_list = coadding.group_coadds(fname_spats)
        for i, arm in enumerate(instrument.arm_prefixes):
            if all(subtable['arm'] == arm):
                coadds = coadding.coadd(grouped_spats_list, args.output_path, instrument.arm_names_pypeit[i], user_config_lines[i])
        assert any([all(subtable['arm'] == arm) for arm in instrument.arm_prefixes]),\
            "Something went wrong with coadding: spec1ds from multiple arms have the same coadd ID!"
        for row in subtable:
            spec1d_table.loc[row['filename']]['coadds'] = coadds
        for i, coadd in enumerate(coadds):
            coadd_to_spec1d[coadd] = list(zip(grouped_spats_list[i]['fnames'], grouped_spats_list[i]['spats']))

    if not args.skip_telluric:
        # telluric correct
        tellcorr_inputs = []
        tell_coadd_fnames = set()
        for i, arm in enumerate(instrument.arm_prefixes):
            if do_arms[i] and instrument.arm_telluric[i]:
                for row in spec1d_table[spec1d_table['arm'] == arm]:
                    if isinstance(row['coadds'], list):
                        for obj in row['coadds']:
                            if not obj in tell_coadd_fnames:
                                tmp = (obj, args.output_path, instrument.arm_names_pypeit[i], user_config_lines[i])
                                tellcorr_inputs.append(tmp)
                                tell_coadd_fnames.add(obj)
        if tellcorr_inputs:
            if args.jobs == 1:
                # do it in series
                for tellcorr_input in tqdm.tqdm(tellcorr_inputs):
                    telluric.telluric_correct(*tellcorr_input)
            else:
                pool = multiprocessing.Pool(args.jobs)
                list(tqdm.tqdm(pool.imap(telluric.picklable_telluric_correct, tellcorr_inputs), total=len(tellcorr_inputs)))
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

    os.makedirs(os.path.join(args.output_path, 'spliced'))

    ## Need to find red + blue fracpos for standards
    # hopefully standards only have one star each?
    # or should i actually try to do matching there
    fracpos_diff_list = []
    stds = spec1d_table['frametype'] == 'standard'

    if sum(do_arms) > 1:
        tol = instrument.calibrate_trace_matching(spec1d_table, args.output_path)

    # when we splice, we are splicing coadds
    # so we maybe need a coadd_table
    # with summed exposure times
    # and mean fracpos
    # and then splicing dict
    # will have to take into account
    # the observation time
    # so that multiple observations throughout the night
    # are kept separate

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
            if sum(do_arms) > 1:
                fracpos = instrument.convert_fracpos(arm, fracpos)
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
                    fracpos_diff_list.append(abs(fracpos_existing - fracpos))
                    if abs(fracpos_existing - fracpos) < tol:
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
    splicing.splice(splicing_dict, roots[0], args.output_path, instrument)

    with open("fracpos_data.pickle", "wb") as f:
        pickle.dump((fracpos_diff_list, instrument.FRACPOS_SUM), f)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))
