"""
Automatic Reduction Pipeline for P200 DBSP.
"""

import argparse
import os
import time
import multiprocessing
from typing import Optional, List
import pickle

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table, Column, Row

from pypeit.pypeitsetup import PypeItSetup
import pypeit.display
from pypeit.spectrographs.util import load_spectrograph

import tqdm

from dbsp_drp import reduction, qa, fluxing, coadding, telluric, splicing
from dbsp_drp import table_edit
from dbsp_drp import fix_headers


def entrypoint():
    main(parser())

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
    argparser.add_argument('-r', '--root', type=os.path.abspath, default=None,
                           required=True,
                           help='File path+root, e.g. /data/DBSP_20200127')

    argparser.add_argument('-d', '--output_path', type=os.path.abspath,
                           default='.',
                           help='Path to top-level output directory.  '
                                'Default is the current working directory.')

    # Argument for specifying only red/blue

    argparser.add_argument('-a', '--arm', default=None, choices=['red', 'blue'],
                           help='[red, blue] to only reduce one arm (null splicing)')

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

    argparser.add_argument('--splicing-interpolate-gaps', default=False, action='store_true',
                           help="Use this option to linearly interpolate across large gaps\n"
                                "in the spectrum during splicing. The default behavior is to\n"
                                "only use data from one detector in these gaps, which results\n"
                                "in a slightly noisier spliced spectrum.")

    return argparser.parse_args() if options is None else argparser.parse_args(options)

def interactive_correction(ps: PypeItSetup) -> None:
    """
    Allows for human correction of FITS headers and frame typing.

    Launches a GUI via dbsp_drp.table_edit, which handles saving updated FITS headers.
    table_edit depends on the current DBSP headers.

    Todo:
        Make table to FITS header mapping mutable

    Args:
        ps (PypeItSetup): PypeItSetup object created in dbsp_drp.reduction.setup
    """
    # function for interactively correcting the fits table
    fitstbl = ps.fitstbl
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

    red_root = os.path.join(args.root, 'red')
    blue_root = os.path.join(args.root, 'blue')
    qa_dict = {}

    if args.parameter_file:
        blue_user_config_lines = reduction.parse_pypeit_parameter_file(args.parameter_file, 'p200_dbsp_blue')
        red_user_config_lines = reduction.parse_pypeit_parameter_file(args.parameter_file, 'p200_dbsp_red')
    else:
        blue_user_config_lines = []
        red_user_config_lines = []

    if args.debug:
        pypeit.display.display.connect_to_ginga(raise_err=True, allow_new=True)

    if do_red:
        red_files = fix_headers.main(red_root, args.no_interactive, args.no_interactive)
        context = reduction.setup(red_files, args.output_path, 'p200_dbsp_red')
        # optionally use interactive correction
        if not args.no_interactive:
            interactive_correction(context[0])
        pypeit_file_red = reduction.write_setup(context, 'all', 'p200_dbsp_red', red_user_config_lines)[0]

    if do_blue:
        blue_files = fix_headers.main(blue_root, args.no_interactive, args.no_interactive)
        context = reduction.setup(blue_files, args.output_path, 'p200_dbsp_blue')
        if not args.no_interactive:
            interactive_correction(context[0])
        pypeit_file_blue = reduction.write_setup(context, 'all', 'p200_dbsp_blue', blue_user_config_lines)[0]


    plt.switch_backend("agg")
    # TODO: parallelize this
    # Would need to look like
    # Splitting up the .pypeit files into bits and pieces
    # Oooh what if I just do the calibration first
    # and then parallelize the reduction
    output_spec1ds_blue = set()
    output_spec1ds_red = set()
    if do_red:
        output_spec1ds_red, output_spec2ds_red = reduction.redux(pypeit_file_red, args.output_path)
        qa_dict = qa.save_2dspecs(qa_dict, output_spec2ds_red, args.output_path, 'p200_dbsp_red')
    if do_blue:
        output_spec1ds_blue, output_spec2ds_blue = reduction.redux(pypeit_file_blue, args.output_path)
        qa_dict = qa.save_2dspecs(qa_dict, output_spec2ds_blue, args.output_path, 'p200_dbsp_blue')

    if do_red or do_blue:
        qa.write_extraction_QA(qa_dict, args.output_path)

    if do_red:
        verification_counter = 0
        red_pypeit_files = reduction.verify_spec1ds(output_spec1ds_red, verification_counter, args.output_path)
        while red_pypeit_files:
            verification_counter += 1

            out_1d, out_2d = reduction.re_redux(red_pypeit_files, args.output_path)
            red_pypeit_files = reduction.verify_spec1ds(out_1d, verification_counter, args.output_path)
            qa_dict = qa.save_2dspecs(qa_dict, out_2d, args.output_path, 'p200_dbsp_red')

            output_spec1ds_red |= out_1d
            output_spec2ds_red |= out_2d
    if do_blue:
        verification_counter = 0
        blue_pypeit_files = reduction.verify_spec1ds(output_spec1ds_blue, verification_counter, args.output_path)
        while blue_pypeit_files:
            verification_counter += 1

            out_1d, out_2d = reduction.re_redux(blue_pypeit_files, args.output_path)
            blue_pypeit_files = reduction.verify_spec1ds(out_1d, verification_counter, args.output_path)
            qa_dict = qa.save_2dspecs(qa_dict, out_2d, args.output_path, 'p200_dbsp_blue')

            output_spec1ds_blue |= out_1d
            output_spec2ds_blue |= out_2d

    # TODO: use a do/while loop to iterate on the manual extraction GUI until user is satisfied
    if args.manual_extraction:
        # wait for user acknowledgement
        input("Ready for manual extraction? If using GNU screen/tmux behind ssh, make sure to check that $DISPLAY is correct.")
        plt.switch_backend("Qt5Agg")

        if do_red:
            red_manual_pypeit_files = reduction.manual_extraction(output_spec2ds_red, pypeit_file_red, args.output_path)
        if do_blue:
            blue_manual_pypeit_files = reduction.manual_extraction(output_spec2ds_blue, pypeit_file_blue, args.output_path)
        if do_red and red_manual_pypeit_files:
            out_1d, out_2d = reduction.re_redux(red_manual_pypeit_files, args.output_path)
            qa.save_2dspecs(qa_dict, out_2d, args.output_path, 'p200_dbsp_red')

            output_spec1ds_red |= out_1d
            output_spec2ds_red |= out_2d
        if do_blue and blue_manual_pypeit_files:
            out_1d, out_2d = reduction.re_redux(blue_manual_pypeit_files, args.output_path)
            qa.save_2dspecs(qa_dict, out_2d, args.output_path, 'p200_dbsp_blue')

            output_spec1ds_blue |= out_1d
            output_spec2ds_blue |= out_2d

    # Find standards and make sensitivity functions
    spec1d_table = Table(names=('filename', 'arm', 'object', 'frametype',
                            'airmass', 'mjd', 'sensfunc', 'exptime'),
                         dtype=(f'U255', 'U4', 'U255', 'U8',
                            float, float, f'U255', float))

    # Ingest spec_1d tables
    spec1ds = output_spec1ds_red | output_spec1ds_blue
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
            sensfunc = fluxing.make_sensfunc(row['filename'], args.output_path, 'p200_dbsp_red', red_user_config_lines)
            if sensfunc == "":
                spec1d_table['frametype'][spec1d_table['filename'] == row['filename']] = 'science'
            else:
                spec1d_table['sensfunc'][spec1d_table['filename'] == row['filename']] = sensfunc
    if do_blue:
        for row in spec1d_table[(spec1d_table['arm'] == 'blue') & (spec1d_table['frametype'] == 'standard')]:
            sensfunc = fluxing.make_sensfunc(row['filename'], args.output_path, 'p200_dbsp_blue', blue_user_config_lines)
            if sensfunc == "":
                spec1d_table['frametype'][spec1d_table['filename'] == row['filename']] = 'science'
            else:
                spec1d_table['sensfunc'][spec1d_table['filename'] == row['filename']] = sensfunc

    if do_red:
        arm = spec1d_table['arm'] == 'red'
        stds = (spec1d_table['frametype'] == 'standard') & arm

        red_arm = load_spectrograph('p200_dbsp_red')
        rawfile = os.path.join(args.root,
            spec1d_table[arm][0]['filename'].split('_')[1].split('-')[0] + '.fits'
        )
        config = '_'.join([
            'red',
            red_arm.get_meta_value(rawfile, 'dispname').replace('/', '_'),
            red_arm.get_meta_value(rawfile, 'dichroic').lower()
        ])
        if np.any(stds):
            for row in spec1d_table[arm]:
                if row['frametype'] == 'science':
                    best_sens = spec1d_table[stds]['sensfunc'][np.abs(spec1d_table[stds]['airmass'] - row['airmass']).argmin()]
                elif row['frametype'] == 'standard':
                    if (stds).sum() == 1:
                        best_sens = spec1d_table[stds]['sensfunc'][np.abs(spec1d_table[stds]['airmass'] - row['airmass']).argmin()]
                    else:
                        best_sens = spec1d_table[stds]['sensfunc'][np.abs(spec1d_table[stds]['airmass'] - row['airmass']).argsort()[1]]
                spec1d_table.loc[row['filename']]['sensfunc'] = best_sens
        else:
            for filename in spec1d_table[arm]['filename']:
                spec1d_table.loc[filename]['sensfunc'] = config
    if do_blue:
        arm = spec1d_table['arm'] == 'blue'
        stds = (spec1d_table['frametype'] == 'standard') & arm

        blue_arm = load_spectrograph('p200_dbsp_blue')
        rawfile = os.path.join(args.root,
            spec1d_table[arm][0]['filename'].split('_')[1].split('-')[0] + '.fits'
        )
        config = '_'.join([
            'blue',
            blue_arm.get_meta_value(rawfile, 'dispname').replace('/', '_'),
            blue_arm.get_meta_value(rawfile, 'dichroic').lower()
        ])
        if np.any(stds):
            for row in spec1d_table[arm]:
                if row['frametype'] == 'science':
                    best_sens = spec1d_table[stds]['sensfunc'][np.abs(spec1d_table[stds]['airmass'] - row['airmass']).argmin()]
                elif row['frametype'] == 'standard':
                    if (stds).sum() == 1:
                        best_sens = spec1d_table[stds]['sensfunc'][np.abs(spec1d_table[stds]['airmass'] - row['airmass']).argmin()]
                    else:
                        best_sens = spec1d_table[stds]['sensfunc'][np.abs(spec1d_table[stds]['airmass'] - row['airmass']).argsort()[1]]
                spec1d_table.loc[row['filename']]['sensfunc'] = best_sens
        else:
            for filename in spec1d_table[arm]['filename']:
                spec1d_table.loc[filename]['sensfunc'] = config

    # build fluxfile
    if do_red:
        spec1d_to_sensfunc = {row['filename']: row['sensfunc'] for row in spec1d_table if row['arm'] == 'red'}
        red_fluxfile = fluxing.build_fluxfile(spec1d_to_sensfunc, args.output_path, 'p200_dbsp_red', red_user_config_lines)
    if do_blue:
        spec1d_to_sensfunc = {row['filename']: row['sensfunc'] for row in spec1d_table if row['arm'] == 'blue'}
        blue_fluxfile = fluxing.build_fluxfile(spec1d_to_sensfunc, args.output_path, 'p200_dbsp_blue', blue_user_config_lines)

    # flux data
    if do_red:
        fluxing.flux(red_fluxfile, args.output_path)
    if do_blue:
        fluxing.flux(blue_fluxfile, args.output_path)

    # coadd - intelligent coadding of multiple files
    # first make a column "coaddID" that is the same for frames to be coadded
    # TODO: when there are multiple exposures of an object, splice/output all of them
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
                if ((row['arm'] == previous_row['arm']) and
                    (row['object'] == previous_row['object']) and
                    ((row['mjd']*S_PER_DAY - previous_row['mjd']*S_PER_DAY
                        - previous_row['exptime']) < previous_row['exptime'])):
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
        if all(subtable['arm'] == 'red'):
            coadds = coadding.coadd(grouped_spats_list, args.output_path, 'p200_dbsp_red', red_user_config_lines)
        if all(subtable['arm'] == 'blue'):
            coadds = coadding.coadd(grouped_spats_list, args.output_path, 'p200_dbsp_blue', blue_user_config_lines)
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
                            tmp = (obj, args.output_path, 'p200_dbsp_red', red_user_config_lines)
                            tellcorr_inputs.append(tmp)
                            tell_coadd_fnames.add(obj)
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
    splicing_dict = {}
    blue_mask = spec1d_table['arm'] == 'blue'
    red_mask = spec1d_table['arm'] == 'red'

    os.makedirs(os.path.join(args.output_path, 'spliced'), exist_ok=True)

    def get_std_trace(std_path: str) -> float:
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

    ## Need to find red + blue fracpos for standards
    # hopefully standards only have one star each?
    # or should i actually try to do matching there
    stds = spec1d_table['frametype'] == 'standard'
    if do_red or do_blue:
        FRACPOS_SUM = 1.0
        FRACPOS_TOL = 0.05
        if do_red and do_blue:
            # real matching + splicing
            std_fracpos_sums = []
            if (stds & blue_mask).any() and (stds & red_mask).any():
                for row in spec1d_table[stds]:
                    # find closest mjd frame of other arm
                    if not row['processed']:
                        other_arm = spec1d_table['arm'] != row['arm']
                        corresponding_row = spec1d_table[other_arm & stds][np.abs(spec1d_table[other_arm & stds]['mjd'] - row['mjd']).argmin()]
                        this_path = os.path.join(args.output_path, 'Science', row['filename'])
                        corresponding_path = os.path.join(args.output_path, 'Science', corresponding_row['filename'])
                        std_fracpos_sums.append(get_std_trace(this_path) + get_std_trace(corresponding_path))
                        spec1d_table.loc[row['filename']]['processed'] = True
                        spec1d_table.loc[corresponding_row['filename']]['processed'] = True
                FRACPOS_SUM = np.mean(std_fracpos_sums)
                FRACPOS_TOL = FRACPOS_SUM * .025

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
        splicing.splice(splicing_dict, args.splicing_interpolate_gaps, red_root, args.output_path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))
