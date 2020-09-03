"""
Automatic Reduction Pipeline for P200 DBSP.
"""

import argparse
import glob
import os
import multiprocessing

from astropy.io import fits
from astropy.table import Table, Column
import numpy as np

from dbsp_drp import p200_arm_redux
from dbsp_drp import table_edit


def parser(options=None):
    """
    Parses command line arguments
    """
    # Define command line arguments.
    argparser = argparse.ArgumentParser(description="Automatic Data Reduction Pipeline for P200 DBSP")

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
                           help='[red, blue] to only reduce one arm')

    argparser.add_argument('--debug', default=False, action='store_true',
                           help='debug')
    argparser.add_argument('-j', '--jobs', type=int, default=1,
                            help='Number of processes to use')

    return argparser.parse_args() if options is None else argparser.parse_args(options)

def interactive_correction(ps):
    # this needs to actually fix the data's FITS headers
    # deleting entire row from .pypeit file is valid though
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
        'show': False,
        'plot': False,
        'do_not_reuse_masters': False,
        'debug': args.debug,
        'qa_dict': {}
    }

    #options_blue['show'] = True
    #options_blue['calib_only'] = True
    #options_blue['do_not_reuse_masters'] = True

    options_red = options_blue.copy()
    options_red['spectrograph'] = 'p200_dbsp_red'
    options_red['root'] = os.path.join(args.root, 'red')
    options_red['qa_dict'] = options_blue['qa_dict']

    if do_red:
        context = p200_arm_redux.setup(options_red)
        # optionally use interactive correction
        if args.no_interactive:
            interactive_correction(context[0])
        pypeit_file_red = p200_arm_redux.write_setup(options_red, context)[0]
        # Grab pypeit file from write_setup
        options_red['pypeit_file'] = pypeit_file_red

    if do_blue:
        context = p200_arm_redux.setup(options_blue)
        if args.no_interactive:
            interactive_correction(context[0])
        pypeit_file_blue = p200_arm_redux.write_setup(options_blue, context)[0]
        # Grab pypeit file from write_setup
        options_blue['pypeit_file'] = pypeit_file_blue




    #options_red['calib_only'] = True
    #options_blue['calib_only'] = True
    if do_red:
        p200_arm_redux.redux(options_red)
        p200_arm_redux.save_2dspecs(options_red)
    if do_blue:
        p200_arm_redux.redux(options_blue)
        p200_arm_redux.save_2dspecs(options_blue)

    if do_red or do_blue:
        p200_arm_redux.write_extraction_QA(options_red)


    # Find standards and make sensitivity functions
    spec1d_table = Table(names=('filename', 'arm', 'object', 'frametype', 'airmass', 'mjd', 'sensfunc'),
                         dtype=('U100', 'U4', 'U20', 'U8', float, float, 'U100'))

    # Ingest spec_1d tables
    paths = glob.glob(os.path.join(args.output_path, "Science/spec1d*.fits"))
    for path in paths:
        with fits.open(path) as hdul:
            arm = 'red' if 'red' in path else 'blue'
            spec1d_table.add_row((path, arm, hdul[0].header['TARGET'], hdul[1].header['OBJTYPE'],
                                  hdul[0].header['AIRMASS'], hdul[0].header['MJD'], ''))
    spec1d_table.add_index('filename')

    if do_red:
        for row in spec1d_table[(spec1d_table['arm'] == 'red') & (spec1d_table['frametype'] == 'standard')]:
            options_red['spec1dfile'] = row['filename']
            spec1d_table['sensfunc'][spec1d_table['filename'] == row['filename']] = p200_arm_redux.make_sensfunc(options_red)
    if do_blue:
        for row in spec1d_table[(spec1d_table['arm'] == 'blue') & (spec1d_table['frametype'] == 'standard')]:
            options_blue['spec1dfile'] = row['filename']
            spec1d_table['sensfunc'][spec1d_table['filename'] == row['filename']] = p200_arm_redux.make_sensfunc(options_blue)

    stds = spec1d_table['frametype'] == 'standard'
    standards_fluxing = []
    for row in spec1d_table:
        arm = spec1d_table['arm'] == row['arm']
        if row['frametype'] == 'science':
            best_sens = spec1d_table[stds & arm]['sensfunc'][np.abs(spec1d_table[stds & arm]['airmass'] - row['airmass']).argmin()]
            standards_fluxing.append(best_sens)
        elif row['frametype'] == 'standard':
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

    # TODO: coadd - intelligent coadding of multiple files
    #options_blue['debug'] = True
    #options_red['debug'] = True
    # figure out where on detector likely target is
    all_spats = []
    # for each spec1d file
    for filename in spec1d_table['filename']:
        with fits.open(filename) as hdul:
            spats = []
            for i in range(1, len(hdul) - 1):
                # grab all of its extensions' spatial positions
                spats.append(int(hdul[i].name.split('-')[0].lstrip('SPAT')))
            all_spats.append(spats)
    # add to table???
    spec1d_table.add_column(all_spats, name="spats")
    spec1d_table.add_column(Column(name="coadds", dtype=object, length=len(spec1d_table))) # need to make this dtype object

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

    # coadd
    # probolem here
    if do_red:
        for row in spec1d_table[spec1d_table['arm'] == 'red']:
            options_red['spec1dfile'] = row['filename']
            spec1d_table.loc[row['filename']]['coadds'] = p200_arm_redux.coadd(options_red)
    if do_blue:
        for row in spec1d_table[spec1d_table['arm'] == 'blue']:
            options_blue['spec1dfile'] = row['filename']
            spec1d_table.loc[row['filename']]['coadds'] = p200_arm_redux.coadd(options_blue)

    # telluric correct
    options_red['debug'] = False
    options_red['plot'] = False
    if do_red:
        tellcorr_inputs = []
        for row in spec1d_table[spec1d_table['arm'] == 'red']:
            if isinstance(row['coadds'], list):
                for obj in row['coadds']:
                    tmp = options_red.copy()
                    tmp['spec1dfile'] = obj
                    tellcorr_inputs.append(tmp)
#                    options_red['spec1dfile'] = obj
#                    p200_arm_redux.telluric_correct(options_red)
        pool = multiprocessing.Pool(args.jobs)
        pool.map(p200_arm_redux.telluric_correct, tellcorr_inputs)
        #for opts in tellcorr_inputs:    
        #    pool.apply_async(p200_arm_redux.telluric_correct, opts, error_callback=lambda e: print("error!"))
        pool.close()
        pool.join()
    # splice data
    splicing_dict = {}
    blue_mask = spec1d_table['arm'] == 'blue'
    red_mask = spec1d_table['arm'] == 'red'
    if do_red and do_blue:
        # make splicing dict
        for target in spec1d_table['object']:
            pass
