"""
Automated reduction steps for P200 DBSP.
"""

import os
import shutil
from typing import Tuple, List

import matplotlib.pyplot as plt

from astropy.io import fits

from configobj import ConfigObj

from pypeit.pypeitsetup import PypeItSetup
import pypeit
import pypeit.pypeit
from pypeit import msgs
from pypeit.spec2dobj import Spec2DObj
from pypeit.specobjs import SpecObjs

from dbsp_drp.manual_tracing import ManualTracingGUI

def parse_pypeit_parameter_file(args: dict) -> None:
    user_config_lines = []
    read_lines = False
    arm = 'red' if 'red' in args['spectrograph'] else 'blue'
    other_arm = 'red' if 'blue' in arm else 'blue'

    if os.path.isfile(args['parameter_file']):
        with open(args['parameter_file']) as f:
            for line in f.readlines():
                if f'[{other_arm}]' in line:
                    read_lines = False
                if read_lines:
                    user_config_lines.append(line)
                if f'[{arm}]' in line:
                    read_lines = True
    args['user_config_lines'] = user_config_lines

def setup(args: dict) -> Tuple[PypeItSetup, str]:
    """Does PypeIt setup, without writing the .pypeit file

    Args:
        args (dict): [description]

    Raises:
        ValueError: [description]
        ValueError: [description]
        IOError: [description]

    Returns:
        Tuple[PypeItSetup, str]: [description]
    """

    # Get the output directory
    output_path = os.getcwd() if args['output_path'] is None else args['output_path']
    sort_dir = os.path.join(output_path, 'setup_files')

    # Initialize PypeItSetup based on the arguments
    if args['root'] is not None:
        ps = PypeItSetup.from_file_root(args['root'], args['spectrograph'],
                                        extension=args['extension'], output_path=sort_dir)
    else:
        # Should never reach here
        raise IOError('Need to set -r !!')

    # Run the setup
    ps.run(setup_only=True, sort_dir=sort_dir, write_bkg_pairs=args['background'])

    return (ps, output_path)

def write_setup(args: dict, context: Tuple[PypeItSetup, str]) -> List[str]:
    """
    Writes the .pypeit file
    """
    ps, output_path = context
    # Use PypeItMetaData to write the complete PypeIt file
    config_list = [item.strip() for item in args['cfg_split'].split(',')]

    ps.user_cfg.append('[calibrations]')
    ps.user_cfg.append('master_dir = Master_' + args['spectrograph'].split('_')[-1])

    user_configobj = ConfigObj(ps.user_cfg)
    user_configobj.merge(ConfigObj(args['user_config_lines']))
    ps.user_cfg = [line + "\n" for line in user_configobj.write()]

    return ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg,
                                   write_bkg_pairs=args['background'], configs=config_list)

def redux(args: dict) -> None:
    """
    Runs the reduction
    """
    splitnm = os.path.splitext(args['pypeit_file'])
    if splitnm[1] != '.pypeit':
        msgs.error("Bad extension for PypeIt reduction file."+msgs.newline()+".pypeit is required")
    logname = splitnm[0] + ".log"

    pypeIt = pypeit.pypeit.PypeIt(args['pypeit_file'], verbosity=2,
                           reuse_masters=not args['do_not_reuse_masters'],
                           overwrite=True, #args['overwrite'],
                           redux_path=args['output_path'], #args['redux_path'], # rethink these
                           calib_only=args['calib_only'],
                           logname=logname, show=args['show'])
    pypeIt.reduce_all()
    msgs.info('Data reduction complete')

    msgs.info('Generating QA HTML')
    pypeIt.build_qa()

    args['output_spec1ds'] |= set(filter(lambda f: os.path.isfile(os.path.join(args['output_path'], 'Science', f)), [
            os.path.basename(pypeIt.spec_output_file(i)) \
            for i in range(len(pypeIt.fitstbl.table)) \
            if pypeIt.fitstbl.table[i]['frametype'] in ['science', 'standard']
        ]))

    args['output_spec2ds'] |= set(filter(lambda f: os.path.isfile(os.path.join(args['output_path'], 'Science', f)), [
            os.path.basename(pypeIt.spec_output_file(i, True)) \
            for i in range(len(pypeIt.fitstbl.table)) \
            if pypeIt.fitstbl.table[i]['frametype'] in ['science', 'standard']
        ]))

def delete_duplicate_hdus_by_name(path: str, base_name: str = ""):
    """
    Removes `SpecObj`s with identical names, leaving one behind.
    """
    specobjs = SpecObjs.from_fitsfile(path)
    if len({sobj['NAME'] for sobj in specobjs}) < specobjs.nobj:
        print(f'{base_name} has duplicate traces, deleting them')
        names = set()
        for i, sobj in enumerate(specobjs):
            if sobj['NAME'] in names:
                specobjs.remove_sobj(i)
            else:
                names.add(sobj['NAME'])
        specobjs.write_to_fits(specobjs.header, path, overwrite=True)
        specobjs.write_info(os.path.splitext(path)[0] + '.txt', "MultiSlit")


def verify_spec1ds(args: dict) -> List[str]:
    """
    Verifies validity of spec1d files, fixes some, and generates and returns a
    list of pypeit files for targets that need to be rerun.
    """
    # TODO: have different reasons files can be flagged for re-reduction, with
    #   a corresponding set of parameters to change
    # TODO: have redux produce and this function consume a set
    #   args['unverified_spec1ds'] so only changed files are re-checked
    # TODO: figure out params to tweak to fix skipped traces
    targets_list = []
    for spec1d in args['output_spec1ds']:
        path = os.path.join(args['output_path'], 'Science', spec1d)
        base_name = spec1d.split('_')[1]

        delete_duplicate_hdus_by_name(path, base_name)

        with fits.open(path) as hdul:
            for hdu in hdul:
                if 'SPAT' in hdu.name:
                    names = hdu.data.dtype.names
                    if ('TRACE_SPAT' in names) and not ('OPT_WAVE' in names) and not ('BOX_WAVE' in names):
                        # only TRACE_SPAT exists, we need to re-run this!
                        # need to investigate how this should be fixed
                        # for now, warn the user and recommend manually tracing
                        print(f'WARNING: trace {hdu.name} in {base_name} was not extracted!')
                        print('It is recommended to re-reduce this object with the -m flag'
                            'and manually place the desired traces.')
    return []

def manual_extraction_GUI(args):
    arm = 'blue' if 'blue' in args['spectrograph'] else 'red'
    spec2ds = args['output_spec2ds']

    gui_dict = {}
    for spec2d in spec2ds:
        path = os.path.join(args['output_path'], 'Science', spec2d)
        # open fits file
        spec = Spec2DObj.from_file(path, 1)
        spec1d_file = spec2d.replace('spec2d', 'spec1d', 1)
        spec1d_file = os.path.join(args['output_path'], 'Science', spec1d_file)

        if os.path.isfile(spec1d_file):
            spec1ds = SpecObjs.from_fitsfile(spec1d_file)
        else:
            spec1ds = None

        base_name = os.path.basename(path).split('_')[1]

        all_left, all_right, _ = spec.slits.select_edges()
        edges = [all_left, all_right]
        traces = [spec1ds[i]['TRACE_SPAT'] for i in range(len(spec1ds))] if spec1ds is not None else None
        fwhms = [spec1ds[i]['FWHM'] for i in range(len(spec1ds))] if spec1ds is not None else None

        gui_dict[base_name] = {
            'spec': spec,
            'edges': edges,
            'traces': traces,
            'fwhms': fwhms
        }

    # call GUI
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gui = ManualTracingGUI(fig, ax, gui_dict)

    manual_traces = gui.manual_dict

    return manual_traces

def write_manual_pypeit_files(args: dict) -> List[str]:
    """
    Writes pypeit files based on the default pypeit file for this reduction,
    but with filtered targets and additional manual parameter lines.

    requires
    args = {
        'pypeit_file': default pypeit file,
        'targets_list': [[target1, target2], [target3, target4]] will result in
            1 & 2 being reduced together and 3 & 4 being reduced together. The
            target names are blueNNNN-ZTF21abcd.
        'manual_lines_fn': function mapping [target1, target2] to the cfg lines
            they need
        'needs_std_fn': function target1 -> bool, True if target1 needs a
            standard star to be reduced alongside it
    }
    """
    old_pypeit_file = args['pypeit_file']

    targets_list = args['targets_list']
    manual_lines_fn = args['manual_lines_fn']
    needs_std_fn = args['needs_std_fn']

    new_pypeit_files = []
    for targets in targets_list:
        if not targets:
            continue
        target_fnames = [target.split('-')[0] for target in targets]
        new_pypeit_file = f'{os.path.splitext(old_pypeit_file)[0]}_{"_".join(targets)}.pypeit'
        shutil.copy(old_pypeit_file, new_pypeit_file)

        manual_lines = manual_lines_fn(targets)

        cfg_lines = []
        setup_lines = []
        setup = False
        with open(old_pypeit_file, 'r') as old_pypeit_fd:
            for line in old_pypeit_fd.readlines():
                if 'science' in line and '|' in line and all([targ_fname not in line for targ_fname in target_fnames]):
                    pass
                elif 'standard' in line and '|' in line and any(needs_std_fn(target) for target in targets):
                    setup_lines.append(line)
                else:
                    if '# Setup' in line:
                        setup = True
                    if setup:
                        setup_lines.append(line)
                    else:
                        cfg_lines.append(line)

        # Note: the user's parameter file is merged into these config lines
        # so any manual settings user gives might override their GUI actions.
        # Should probably caution user re: this.
        final_cfg = ConfigObj(manual_lines)
        final_cfg.merge(ConfigObj(cfg_lines))
        final_lines = [line + "\n" for line in final_cfg.write()] + setup_lines
        with open(new_pypeit_file, 'w') as new_pypeit_fd:
            new_pypeit_fd.writelines(final_lines)

        new_pypeit_files.append(new_pypeit_file)
    return new_pypeit_files


def manual_extraction(args: dict) -> list:
    manual_dict = manual_extraction_GUI(args)

    args['targets_list'] = [[key] for key in manual_dict.keys()]
    args['manual_lines_fn'] = lambda targ: ['# Added by DBSP_DRP for manual extraction\n',
            '[reduce]\n',
            '[[extraction]]\n',
            '[[[manual]]]\n',
            f"spat_spec = {str(manual_dict[targ[0]]['spat_spec']).strip('[]')}\n",
            f"det = {str([1 for trace in manual_dict[targ[0]]['spat_spec']]).strip('[]')}\n",
            f"fwhm = {str(manual_dict[targ[0]]['fwhm']).strip('[]')}\n",
            "[[skysub]]\n",
            f"user_regions = {str(manual_dict[targ[0]]['bgs']).strip('[]')}\n"
        ]
    args['needs_std_fn'] = lambda targ: manual_dict[targ]['needs_std']
    ret = write_manual_pypeit_files(args)
    # cleanup args so that multiprocessing can pickle it
    del args['targets_list']
    del args['manual_lines_fn']
    del args['needs_std_fn']
    return ret

def re_redux(args: dict, pypeit_files: list) -> None:
    for pypeit_file in pypeit_files:
        these_opts = args.copy()
        these_opts['pypeit_file'] = pypeit_file
        print(f"Using pypeit file {pypeit_file}")
        redux(these_opts)
