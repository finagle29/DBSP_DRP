"""
Automated reduction steps for P200 DBSP.
"""

import os
import shutil
from typing import Tuple, List, Callable

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

def parse_pypeit_parameter_file(parameter_file: str,
        spectrograph: str) -> List[str]:
    user_config_lines = []
    read_lines = False
    arm = 'red' if 'red' in spectrograph else 'blue'
    other_arm = 'red' if 'blue' in arm else 'blue'

    if os.path.isfile(parameter_file):
        with open(parameter_file) as f:
            for line in f.readlines():
                if f'[{other_arm}]' in line:
                    read_lines = False
                if read_lines:
                    user_config_lines.append(line)
                if f'[{arm}]' in line:
                    read_lines = True
    else:
        print(f"ERROR: Parameter file {parameter_file} does not exist!!!")
    return user_config_lines

def setup(root: str, extension: str, output_path: str,
        spectrograph: str) -> Tuple[PypeItSetup, str]:
    """Does PypeIt setup, without writing the .pypeit file
    """

    # Get the output directory
    output_path = os.getcwd() if output_path is None else output_path
    sort_dir = os.path.join(output_path, 'setup_files')

    # Initialize PypeItSetup based on the arguments
    ps = PypeItSetup.from_file_root(root, spectrograph, extension=extension, output_path=sort_dir)

    # Run the setup
    ps.run(setup_only=True, sort_dir=sort_dir)

    return (ps, output_path)

def write_setup(context: Tuple[PypeItSetup, str], cfg_split: str,
        spectrograph: str, user_config_lines: List[str]) -> List[str]:
    """
    Writes the .pypeit file
    """
    ps, output_path = context
    # Use PypeItMetaData to write the complete PypeIt file
    config_list = [item.strip() for item in cfg_split.split(',')]

    ps.user_cfg.append('[calibrations]')
    ps.user_cfg.append('master_dir = Master_' + spectrograph.split('_')[-1])

    user_configobj = ConfigObj(ps.user_cfg)
    user_configobj.merge(ConfigObj(user_config_lines))
    ps.user_cfg = [line + "\n" for line in user_configobj.write()]

    return ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg, configs=config_list)

def redux(pypeit_file: str, output_path: str, reuse_masters: bool = True,
        show: bool = False, calib_only: bool = False) -> Tuple[set, set]:
    """
    Runs the reduction
    """
    splitnm = os.path.splitext(pypeit_file)
    if splitnm[1] != '.pypeit':
        msgs.error("Bad extension for PypeIt reduction file."+msgs.newline()+".pypeit is required")
    logname = splitnm[0] + ".log"

    pypeIt = pypeit.pypeit.PypeIt(pypeit_file, verbosity=2,
                           reuse_masters=reuse_masters,
                           overwrite=True,
                           redux_path=output_path, #args['redux_path'], # rethink these
                           calib_only=calib_only,
                           logname=logname, show=show)
    pypeIt.reduce_all()
    msgs.info('Data reduction complete')

    msgs.info('Generating QA HTML')
    pypeIt.build_qa()

    output_spec1ds = set(filter(lambda f: os.path.isfile(os.path.join(output_path, 'Science', f)), [
            os.path.basename(pypeIt.spec_output_file(i)) \
            for i in range(len(pypeIt.fitstbl.table)) \
            if pypeIt.fitstbl.table[i]['frametype'] in ['science', 'standard']
        ]))

    output_spec2ds = set(filter(lambda f: os.path.isfile(os.path.join(output_path, 'Science', f)), [
            os.path.basename(pypeIt.spec_output_file(i, True)) \
            for i in range(len(pypeIt.fitstbl.table)) \
            if pypeIt.fitstbl.table[i]['frametype'] in ['science', 'standard']
        ]))

    return output_spec1ds, output_spec2ds

def delete_duplicate_hdus_by_name(path: str, base_name: str = ""):
    """
    Removes ``SpecObj`` s with identical names, leaving one behind.
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


def verify_spec1ds(output_spec1ds: List[str], verification_counter: int, output_path: str) -> List[str]:
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
    for spec1d in output_spec1ds:
        path = os.path.join(output_path, 'Science', spec1d)
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

def manual_extraction_GUI(output_spec2ds: List[str], output_path: str) -> dict:
    spec2ds = output_spec2ds

    gui_dict = {}
    for spec2d in spec2ds:
        path = os.path.join(output_path, 'Science', spec2d)
        # open fits file
        spec = Spec2DObj.from_file(path, 1)
        spec1d_file = spec2d.replace('spec2d', 'spec1d', 1)
        spec1d_file = os.path.join(output_path, 'Science', spec1d_file)

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

def write_manual_pypeit_files(old_pypeit_file: str, targets_list: List[List[str]],
        manual_lines_fn: Callable[[List[str]], List[str]],
        needs_std_fn: Callable[[str], bool]) -> List[str]:
    """
    Writes pypeit files based on the default pypeit file for this reduction,
    but with filtered targets and additional manual parameter lines.

    Arguments:
        old_pypeit_file: default pypeit file,
        targets_list: [[target1, target2], [target3, target4]] will result in
            1 & 2 being reduced together and 3 & 4 being reduced together. The
            target names are blueNNNN-ZTF21abcd.
        manual_lines_fn: function mapping [target1, target2] to the cfg lines
            they need
        needs_std_fn: function target1 -> bool, True if target1 needs a
            standard star to be reduced alongside it
    """

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


def manual_extraction(output_spec2ds: List[str], pypeit_file: str,
        output_path: str) -> list:
    manual_dict = manual_extraction_GUI(output_spec2ds, output_path)

    targets_list = [[key] for key in manual_dict.keys()]
    manual_lines_fn = lambda targ: ['# Added by DBSP_DRP for manual extraction\n',
            '[reduce]\n',
            '[[extraction]]\n',
            '[[[manual]]]\n',
            f"spat_spec = {str(manual_dict[targ[0]]['spat_spec']).strip('[]')}\n",
            f"det = {str([1 for trace in manual_dict[targ[0]]['spat_spec']]).strip('[]')}\n",
            f"fwhm = {str(manual_dict[targ[0]]['fwhm']).strip('[]')}\n",
            "[[skysub]]\n",
            f"user_regions = {str(manual_dict[targ[0]]['bgs']).strip('[]')}\n"
        ]
    needs_std_fn = lambda targ: manual_dict[targ]['needs_std']
    return write_manual_pypeit_files(pypeit_file, targets_list, manual_lines_fn, needs_std_fn)

def re_redux(pypeit_files: list, output_path: str) -> Tuple[set, set]:
    output_spec1ds = set()
    output_spec2ds = set()
    for pypeit_file in pypeit_files:
        print(f"Using pypeit file {pypeit_file}")
        out_1d, out_2d = redux(pypeit_file, output_path)
        output_spec1ds |= out_1d
        output_spec2ds |= out_2d
    return output_spec1ds, output_spec2ds
