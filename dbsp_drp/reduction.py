"""
Automated reduction steps for P200 DBSP.
"""

import datetime
import os
import shutil
from typing import Tuple, List, Callable

import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table, Row
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u

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
    """
    Grab user-provided PypeIt configuration for ``spectrograph`` from
    ``parameter_file``, which contains user-provided PypeIt configuration for
    both arms of DBSP.

    Args:
        parameter_file (str): User-provided file with their PypeIt config.
        spectrograph (str): PypeIt name of spectrograph.

    Returns:
        List[str]: User-provided PypeIt configuration for ``spectrograph``
    """
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

def search_table_for_arc(row: Row, i: int, table: Table, step: int, max_sep: Angle) -> Tuple[int, int]:
    """
    Searches table starting at row ``i`` in steps of ``step`` for frames within
    ``max_separation`` of ``row``. If an arc frame is found, the arc's ``calib``
    is assigned to ``row``.

    Args:
        row (Row): starting row
        i (int): index of starting row
        table (Table): table to search in and modify
        step (int): direction to look in. 1 for forwards, -1 for backwards
        max_sep (Angle): maximum separation allowed between frames pointing
            at the same object

    Returns:
        Tuple[int, int]: (distance to arc frame, calib ID) or (-1, -1) if
            no arc found.
    """
    i_coord = SkyCoord(row['ra'], row['dec'], unit=u.deg)
    j = step
    keepLooking = True
    while keepLooking:
        if (i + j < 0) or (i + j >= len(table)):
            keepLooking = False
            continue
        j_coord = SkyCoord(table[i+j]['ra'], table[i+j]['dec'], unit=u.deg)
        sep = i_coord.separation(j_coord)
        if sep < max_sep:
            if 'arc' in table[i+j]['frametype']:
                table[i]['calib'] = table[i+j]['calib']
                return (abs(j), table[i]['calib'])
            j += step
        else:
            keepLooking = False
    return (-1, -1)

def set_calibs(table: Table):
    """
    Automagically set 'calib' column for .pypeit file.

    Bias and flat frames get calib 'all'
    Arcs at airmass 1.0 get calib 0 (default for science/standards)
    A consecutive set of arcs at airmass > 1.0 get the same calib, starting at 1.
    Science and standards get the calib of the nearest arc frame, iff the
    telescope hasn't changed pointing between the arc and the science/standard.
    If there is no arc frame with the same pointing, science and standards are
    assigned calib 0.

    Args:
        table (Table): PypeItSetup.fitstbl.table
    """
    # bias and flats get calib 'all'
    for row in table:
        if 'bias' in row['frametype'] or 'trace' in row['frametype']:
            table.loc[row['filename']]['calib'] = 'all'
        if 'arc' in row['frametype'] and row['airmass'] == 1.0:
            table.loc[row['filename']]['calib'] = '0'

    calib_ID = 1
    max_separation = Angle(20 * u.arcsec)
    prev_row = None
    prev_coord = None
    for row in table:
        if ('arc' in row['frametype']) and (row['airmass'] != 1.0) and (row['ra'] is not None) and (row['dec'] is not None):
            coord = SkyCoord(row['ra'], row['dec'], unit=u.deg)
            if calib_ID == 1:
                table.loc[row['filename']]['calib'] = calib_ID
                calib_ID += 1
            else:
                sep = coord.separation(prev_coord)
                if ((sep < max_separation) and
                    ('arc' in prev_row['frametype'])):
                    table.loc[row['filename']]['calib'] = calib_ID - 1
                else:
                    table.loc[row['filename']]['calib'] = calib_ID
                    calib_ID += 1
            prev_coord = coord
        prev_row = row

    # assign calib_IDs to science / standard objects
    # matching to nearest arc by number / filename such that the telescope pointing didn't move
    # between the arc and the science / standard

    for i, row in enumerate(table):
        # or instead of matching, just check neighbors and see if neighbors match
        if ('science' in row['frametype']) or ('standard' in row['frametype']):
            dist_forward, calib_forward = search_table_for_arc(row, i, table, 1, max_separation)
            dist_back, _ = search_table_for_arc(row, i, table, -1, max_separation)
            if (dist_forward == -1) and (dist_back == -1):
                table[i]['calib'] = 0
            elif ((dist_forward != -1) and (dist_back != -1) and
                (dist_forward < dist_back)):
                    table[i]['calib'] = calib_forward

def setup(file_list: List[str], output_path: str, spectrograph: str) -> Tuple[PypeItSetup, str]:
    """
    Does PypeIt setup, without writing the .pypeit file

    Args:
        file_list (List[str]): List of raw data files to reduce.
        output_path (str): reduction output path
        spectrograph (str): PypeIt name of spectrograph.

    Returns:
        Tuple[PypeItSetup, str]: PypeItSetup object, reduction output path
    """
    # TODO: do we need to check if output_path is None or return it?
    # Get the output directory
    output_path = os.getcwd() if output_path is None else output_path
    sort_dir = os.path.join(output_path, 'setup_files')
    os.makedirs(sort_dir, exist_ok=True)

    # Initialize PypeItSetup based on the arguments
    cfg_lines = ['[rdx]', f'spectrograph = {spectrograph}']
    fname = os.path.join(sort_dir, f'{spectrograph}_{datetime.date.today().strftime("%Y-%m-%d")}.pypeit')
    ps = PypeItSetup(file_list, setups=[], cfg_lines=cfg_lines, pypeit_file=fname)

    # Run the setup
    ps.run(setup_only=True, sort_dir=sort_dir)


    table = ps.fitstbl.table
    table.sort('filename')
    table.add_index('filename')

    # Make the target field accept long names
    table['target'] = table['target'].astype('U255')

    # now we guess the calib ids
    set_calibs(table)

    return (ps, output_path)

def write_setup(context: Tuple[PypeItSetup, str], cfg_split: str,
        spectrograph: str, user_config_lines: List[str]) -> List[str]:
    """
    Writes the .pypeit file

    Args:
        context (Tuple[PypeItSetup, str]): PypeItSetup object, reduction output path
        cfg_split (str): [description]
        spectrograph (str): PypeIt name of spectrograph.
        user_config_lines (List[str]): User-provided PypeIt configuration.

    Raises:
        RuntimeError: Raised if files were not assigned a frame type.

    Returns:
        List[str]: List of PypeIt files generated.
    """
    ps, output_path = context

    if 'None' in ps.fitstbl['frametype']:
        untyped_files = list(ps.fitstbl['filename'][ps.fitstbl['frametype'] == 'None'])
        err_msg = "File " if len(untyped_files) == 1 else "Files "
        err_msg += f"{str(untyped_files).strip('[]')} were not automatically assigned a frame type, please re-run without -i argument and ensure that all files have frametypes."
        raise RuntimeError(err_msg)

    # remove comb_id and bkg_id columns because they are not exposed to
    # DBSP_DRP user, so PypeIt's defaults are good enough AFTER frametypes are
    # corrected.
    ps.fitstbl.table.remove_columns(['comb_id', 'bkg_id'])
    # Use PypeItMetaData to write the complete PypeIt file
    config_list = [item.strip() for item in cfg_split.split(',')]

    ps.user_cfg.append('[calibrations]')
    ps.user_cfg.append('master_dir = Master_' + spectrograph.split('_')[-1])

    user_configobj = ConfigObj(ps.user_cfg)
    user_configobj.merge(ConfigObj(user_config_lines))
    ps.user_cfg = [line + "\n" for line in user_configobj.write()]

    return ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg, configs=config_list, write_bkg_pairs=True)

def redux(pypeit_file: str, output_path: str, reuse_masters: bool = True,
        show: bool = False, calib_only: bool = False) -> Tuple[set, set]:
    """
    Runs the reduction

    Args:
        pypeit_file (str): Path to PypeIt reduction file.
        output_path (str): reduction output path
        reuse_masters (bool, optional): Reuse master calibration files (if they exist). Defaults to True.
        show (bool, optional): Show debugging/intermediate plots. Defaults to False.
        calib_only (bool, optional): Only perform calibration? Defaults to False.

    Returns:
        Tuple[set, set]: set of filenames of reduced (spec1d, spec2d) files.
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
    Removes ``SpecObj`` s with identical names from a FITS file containing a
    ``SpecObjs`` object, leaving one ``SpecObj`` with the duplicate name behind.

    Args:
        path (str): Path to FITS file.
        base_name (str, optional): Name of file used for logging. Defaults to "".
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

def delete_completely_masked_hdus(path: str, base_name: str = ""):
    """
    Removes ``SpecObj`` s that are completely masked, i.e. no good pixels/data,
    from a FITS file containing a ``SpecObjs`` file.

    Args:
        path (str): Path to FITS file.
        base_name (str, optional): Name of file used for logging. Defaults to "".
    """
    specobjs = SpecObjs.from_fitsfile(path)
    changed = False
    for i, specobj in enumerate(specobjs):
        if ((specobj.OPT_MASK is not None and not specobj.OPT_MASK.any()) or
            (specobj.BOX_MASK is not None and not specobj.BOX_MASK.any())):
            print(f'{base_name} has a completely masked trace, deleting that trace')
            specobjs.remove_sobj(i)
            changed = True

    if changed:
        specobjs.write_to_fits(specobjs.header, path, overwrite=True)
        specobjs.write_info(os.path.splitext(path)[0] + '.txt', "MultiSlit")

def verify_spec1ds(output_spec1ds: List[str], verification_counter: int,
        output_path: str) -> List[str]:
    """
    Verifies validity of spec1d files, fixes some, and generates and returns a
    list of pypeit files for targets that need to be rerun.

    Args:
        output_spec1ds (List[str]): List of spec1d filenames
        verification_counter (int): number of times verification has been run
        output_path (str): reduction output path

    Returns:
        List[str]: List of PypeIt files for targets that need to be rereduced.
    """
    # TODO: have different reasons files can be flagged for re-reduction, with
    #   a corresponding set of parameters to change
    # TODO: have redux produce and this function consume a set
    #   of `unverified_spec1ds` so only changed files are re-checked
    targets_list = []
    for spec1d in output_spec1ds:
        path = os.path.join(output_path, 'Science', spec1d)
        base_name = spec1d.split('_')[1]

        delete_duplicate_hdus_by_name(path, base_name)
        delete_completely_masked_hdus(path, base_name)

        with fits.open(path) as hdul:
            for hdu in hdul:
                if 'SPAT' in hdu.name:
                    names = hdu.data.dtype.names
                    if ('TRACE_SPAT' in names) and not ('OPT_WAVE' in names) and not ('BOX_WAVE' in names):
                        # only TRACE_SPAT exists, object was traced but not extracted
                        # we need to re-run this!
                        # The underlying bug in PypeIt has been fixed as of v1.4.1
                        # but this may still occur so
                        # warn the user and recommend manually tracing
                        print(f'WARNING: trace {hdu.name} in {base_name} was not extracted!')
                        print('It is recommended to re-reduce this object with the -m flag'
                            'and manually place the desired traces.')
    return []

def manual_extraction_GUI(output_spec2ds: List[str], output_path: str) -> dict:
    """
    Runs the Manual Tracing GUI.

    Args:
        output_spec2ds (List[str]): List of spec2d files to inspect in GUI.
        output_path (str): reduction output path

    Returns:
        dict: Maps target names needing manual tracing to a dict specifying how
            they should be manually traced.
    """
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
        fwhmfit = [spec1ds[i]['FWHMFIT'] for i in range(len(spec1ds))] if spec1ds is not None else None

        gui_dict[base_name] = {
            'spec': spec,
            'edges': edges,
            'traces': traces,
            'fwhms': fwhms,
            'fwhmfit': fwhmfit
        }

    # call GUI
    gui = ManualTracingGUI(gui_dict)

    manual_traces = gui.manual_dict

    return manual_traces

def write_manual_pypeit_files(old_pypeit_file: str, targets_list: List[List[str]],
        manual_lines_fn: Callable[[List[str]], List[str]],
        needs_std_fn: Callable[[str], bool]) -> List[str]:
    """
    Writes pypeit files based on the default pypeit file for this reduction,
    but with filtered targets and additional manual parameter lines.

    Args:
        old_pypeit_file (str): default pypeit file,
        targets_list (List[List[str]]): [[target1, target2], [target3, target4]]
            will result in 1 & 2 being reduced together and 3 & 4 being reduced
            together. The target names are blueNNNN-ZTF21abcd.
        manual_lines_fn (Callable[[List[str]], List[str]]): function mapping
            [target1, target2] to the cfg lines they need.
        needs_std_fn (Callable[[str], bool]): function mapping target1 -> bool.
            True if target1 needs a standard star to be reduced alongside it.

    Returns:
        List[str]: List of new PypeIt reduction files for the manually traced
            targets.
    """

    new_pypeit_files = []
    for targets in targets_list:
        if not targets:
            continue
        target_fnames = [target.split('-')[0] for target in targets]
        new_pypeit_file = f'{os.path.splitext(old_pypeit_file)[0]}_{"_".join(targets)}.pypeit'
        shutil.copyfile(old_pypeit_file, new_pypeit_file)

        manual_lines = manual_lines_fn(targets)

        cfg_lines = []
        setup_lines = []
        setup = False
        with open(old_pypeit_file, 'r') as old_pypeit_fd:
            for line in old_pypeit_fd.readlines():
                if 'science' in line and '|' in line and all([targ_fname not in line for targ_fname in target_fnames]):
                    pass
                elif 'standard' in line and '|' in line and any(needs_std_fn(target) for target in targets):
                    pass
                    #setup_lines.append(line)
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
    """
    wglwk

    Args:
        output_spec2ds (List[str]): List of spec2d files to potentially
            manually extract.
        pypeit_file (str): PypeIt reduction file for initial reduction.
        output_path (str): reduction output path.

    Returns:
        list: List of new PypeIt reduction files for the manually traced targets.
    """
    manual_dict = manual_extraction_GUI(output_spec2ds, output_path)

    targets_list = [[key] for key in manual_dict.keys()]
    manual_lines_fn = lambda targ: ['# Added by DBSP_DRP for manual extraction\n',
            '[reduce]\n',
            '[[extraction]]\n',
            'use_user_fwhm = True\n',
            '[[[manual]]]\n',
            f"spat_spec = {str(manual_dict[targ[0]]['spat_spec']).strip('[]')}\n",
            f"det = {str([1 for trace in manual_dict[targ[0]]['spat_spec']]).strip('[]')}\n",
            f"fwhm = {str(manual_dict[targ[0]]['fwhm']).strip('[]')}\n",
            "[[skysub]]\n",
            f"user_regions = {str(manual_dict[targ[0]]['bgs']).strip('[]')}\n"
        ]
    needs_std_fn = lambda targ: manual_dict[targ]['needs_std']
    return write_manual_pypeit_files(pypeit_file, targets_list, manual_lines_fn, needs_std_fn)

def re_redux(pypeit_files: List[str], output_path: str) -> Tuple[set, set]:
    """
    Runs multiple reductions, returns the combined sets of output spec1d,
    spec2d files.

    Args:
        pypeit_files (List[str]): List of PypeIt reduction files to run.
        output_path (str): reduction output path

    Returns:
        Tuple[set, set]: set of filenames of reduced (spec1d, spec2d) files.
    """
    output_spec1ds = set()
    output_spec2ds = set()
    for pypeit_file in pypeit_files:
        print(f"Using pypeit file {pypeit_file}")
        out_1d, out_2d = redux(pypeit_file, output_path)
        output_spec1ds |= out_1d
        output_spec2ds |= out_2d
    return output_spec1ds, output_spec2ds
