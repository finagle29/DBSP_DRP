"""
Automated coadding for P200 DBSP.
"""

import os
from typing import List

from astropy.io import fits

from pypeit import coadd1d
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import pypeitpar

def coadd(args: dict) -> List[str]:
    """
    takes in args['grouped_spats_list'], a list of dicts mapping 'fnames' to a
        list of filenames and 'spats' to a list of integer spatial pixel
        positions.
    Returns a list of filenames of coadded spectra.
    """
    outfiles = []
    for d in args['grouped_spats_list']:
        fnames = d['fnames']
        spats = d['spats']
        basename = '_'.join([fname.split("_")[1].split("-")[0]
            for fname in fnames]) + "_" + \
                fnames[0].split("_")[1].split("-")[1]

        objnames = []
        for spat, fname in zip(spats, fnames):
            path = os.path.join(args['output_path'], 'Science', fname)
            hdul = fits.open(path)
            for hdu in hdul:
                if f'SPAT{spat:04d}' in hdu.name:
                    objnames.append(hdu.name)
                    break

        outfile = os.path.join(args['output_path'], "Science", f"{basename}_{'_'.join(objnames)}.fits")
        coadd_one_object([os.path.join(args['output_path'], 'Science', fname) for fname in fnames],
            objnames, outfile, args)
        outfiles.append(os.path.basename(outfile))
    return outfiles

def coadd_one_object(spec1dfiles: List[str], objids: List[str], coaddfile: str, args: dict):
    par = load_spectrograph(args['spectrograph']).default_pypeit_par()
    default_cfg_lines = par.to_config()
    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines = default_cfg_lines, merge_with=args['user_config_lines'])
    # Instantiate
    coAdd1d = coadd1d.CoAdd1D.get_instance(spec1dfiles, objids, par=par['coadd1d'], debug=args['debug'], show=args['debug'])
    # Run
    coAdd1d.run()
    # Save to file
    coAdd1d.save(coaddfile)

def group_coadds(fname_to_spats: dict):
    """
    Groups coadds. Destroys input.

    Takes in dict mapping filenames to a list of integer spatial positions
    Returns list of dicts mapping 'fnames' to a list of filenames and 'spats'
        to a list of integer spatial positions.
    """
    # input is dict mapping fname to spats
    # end result is mapping from arb. label of trace -> spats list and fnames list
    THRESHOLD = 2
    result = []
    while any(fname_to_spats.values()):
        potential_group = [(spats[0], fname) for fname, spats in fname_to_spats.items()]
        potential_group.sort(key=lambda x: x[0])
        min_spat, its_fname = potential_group.pop()
        fname_to_spats[its_fname].remove(min_spat)
        # new group!
        result.append({'spats': [min_spat], 'fnames': [its_fname]})
        # see if any of the others in the potential group are in:
        for spat, fname in potential_group:
            if spat - min_spat < THRESHOLD: # might want to abs this and double check the sorting
                result[-1]['spats'].append(spat)
                result[-1]['fnames'].append(fname)
                fname_to_spats[fname].remove(spat)

        # filter dict to remove fnames with no spats left
        fname_to_spats = {fname: spats for fname, spats in fname_to_spats.items() if spats}

    return result
