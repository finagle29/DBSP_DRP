"""
Automated Reduction Pipeline for one arm of P200 DBSP
"""
import os

from pypeit.pypeitsetup import PypeItSetup
from pypeit import defs
from pypeit import pypeit
from pypeit import msgs
from pypeit import sensfunc
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts.flux_calib import read_fluxfile
from astropy.io import fits

from pypeit import fluxcalibrate
from pypeit.par import pypeitpar

def setup(args):
    """
    Does PypeIt setup, without writing the .pypeit file
    """
    # Check that the spectrograph is provided if using a file root
    if args['root'] is not None:
        if args['spectrograph'] is None:
            raise ValueError('Must provide spectrograph identifier with file root.')
        # Check that input spectrograph is supported
        instruments_served = defs.pypeit_spectrographs
        if args['spectrograph'] not in instruments_served:
            raise ValueError('Instrument \'{0}\' unknown to PypeIt.\n'.format(args['spectrograph'])
                             + '\tOptions are: {0}\n'.format(', '.join(instruments_served))
                             + '\tSelect an available instrument or consult the documentation '
                             + 'on how to add a new instrument.')

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

def write_setup(args, context):
    """
    Writes the .pypeit file
    """
    ps, output_path = context
    # Use PypeItMetaData to write the complete PypeIt file
    pypeit_file = os.path.join(output_path, '{0}.pypeit'.format(args['spectrograph']))
    config_list = [item.strip() for item in args['cfg_split'].split(',')]

    #ps.user_cfg += par['calibrations']['master_dir']
    ps.user_cfg.append('[calibrations]')
    ps.user_cfg.append('master_dir = Master_' + args['spectrograph'].split('_')[-1])

    return ps.fitstbl.write_pypeit(pypeit_file, cfg_lines=ps.user_cfg,
                                   write_bkg_pairs=args['background'], configs=config_list)

def redux(args):
    """
    Runs the reduction
    """
    splitnm = os.path.splitext(args['pypeit_file'])
    if splitnm[1] != '.pypeit':
        msgs.error("Bad extension for PypeIt reduction file."+msgs.newline()+".pypeit is required")
    logname = splitnm[0] + ".log"

    pypeIt = pypeit.PypeIt(args['pypeit_file'], verbosity=2,
                           reuse_masters=not args['do_not_reuse_masters'],
                           overwrite=True, #args['overwrite'],
                           redux_path=args['output_path'], #args['redux_path'], # rethink these
                           calib_only=args['calib_only'],
                           logname=logname, show=args['show'])
    pypeIt.reduce_all()
    msgs.info('Data reduction complete')

    msgs.info('Generating QA HTML')
    pypeIt.build_qa()

def make_sensfunc(args):
    """
    Makes a sensitivity function
    """
    par = load_spectrograph(args['spectrograph']).default_pypeit_par()
    outfile = os.path.join(args['output_path'], os.path.basename(args['spec1dfile']).replace('spec1d', 'sens'))
    sensobj = sensfunc.SensFunc.get_instance(args['spec1dfile'], outfile, par=par['sensfunc'], debug=args['debug'])
    sensobj.run()
    sensobj.save()
    return outfile

def build_fluxfile(args):
    """
    Writes the fluxfile for fluxing.
    """
    # args['spec1dfiles'] is a dict mapping spec1d files to the sensitivity function file they should use
    cfg_lines = []
    cfg_lines.append('[fluxcalib]')
    if 'red' in args['spectrograph']:
        cfg_lines.append('extinct_correct = False')
    
    # data section
    cfg_lines.append('flux read')
    for spec1d, sensfun in args['spec1dfiles'].items():
        cfg_lines.append(f'  {spec1d} {sensfun}')
    cfg_lines.append('flux end')

    ofile = os.path.join(args['output_path'], f'{args["spectrograph"]}.flux')

    with open(ofile, mode='wt') as f:
        for line in cfg_lines:
            f.write(line)
            f.write("\n")

    return ofile

def flux(args):
    """
    Fluxes spectra.
    """
    # Load the file
    config_lines, spec1dfiles, sensfiles = read_fluxfile(args['flux_file'])
    # Read in spectrograph from spec1dfile header
    header = fits.getheader(spec1dfiles[0])
    spectrograph = load_spectrograph(header['PYP_SPEC'])

    # Parameters
    spectrograph_def_par = spectrograph.default_pypeit_par()
    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(), merge_with=config_lines)
    # Write the par to disk
    par_outfile = os.path.join(args['output_path'], f"{args['flux_file']}.par")
    print(f"Writing the parameters to {par_outfile}")
    par.to_config(par_outfile)

    # Instantiate
    FxCalib = fluxcalibrate.FluxCalibrate.get_instance(spec1dfiles, sensfiles, par=par['fluxcalib'], debug=args['debug'])
    msgs.info('Flux calibration complete')

