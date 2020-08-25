"""
Automated Reduction Pipeline for one arm of P200 DBSP.
"""

import os
import glob
import shutil

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.visualization import ZScaleInterval, ImageNormalize


from yattag import Doc, indent

from pypeit.pypeitsetup import PypeItSetup
from pypeit import defs
from pypeit import pypeit
from pypeit import msgs
from pypeit import sensfunc
from pypeit.spec2dobj import Spec2DObj
from pypeit.specobjs import SpecObjs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts.flux_calib import read_fluxfile
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

def splice(args):
    """
    Splices red and blue spectra together.
    """
    for target, (bluefile, redfile) in args['splicing_dict'].items():
        with open(bluefile) as spec_b:
            with open(redfile) as spec_r:
                final_wvs, final_flam, final_flam_sig = adjust_and_combine_overlap(spec_b, spec_r)
                # now write this to disk
                t = Table([final_wvs, final_flam, final_flam_sig], names=('OPT_WAVE_SPLICED', 'OPT_FLAM_SPLICED', 'OPT_FLAM_SPLICED_SIG'))
                # TODO: make output fits file nicer, possibly copy header from one of the blue/red files
                t.write(f'{target}.fits', format='fits')


def adjust_and_combine_overlap(spec_b, spec_r):
    """

    """
    # combination steps
    overlap_lo = spec_r[1].data['OPT_WAVE'][0]
    overlap_hi = spec_b[1].data['OPT_WAVE'][-1]
    # maybe need to fix the overlaps?
    # need more finely spaced grid to be completely contained within coarser grid

    lw = 0.5
    plt.figure(figsize=(20,10))

    mask_b = ~((spec_b[1].data['OPT_FLAM_SIG'] > 3) & (spec_b[1].data['OPT_WAVE'] > overlap_lo))
    mask_r = ~((spec_r[1].data['OPT_FLAM_SIG'] > 3) & (spec_r[1].data['OPT_WAVE'] < overlap_hi))

    olap_r = (spec_r[1].data['OPT_WAVE'] < overlap_hi)
    olap_b = (spec_b[1].data['OPT_WAVE'] > overlap_lo)

    red_mult=1
    red_mult = sigma_clipped_stats(spec_b[1].data['OPT_FLAM'][olap_b])[1]/sigma_clipped_stats(spec_r[1].data['OPT_FLAM'][olap_r])[1]
    #red_mult = np.average(spec_aag[1].data['OPT_FLAM'][olap_b], weights=spec_aag[1].data['OPT_FLAM_SIG'][olap_b] ** -2.0)/np.average(spec_aag_red[1].data['OPT_FLAM'][olap_r], weights=spec_aag_red[1].data['OPT_FLAM_SIG'][olap_r] ** -2.0)

    # different dispersion.
    wvs_b = spec_b[1].data['OPT_WAVE'][~olap_b]
    wvs_r = spec_r[1].data['OPT_WAVE'][~olap_r]
    flam_b = spec_b[1].data['OPT_FLAM'][~olap_b]
    flam_r = spec_r[1].data['OPT_FLAM'][~olap_r]
    flam_sig_b = spec_b[1].data['OPT_FLAM_SIG'][~olap_b]
    flam_sig_r = spec_r[1].data['OPT_FLAM_SIG'][~olap_r]


    olap_wvs_r = spec_r[1].data['OPT_WAVE'][olap_r]
    olap_flam_r = spec_r[1].data['OPT_FLAM'][olap_r]
    olap_flam_sig_r = spec_r[1].data['OPT_FLAM_SIG'][olap_r]
    olap_wvs_b = spec_b[1].data['OPT_WAVE'][olap_b][:-1]
    olap_flam_b = spec_b[1].data['OPT_FLAM'][olap_b][:-1]
    olap_flam_sig_b = spec_b[1].data['OPT_FLAM_SIG'][olap_b][:-1]

    olap_flam_r_interp = np.interp(olap_wvs_b, olap_wvs_r, olap_flam_r)
    olap_flam_r_interp, olap_flam_sig_r_interp = interp_w_error(olap_wvs_b, olap_wvs_r, olap_flam_r, olap_flam_sig_r)

    olap_flams = np.array([olap_flam_r_interp, olap_flam_b])
    sigs = np.array([olap_flam_sig_r_interp, olap_flam_sig_b])
    weights = sigs ** -2.0

    olap_flam_avgd = np.average(olap_flams, axis=0, weights=weights)
    olap_flam_sig_avgd = 1.0 / np.sqrt(np.mean(weights, axis=0))

    final_wvs = np.concatenate((wvs_b, olap_wvs_b, wvs_r))
    final_flam = np.concatenate((flam_b, olap_flam_avgd, flam_r))
    final_flam_sig = np.concatenate((flam_sig_b, olap_flam_sig_avgd, flam_sig_r))

    plt.errorbar(spec_b[1].data['OPT_WAVE'][mask_b], spec_b[1].data['OPT_FLAM'][mask_b], yerr=spec_b[1].data['OPT_FLAM_SIG'][mask_b], lw=lw)
    plt.errorbar(spec_r[1].data['OPT_WAVE'][mask_r], red_mult*spec_r[1].data['OPT_FLAM'][mask_r], yerr=red_mult*spec_r[1].data['OPT_FLAM_SIG'][mask_r], lw=lw)
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Flux (erg/s/cm${}^2$/$\mathring{A}$)')
    plt.title('Fluxed spectrum of ZTF20aagzdjp')
    #plt.ylim(-100,600)
    plt.ylim(0,20)
    #plt.xlim(5000,6000)
    #plt.savefig('ZTF20aagzdjp_fluxed_spectrum.png')


    mask = final_flam_sig < 3
    plt.figure(figsize=(20,10))
    plt.errorbar(final_wvs[mask], final_flam[mask], yerr=final_flam_sig[mask])
    plt.grid()
    plt.ylim(0, 20)
    plt.xlim(3000, 10500)
    #plt.xlim(7000,8000)
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Flux (erg/s/cm${}^2$/$\mathring{A}$)')
    plt.title('Fluxed spectrum of ZTF20aagzdjp')
    #plt.savefig('ZTF20aagzdjp_fluxed_spectrum.png')
    
    return final_wvs, final_flam, final_flam_sig

def interp_w_error(x, xp, yp, err_yp):
    y = np.zeros_like(x)
    yerr = np.zeros_like(x)
    slopes = np.zeros(xp.shape[0] - 1)
    for i in range(len(slopes)):
        slopes[i] = (yp[i+1] - yp[i])/(xp[i+1] - xp[i])
    
    for i in range(len(x)):
        # find the index j into xp such that xp[j-1] < x[i] < xp[j]
        j = np.searchsorted(xp, x[i], side='right')
        if (x[i] == xp[j-1]):
            y[i] = yp[j-1]
            yerr[i] = err_yp[j-1]
        else:
        # now y[i] = xp[j] + slopes[j]*(x[i] - xp[j])
            y[i] = yp[j-1] + slopes[j-1]*(x[i] - xp[j-1])
            yerr[i] = np.sqrt((((x[i] - xp[j])*err_yp[j-1]) ** 2 + ((x[i] - xp[j-1])*err_yp[j]) ** 2) / ((xp[j-1] - xp[j]) ** 2))
    return y, yerr

def save_one2dspec(ax, spec, edges, traces):
    all_left, all_right = edges
    
    norm = ImageNormalize(spec, interval=ZScaleInterval())
    im = ax.imshow(spec, origin='upper', norm=norm, cmap='gray')
    #im, norm = imshow_norm(spec, interval=ZScaleInterval(), cmap='gray')

    plt.axis('off')

    xs = np.arange(spec.shape[0])

    for i in range(all_left.shape[1]):
        ax.plot(all_left[:, i], xs, 'green', lw=1)
        ax.plot(all_right[:, i], xs, 'red', lw=1)

    if traces is not None:
        for trace in traces:
            ax.plot(trace, xs, 'orange', lw=1)

def save_2dspecs(args):
    arm = 'blue' if 'blue' in args['spectrograph'] else 'red'
    paths = glob.glob(os.path.join(args['output_path'], 'Science', f'spec2d_{arm}*.fits'))
    
    out_path = os.path.join(args['output_path'], 'QA', 'PNGs', 'Extraction')
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    for path in paths:
        # open fits file
        spec = Spec2DObj.from_file(path, 1)
        spec1d_file = path.replace('spec2d', 'spec1d')
        if os.path.isfile(spec1d_file):
            spec1ds = SpecObjs.from_fitsfile(spec1d_file)
        else:
            spec1ds = None
            msgs.warn('Could not find spec1d file: {:s}'.format(spec1d_file) + msgs.newline() +
                  '             No objects were extracted.')
        
        base_name = os.path.join(out_path, os.path.basename(path).split('_')[1])

        all_left, all_right, _ = spec.slits.select_edges()
        edges = [all_left, all_right]
        traces = [spec1ds[i]['TRACE_SPAT'] for i in range(len(spec1ds))] if spec1ds is not None else None
        mask = spec.bpmmask == 0

        out_shape = spec.sciimg.shape
        fig_width = 0.15 + 4*out_shape[1]/100.
        subfig_width = out_shape[1]/100. / fig_width
        padding_width = 0.15 / fig_width
        fig = plt.figure(figsize=(fig_width, out_shape[0]/100.), dpi=100)

        # science image
        #ax = fig.add_subplot(1, 4, 1)
        ax = fig.add_axes([0, 0, subfig_width, 1])
        save_one2dspec(ax, spec.sciimg, edges, traces)

        # sky subtracted science image
        #ax = fig.add_subplot(1, 4, 2)
        ax = fig.add_axes([subfig_width + padding_width, 0, subfig_width, 1])
        save_one2dspec(ax, (spec.sciimg - spec.skymodel) * mask, edges, traces)

        # sky subtracted images divided by noise
        #ax = fig.add_subplot(1, 4, 3)
        ax = fig.add_axes([2*(subfig_width + padding_width), 0, subfig_width, 1])
        save_one2dspec(ax, (spec.sciimg - spec.skymodel) * np.sqrt(spec.ivarmodel) * mask, edges, traces)

        # total residauls
        #ax = fig.add_subplot(1, 4, 4)
        ax = fig.add_axes([3*(subfig_width + padding_width), 0, subfig_width, 1])
        save_one2dspec(ax, (spec.sciimg - spec.skymodel - spec.objmodel) * np.sqrt(spec.ivarmodel) * mask, edges, traces)

        # save figure
        fig.savefig(f'{base_name}_extraction.png', bbox_inches='tight', pad_inches=0)
        plt.close(fig)

def write_extraction_QA(args):
    out_path = os.path.join(args['output_path'], 'QA')

    doc, tag, text = Doc().tagtext()

    doc.asis('<!DOCTYPE html>')
    
    pngs = [os.path.basename(fullpath) for fullpath in glob.glob(os.path.join(out_path, 'PNGs', 'Extraction', '*.png'))]
    pngs = sorted(pngs)
    relpath = [os.path.join('PNGs', 'Extraction', png) for png in pngs]

    expnames = [png.split('-')[0] for png in pngs]
    objnames = [png.split('-')[1].split('_')[0] for png in pngs]


    with tag('html'):
        with tag('head'):
            with tag('title'):
                text("DBSP Extraction QA")
            doc.stag('link', rel='stylesheet', href='dbsp_qa.css')
            with tag('script', src='./dbsp_qa.js'):
                pass
        with tag('body'):
            with tag('div'):
                doc.attr(klass='tab')
                for expname in expnames:
                    with tag('button'):
                        doc.attr(klass='tablinks', onclick=f"openExposure(event, '{expname}')")
                        text(expname)
            for i in range(len(pngs)):
                with tag('div'):
                    doc.attr(klass='tabcontent', id=expnames[i])
                    with tag('h1'):
                        text(f'{expnames[i]} {objnames[i]}')
                    doc.stag('img', src=relpath[i])
    
    result = indent(doc.getvalue())
    with open(os.path.join(out_path, 'Extraction.html'), mode='wt') as f:
        f.write(result)
    
    shutil.copy("dbsp_qa.js", os.path.join(out_path, "dbsp_qa.js"))
    shutil.copy("dbsp_qa.css", os.path.join(out_path, "dbsp_qa.css"))
