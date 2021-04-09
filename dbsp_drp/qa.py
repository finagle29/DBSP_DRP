"""
Quality Assurance for automated reduction of P200 DBSP.
"""

import os
import glob
import shutil
from typing import Tuple, List

from pkg_resources import resource_filename

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.visualization import ZScaleInterval, ImageNormalize

from yattag import Doc, indent

from pypeit import msgs
from pypeit.spec2dobj import Spec2DObj
from pypeit.specobjs import SpecObjs

def save_one2dspec(ax: plt.Axes, spec: np.ndarray, edges: Tuple[np.ndarray, np.ndarray], traces: List[np.ndarray]) -> None:
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

def save_2dspecs(qa_dict: dict, output_spec2ds: List[str], output_path: str) -> None:
    obj_png_dict = qa_dict

    out_path = os.path.join(output_path, 'QA', 'PNGs', 'Extraction')
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    for fname in output_spec2ds:
        path = os.path.join(output_path, 'Science', fname)
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
        fig_width = 0.15 + 5*out_shape[1]/100.
        subfig_width = out_shape[1]/100. / fig_width
        padding_width = 0.15 / fig_width
        fig = plt.figure(figsize=(fig_width, out_shape[0]/100.), dpi=100)

        # science image
        ax = fig.add_axes([0, 0, subfig_width, 1])
        save_one2dspec(ax, spec.sciimg, edges, traces)

        # sky model
        ax = fig.add_axes([subfig_width + padding_width, 0, subfig_width, 1])
        save_one2dspec(ax, spec.skymodel * mask, edges, traces)

        # sky subtracted science image
        ax = fig.add_axes([2*(subfig_width + padding_width), 0, subfig_width, 1])
        save_one2dspec(ax, (spec.sciimg - spec.skymodel) * mask, edges, traces)

        # sky subtracted images divided by noise
        ax = fig.add_axes([3*(subfig_width + padding_width), 0, subfig_width, 1])
        save_one2dspec(ax, (spec.sciimg - spec.skymodel) * np.sqrt(spec.ivarmodel) * mask, edges, traces)

        # total residauls
        ax = fig.add_axes([4*(subfig_width + padding_width), 0, subfig_width, 1])
        save_one2dspec(ax, (spec.sciimg - spec.skymodel - spec.objmodel) * np.sqrt(spec.ivarmodel) * mask, edges, traces)

        # save figure
        outname = f'{base_name}_extraction.png'
        msgs.info(f"Saving {os.path.basename(outname)}")
        fig.savefig(outname, bbox_inches='tight', pad_inches=0)
        plt.close(fig)

        airmass = spec.head0['AIRMASS']
        ut = spec.head0['UTSHUT']
        fname = spec.head0['FILENAME']

        if obj_png_dict.get(spec.head0['OBJECT']):
            obj_png_dict[spec.head0['OBJECT']]['pngs'].append(outname)
            obj_png_dict[spec.head0['OBJECT']]['airmass'].append(airmass)
            obj_png_dict[spec.head0['OBJECT']]['time'].append(ut)
            obj_png_dict[spec.head0['OBJECT']]['fname'].append(fname)
        else:
            obj_png_dict[spec.head0['OBJECT']] = {
                'pngs': [outname],
                'airmass': [airmass],
                'time': [ut],
                'fname': [fname]
            }

    return obj_png_dict

def write_extraction_QA(qa_dict: dict, output_path: str, spectrograph: str) -> None:
    out_path = os.path.join(output_path, 'QA')

    doc, tag, text = Doc().tagtext()
    doc.asis('<!DOCTYPE html>')

    with tag('html'):
        with tag('head'):
            with tag('title'):
                text(f"{spectrograph} Extraction QA")
            doc.stag('link', rel='stylesheet', href='dbsp_qa.css')
            with tag('script', src='./dbsp_qa.js'):
                pass
        with tag('body'):
            with tag('div'):
                doc.attr(klass='tab')
                for target in qa_dict:
                    with tag('button'):
                        doc.attr(klass='tablinks', onclick=f"openExposure(event, '{target}')")
                        text(target)
            for target, obj_dict in qa_dict.items():
                with tag('div'):
                    doc.attr(klass='tabcontent', id=target)
                    with tag('h1'):
                        text(target)
                    for i, png in enumerate(obj_dict['pngs']):
                        with tag('h2'):
                            text(obj_dict['fname'][i])
                            text('\t')
                            text(obj_dict['time'][i])
                            text('\t')
                            text(obj_dict['airmass'][i])
                        doc.stag('img', src=os.path.join('PNGs', 'Extraction', os.path.basename(png)))

    msgs.info("Writing Extraction QA page")
    result = indent(doc.getvalue())
    extraction_page = os.path.join(out_path, 'Extraction.html')
    with open(extraction_page, mode='wt') as f:
        f.write(result)

    shutil.copy(resource_filename("dbsp_drp", "/data/dbsp_qa.js"), os.path.join(out_path, "dbsp_qa.js"))
    shutil.copy(resource_filename("dbsp_drp", "/data/dbsp_qa.css"), os.path.join(out_path, "dbsp_qa.css"))
    msgs.info(f"Extraction QA page available at {extraction_page}")
