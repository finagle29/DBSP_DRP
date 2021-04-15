import os
import warnings

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

from dbsp_drp import splicing

def interp_w_error_tester(xmin, xmax, xsize, xpmin, xpmax, xpsize):
    x = np.linspace(xmin, xmax, xsize)
    xp = np.linspace(xpmin, xpmax, xpsize)
    yp = 1 - xp ** 2 + 1.1 * xp + np.sin(xp)
    err_yp = np.random.rand(xpsize)

    y, yerr = splicing.interp_w_error(x, xp, yp, err_yp)
    y_numpy = np.interp(x, xp, yp)
    assert np.all(y == y_numpy)

def test_interp_w_error():
    interp_w_error_tester(0, 10, 100, -1, 11, 23)
    interp_w_error_tester(0, 10, 100, -1, 11, 1023)
    interp_w_error_tester(0, 10, 100, 2, 9, 23)
    interp_w_error_tester(0, 10, 100, 2, 9, 1023)

    interp_w_error_tester(0, 10, 100, -1, 11, 1)
    interp_w_error_tester(0, 10, 1, -1, 11, 100)
    interp_w_error_tester(0, 10, 100, 2, 9, 1)
    interp_w_error_tester(0, 10, 1, 2, 9, 100)

def test_aco_all(tmp_path):
    def write_test_data(wv_min, wv_max, mul, name):
        wv = np.linspace(wv_min, wv_max, int(1000*mul)) + np.random.random(int(1000*mul))
        flux = (0.001 * wv + 1.0) * mul + np.random.random(int(1000*mul))
        ivar = np.ones_like(wv)

        col_wvs = fits.Column(name='wave', array=wv, format='D')
        col_flux = fits.Column(name='flux', array=flux, format='D')
        col_error = fits.Column(name='ivar', array=ivar, format='D')
        table_hdu = fits.BinTableHDU.from_columns([col_wvs, col_flux, col_error])
        hdulist = fits.HDUList([fits.PrimaryHDU(), table_hdu])
        hdulist.writeto(os.path.join(tmp_path, f"{name}.fits"))

    write_test_data(2000, 4500, 1, 'u')
    write_test_data(4000, 6500, 1.1, 'g')
    write_test_data(6000, 8500, 1.2, 'r')
    write_test_data(8000, 11500, 1.3, 'i')

    paths = [os.path.join(tmp_path, f"{name}.fits") for name in ['u', 'g', 'r', 'i']]

    ((wv, flux, sig), _) = splicing.adjust_and_combine_overlap_all(paths, target="test")

    paths[1] = None
    ((wv, flux, sig), _) = splicing.adjust_and_combine_overlap_all(paths, target="test")

    paths[2] = None
    ((wv, flux, sig), _) = splicing.adjust_and_combine_overlap_all(paths, target="test")

    paths[0] = None
    ((wv, flux, sig), _) = splicing.adjust_and_combine_overlap_all(paths, target="test")
