import numpy as np

from dbsp_drp import splicing

def interp_w_error_tester(xmin, xmax, xsize, xpmin, xpmax, xpsize):
    x = np.linspace(xmin, xmax, xsize)
    xp = np.linspace(xpmin, xpmax, xpsize)
    yp = 1 - xp ** 2 + 1.1 * xp + np.sin(xp)
    err_yp = np.random.rand(xpsize)

    y, yerr = splicing.interp_w_error(x, xp, yp, err_yp, False)
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
