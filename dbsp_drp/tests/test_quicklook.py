import os

import pytest

from dbsp_drp import quicklook

def quicklook_tester(wdir, path, fname):
    os.chdir(wdir)
    args = quicklook.parse([os.path.join(path, 'ql_test_data', fname), '--no-show'])
    quicklook.main(args)

@pytest.fixture(scope='session')
def working_dir(tmp_path_factory):
    return tmp_path_factory.mktemp('ql_out')

def test_quicklook_red_calib(files_path, working_dir):
    quicklook_tester(working_dir, files_path, 'red')

def test_quicklook_blue_calib(files_path, working_dir):
    quicklook_tester(working_dir, files_path, 'blue')

def test_quicklook_red(files_path, working_dir):
    quicklook_tester(working_dir, files_path, 'red0052.fits')

def test_quicklook_blue(files_path, working_dir):
    quicklook_tester(working_dir, files_path, 'blue0050.fits')
