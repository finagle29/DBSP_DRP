import os
import shutil
from typing import Callable

import pytest
import astropy.table
from astropy.io import fits
import pypeit.pypeitsetup
from PySide2.QtCore import Qt
from PySide2 import QtGui

from dbsp_drp import table_edit

DEFAULT_COLS = ('filename', 'frametype', 'ra', 'dec', 'target', 'dispname', 'binning', 'mjd',
                'airmass', 'exptime', 'dispangle', 'dichroic', 'slitwid')


@pytest.fixture(scope="module")
def pypeit_setup(tmp_path_factory, files_path) -> pypeit.pypeitsetup.PypeItSetup:
    tmp_path = tmp_path_factory.mktemp("table_edit_test_data_copy")
    orig_files = os.path.join(files_path, 'table_edit_test_data')
    tmp_root_dir = shutil.copytree(orig_files, tmp_path, dirs_exist_ok=True)
    ps = pypeit.pypeitsetup.PypeItSetup.from_file_root(tmp_root_dir, 'p200_dbsp_red',
                                        extension='.fits', output_path=tmp_root_dir)
    ps.run(setup_only=True, sort_dir=tmp_root_dir)
    return ps

@pytest.fixture
def astropy_table(pypeit_setup: pypeit.pypeitsetup.PypeItSetup) -> astropy.table.Table:
    table = pypeit_setup.fitstbl.table.copy()
    table.sort('filename')
    return table

@pytest.fixture
def table_model(astropy_table: astropy.table.Table) -> Callable[[tuple, list], table_edit.TableModel]:
    def _table_maker(cols: tuple, del_files: list) -> table_edit.TableModel:
        return table_edit.TableModel(astropy_table, cols, del_files)

    return _table_maker

def test_table_model_data(table_model: table_edit.TableModel) -> None:
    model = table_model(DEFAULT_COLS, [])
    ix = model.createIndex(0, 4)
    assert model.data(ix, Qt.DisplayRole) == "bias"
    assert model.data(ix, Qt.EditRole) == "bias"
    ix = model.createIndex(0, 2)
    assert model.data(ix, Qt.DisplayRole) == "--"
    assert model.data(ix, Qt.BackgroundRole) == QtGui.QColor(Qt.red)
    ix = model.createIndex(0, 10)
    assert model.data(ix, Qt.DisplayRole) == "24d38m00s"
    assert model.data(ix, Qt.EditRole) == "24d38m00s"

def test_table_model_setdata(table_model: table_edit.TableModel) -> None:
    model = table_model(DEFAULT_COLS, [])
    ix = model.createIndex(0, 2)
    assert not model.setData(ix, "abc")
    assert model.setData(ix, "4h20m")
    assert model.data(ix, Qt.DisplayRole) == "4h20m00s"
    ix_am = model.createIndex(0, 8)
    assert model.data(ix_am, Qt.BackgroundRole) == QtGui.QColor(Qt.red)
    ix = model.createIndex(0, 3)
    assert not model.setData(ix, "abc")
    assert model.setData(ix, "4d20m")
    assert model.data(ix, Qt.DisplayRole) == "4d20m00s"

    assert model.data(ix_am, Qt.BackgroundRole) is None

    assert not model.setData(ix_am, "abc")
    assert model.setData(ix_am, "1.234")
    assert float(model.data(ix_am, Qt.DisplayRole)) == float("1.234")

def test_updating_fits_header(table_model: table_edit.TableModel) -> None:
    model = table_model(DEFAULT_COLS, [])

    ix = model.createIndex(0, 4)
    model.setData(ix, "foo")

    ix = model.createIndex(0, 2)
    model.setData(ix, "4h20m")

    ix = model.createIndex(0, 3)
    model.setData(ix, "4d20m")

    ix_fname = model.createIndex(0,0)
    fname = model.data(ix_fname, Qt.DisplayRole)
    # secret knowledge!!
    path = os.path.join(model._data['directory'][0], fname)

    model._update_fits()

    hdul = fits.open(path)
    assert hdul[0].header['OBJECT'] == "foo"
    assert hdul[0].header['RA'] == "4:20:00"
    assert hdul[0].header['DEC'] == "+4:20:00"
    hdul.close()
