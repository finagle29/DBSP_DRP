# Interactive Table Editor
import os
from datetime import datetime
from typing import Union, Tuple, List

from astropy.table import Table
import matplotlib.pyplot as plt
from numpy import ma

from PySide2 import QtCore, QtGui, QtWidgets
from PySide2.QtCore import Qt

from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz, ICRS
from astropy.time import Time
import astropy
from astropy import units as u
from astropy.io import fits


loc = EarthLocation.of_site('Palomar')

def update_airmass(row: astropy.table.Row) -> None:
    """
    Updates the ``airmass`` entry of ``row`` based on the ``ra`` ``dec`` and
    ``mjd`` entries.

    Args:
        row (astropy.table.Row): Row to be updated
    """
    skycoord = SkyCoord(row['ra'], row['dec'], unit=u.deg)
    time = Time(row['mjd'], format='mjd')
    altaz = AltAz(obstime=time, location=loc)
    row['airmass'] = skycoord.transform_to(altaz).secz

def get_zenith_ra_dec(time: str) -> astropy.coordinates.SkyCoord:
    """
    Returns RA and Dec of the zenith at Palomar at the given time

    Args:
        time (str): MJD time

    Returns:
        astropy.coordinates.SkyCoord: object containing the zenith RA and Dec
    """
    time = Time(time, format='mjd')
    altaz = AltAz(alt=Angle(90, unit=u.deg), az=Angle(0, unit=u.deg), obstime=time, location=loc)
    return altaz.transform_to(ICRS)

class TableModel(QtCore.QAbstractTableModel):
    def __init__(self, data: Table, cols: tuple = None, del_files: list = None):
        """
        Create TableModel object.

        Args:
            data (Table): data table represented by this object.
            cols (tuple, optional): Columns to display. Defaults to None, in
                which case all columns are displayed.
            del_files (list, optional): List of deleted files/rows. Defaults to None.
        """
        super(TableModel, self).__init__()
        # would use
        # self._data = Table(data, masked = True, copy = False)
        # but it actually copies for some reason? --MSR
        self._data = data
        self._cols = cols if cols is not None else data.colnames

        if data.masked:
            self._mask = self._data.mask
        else:
            self._mask = Table()


        for col in self._cols:
            self._mask.add_column(self._data[col] == None, name=col)

        self._col_count = len(self._cols)

        self._modified_files = set()
        self._deleled_files = del_files if del_files is not None else []

    def data(self, index: QtCore.QModelIndex, role) -> Union[str, QtGui.QColor, None]:
        if role in (Qt.DisplayRole, Qt.EditRole):
            col = self._cols[index.column()]
            row = index.row()
            value = self._data[col][row]
            mask = self._mask[col][row]
            if mask:
                return str(ma.masked)
            if col == 'ra':
                ra = Angle(value, unit=u.deg)
                value = ra.to_string(unit=u.hour)
            elif col == 'dec':
                dec = Angle(value, unit=u.deg)
                value = dec.to_string()
            elif col == 'dispangle':
                ang = Angle(value, unit=u.deg)
                value = ang.to_string()
            elif col == 'mjd' or col == 'airmass':
                value = f'{value:.4f}'
            return str(value)
        if role == Qt.BackgroundRole:
            col = self._cols[index.column()]
            masked = self._mask[col][index.row()] or (self._data[col][index.row()] == 'None')
            if masked:
                return QtGui.QColor(Qt.red)

    def setData(self, index: QtCore.QModelIndex, value: object,
                role: Qt.ItemDataRole=Qt.EditRole)-> bool:
        if role != Qt.EditRole:
            return False
        try:
            col = self._cols[index.column()]
            row = index.row()
            am = False
            if col == 'ra':
                ra = Angle(value, unit=u.hour)
                value = ra.degree
                if self._data['dec'][row]:
                    am = True
            elif col == 'dec':
                dec = Angle(value, unit=u.deg)
                value = dec.degree
                if self._data['ra'][row]:
                    am = True
            elif col =='dispangle':
                value = Angle(value, unit=u.deg).degree
            elif col == 'mjd' or col == 'airmass':
                value = float(value)
            self._data[col][row] = value
            self._mask[col][row] = False
            if am:
                update_airmass(self._data[row])
                self._mask['airmass'][row] = False
                # TODO: emit a dataChanged event for this row
            self.dataChanged.emit(index, index)
            self._modified_row(index)
            return True
        except ValueError:
            print(f"Error: could not parse {value} as the required type for column {col}.")
            return False

    def flags(self, index: QtCore.QModelIndex) -> Qt.ItemFlag:
        if not index.isValid():
            return Qt.ItemIsEnabled
        return QtCore.QAbstractTableModel.flags(self, index) | Qt.ItemIsEditable

    def rowCount(self, index: QtCore.QModelIndex) -> int:
        return len(self._data)

    def columnCount(self, index: QtCore.QModelIndex) -> int:
        return self._col_count

    def headerData(self, section, orientation: Qt.Alignment, role: Qt.ItemDataRole) -> str:
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self._cols[section]
            if orientation == Qt.Vertical:
                return ''

    def removeRows(self, position: int, rows: int, index: QtCore.QModelIndex = QtCore.QModelIndex()) -> bool:
        self.beginRemoveRows(index.parent(), position, position + rows - 1)

        for row in range(rows):
            self._deleled_files.append(self._data['filename'][position])
            self._data.remove_row(position)
            self._mask.remove_row(position)

        self.endRemoveRows()
        return True

    def frametype(self, index: QtCore.QModelIndex) -> str:
        return self._data['frametype'][index.row()]

    def _set_data_to_zenith(self, index: QtCore.QModelIndex) -> None:
        row = index.row()
        time = self._data['mjd'][row]
        ra_dec = get_zenith_ra_dec(time)
        self._data['ra'][row] = ra_dec.ra.deg
        self._data['dec'][row] = ra_dec.dec.deg
        self._data['airmass'][row] = 1.0

        self._modified_files.add(self._data['filename'][row])

    def _modified_row(self, index: QtCore.QModelIndex) -> None:
        row = index.row()
        self._modified_files.add(self._data['filename'][row])

    def _update_fits(self) -> None:
        def update_header(header: fits.Header, keyword, new, comment):
            old = header.get(keyword)
            print(f'old {keyword}: {old}\tnew {keyword}: {new}')
            if old != new:
                header[keyword] = new
                header.comments[keyword] = comment

        dt = datetime.now()
        now_str = dt.isoformat(timespec='seconds')
        for fname in self._modified_files:
            row = self._data[self._data['filename'] == fname]
            maskrow = self._mask[self._data['filename'] == fname]
            path = os.path.join(row['directory'][0], fname)
            with fits.open(path, mode='update') as hdul:
                if not maskrow['ra']:
                    update_header(hdul[0].header, 'RA',
                        Angle(row['ra'][0], unit=u.deg).to_string(unit=u.hour, sep=':'),
                        f'Right Ascension. Modified {now_str}')
                if not maskrow['dec']:
                    update_header(hdul[0].header, 'DEC',
                        Angle(row['dec'][0], unit=u.deg).to_string(unit=u.deg, sep=':',
                            alwayssign=True),
                        f'Declination. Modified {now_str}')
                if not maskrow['airmass']:
                    update_header(hdul[0].header, 'AIRMASS', f"{row['airmass'][0]:.3f}",
                        f'Airmass. Modified {now_str}')
                if not maskrow['target']:
                    update_header(hdul[0].header, 'OBJECT', row['target'][0],
                        f'object title. Modified {now_str}')
                if not maskrow['frametype'] and row['frametype'][0] in ['science', 'standard']:
                    update_header(hdul[0].header, 'IMGTYPE', 'object',
                        f'frame type. Modified {now_str}')
                if not maskrow['dispangle']:
                    update_header(hdul[0].header, 'ANGLE',
                        Angle(row['dispangle'][0], unit=u.deg).to_string(unit=u.deg, sep=(' deg ', ' min'), fields=2),
                        f'Grating Angle. Modified {now_str}')
                if not maskrow['dispname']:
                    update_header(hdul[0].header, 'GRATING', row['dispname'][0], f'lines/mm & Blaze. Modified {now_str}')


class TableView(QtWidgets.QTableView):
    def __init__(self, model: QtCore.QAbstractTableModel, parent=None):
        super().__init__(parent)
        self.setModel(model)

        self.delegate = Delegate(self, model._cols)
        self.setItemDelegate(self.delegate)

        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.setEditTriggers(QtWidgets.QAbstractItemView.DoubleClicked)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.show_context_menu)
        self.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding)
        self.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)

    def show_context_menu(self, point: QtCore.QPoint) -> None:
        index = self.indexAt(point)
        if index.isValid():
            menu = QtWidgets.QMenu(self)
            menu.setAttribute(Qt.WA_DeleteOnClose)

            frametype = self.model().frametype(index)
            if any([x in frametype for x in ['bias', 'arc', 'flat']]):
                action = menu.addAction('Set RA/Dec and Airmass to Zenith')
                action.triggered.connect(lambda: self.model()._set_data_to_zenith(index))

            selected_rows = sorted(set(ix.row() for ix in self.selectionModel().selectedIndexes()), reverse=True)
            action = menu.addAction(f'Delete row{"s" if len(selected_rows) > 1 else ""}')
            action.triggered.connect(lambda: [self.model().removeRow(row) for row in selected_rows] +
                [self.selectionModel().clear()])

            menu.popup(self.viewport().mapToGlobal(point))


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, data: Table, cols: tuple = None, del_files: list = None):
        super().__init__()

        self.table = TableView(TableModel(data, cols, del_files))

        self.setCentralWidget(self.table)


        self.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding)

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self.table.model()._update_fits()


class Delegate(QtWidgets.QStyledItemDelegate):
    def __init__(self, parent, cols: tuple):
        super().__init__(parent)
        self._cols = cols

    def createEditor(self, parent, option, index: QtCore.QModelIndex):
        editor = QtWidgets.QLineEdit(parent)
        editor.setFrame(False)
        return editor

    def setEditorData(self, editor: QtWidgets.QLineEdit, index: QtCore.QModelIndex) -> None:
        value = index.model().data(index, Qt.EditRole)
        editor.setText(value)

    def setModelData(self, editor: QtWidgets.QLineEdit, model: TableModel, index: QtCore.QModelIndex) -> None:
        value = editor.text()
        model.setData(index, value)

    def updateEditorGeometry(self, editor: QtWidgets.QLineEdit, option: QtWidgets.QStyleOptionViewItem, index: QtCore.QModelIndex) -> None:
        editor.setGeometry(option.rect)

def main(table: Table, del_files: List[str]):
    """
    Opens header/metadata editing table GUI.

    Args:
        table (Table): table containing metadata
        del_files (List[str]): List of files to be deleted. Mutated!
    """
    cols = ('filename', 'frametype', 'ra', 'dec', 'target', 'dispname',
        'binning', 'mjd', 'airmass', 'exptime', 'dispangle', 'dichroic',
        'slitwid', 'calib')
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication([])
    else:
        app = QtWidgets.QApplication.instance()
    window = MainWindow(table, cols, del_files)
    window.show()
    app.exec_()

# TODO:
#   - default size / column sizing
