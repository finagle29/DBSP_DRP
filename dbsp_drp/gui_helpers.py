"""
Module for classes and functions shared between GUIs.
"""

from matplotlib.backend_tools import ToolBase
from PySide2 import QtWidgets

class HelpTool(ToolBase):
    """
    Print help text.
    """
    # keyboard shortcut
    default_keymap = 'h'
    description = 'Help'
    image = 'help'

    def __init__(self, *args, helptext="", **kwargs):
        self.helptext = helptext
        super().__init__(*args, **kwargs)

    def trigger(self, *args, **kwargs):
        QtWidgets.QMessageBox.information(None, "Help", self.helptext)
