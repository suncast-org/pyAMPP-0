#!/usr/bin/env python3

import numpy as np
import astropy.time
import astropy.units as u

from PyQt5.QtWidgets import QApplication
from pyampp.gxbox.magfield_viewer import MagFieldViewer
import h5py

class Box:
    def __init__(self, file_name, res):
        h5f = h5py.File(file_name, "r")

        res = res.to(u.Mm)
        bshape = np.array(h5f["nlfff"]["bx"][:].shape)
        print(bshape)

        self.b3d = {
                "nlfff": {
                    "bx": h5f["nlfff"]["bx"][:].transpose(2, 1, 0),
                    "by": h5f["nlfff"]["by"][:].transpose(2, 1, 0),
                    "bz": h5f["nlfff"]["bz"][:].transpose(2, 1, 0),
                    }
                }
        self.grid_coords = {
                "x": np.linspace(0, bshape[0], bshape[0])*res,
                "y": np.linspace(0, bshape[1], bshape[1])*res,
                "z": np.linspace(0, bshape[2], bshape[2])*res,
                }
        print(self.grid_coords)
        h5f.close()

    def open_window(self):
        app = QApplication([])
        gxbox = MagFieldViewer(self)
        gxbox.show()
        app.exec_()