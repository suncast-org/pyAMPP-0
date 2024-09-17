#!/usr/bin/env python3

import numpy as np
import astropy.units as u
from astropy.time import Time
from pathlib import Path
import h5py
import os
from gximagecomputing import GXRadioImageComputing

class RadioImageComputing(GXRadioImageComputing):
    def __init__(self, libname=None):
        super().__init__(libname)

    def load_model_hdf(self, file_name):
        model_f = h5py.File(file_name, "r")
        chromo_box = model_f["chromo"]
        header = dict(chromo_box.attrs)

        chromo_box_loaded = {}
        for k in chromo_box.keys():
            try:
                chromo_box_loaded[k] = chromo_box[k][:]
            except ValueError:
                chromo_box_loaded[k] = chromo_box[k][()]
        model_f.close()
        header["obs_time"] = Time(header["obs_time"])

        return self.load_model_dict(chromo_box_loaded, header)