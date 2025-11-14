"""Centralized default inputs and grids for OCRA simulations."""

import numpy as np


PH_DEFAULTS = {
    "DIV": "Ca",
    "totP": 1,
    "Temp": 288,
    "nDIV": 1,
    "nSiO2": 1,
    "totnum": 100,
    "plot_flag": True,
    "table_flag": True,
    "analytical_flag": False,
    "comparison": "PCO2",
}


CCD_DEFAULTS = {
    "beta": 0.3,
    "nSiO2": 1,
    "nDIV": 1,
    "totnum": 10,
    "numQ1": 10,
    "numQ2": 10,
    "plot_flag": True,
    "table_flag": True,
}


PHASE_DEFAULTS = {
    "DIV": "Ca",
    "nSiO2": 1,
    "Temp": 310,
    "beta": 0.3,
    "totnum": 20,
    "plot_flag": True,
    "table_flag": True,
}


GRID_DEFAULTS = {
    "temps": lambda size: np.linspace(273.16, 372.16, num=size),
    "pco2s": lambda size: np.logspace(-8, -0.5, num=size),
    "totps": lambda size: np.logspace(0, np.log10(5000), num=size),
}
