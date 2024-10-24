#########################################################################################################################################
## ######################################################################################################################################
## CODE: ORBIT CALCULATOR DELOREAN (v1)
## AUTHOR: MATIAS A. BLANA D.
## CHILE, SANTIAGO JULY 2023
## VERSION : MW DWARF SATELLITES INCLUDING MW MULTIPLE POTENTIALS and M31, AND Fornax Cluster potentials, WITH ORBITS WITH COSMIC EXPANSION OPTIONS
## SCRIPT  : LOADS PYTHON LIBRARIES 
## #######################################################################################################################################
##########################################################################################################################################
import sys
import h5py
import copy
import pickle 
import os as os
import numpy as np
import time as timer
from matplotlib import cm
from matplotlib import image
import multiprocessing as mp
from astropy import units as u
import matplotlib.pyplot as plt
import astropy.coordinates as coord
from scipy.interpolate import interp1d
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from astropy.coordinates import SkyCoord
from matplotlib.collections import LineCollection
from astropy.cosmology import Planck15, z_at_value

# Same and other librearies, just in case
# import scipy as sc
# from numpy import linalg as LA
# from matplotlib import rc
# # rc('text', usetex=True)
# from decimal import Decimal
# from matplotlib import rc
# from scipy.optimize import curve_fit
# from astropy.coordinates import Galactocentric
# from IPython.core.display import display, HTML
# display(HTML("<style>.container { width:100% !important; }</style>"))
# import commah
# import pickle 

# import shutil
# import import_ipynb