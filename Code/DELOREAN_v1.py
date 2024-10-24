#########################################################################################################################################
## ######################################################################################################################################
## CODE: ORBIT CALCULATOR DELOREAN (v2)
## AUTHOR: MATIAS A. BLAÑA D. (main developer)
## CO-AUTHORS: ROY BUSTOS (Malin1 Setup developer)
## CHILE, SANTIAGO NOVEMBER 2023
## VERSION : 1.2 MW DWARF SATELLITES INCLUDING MW MULTIPLE POTENTIALS and M31, AND Fornax Cluster potentials, WITH ORBITS WITH COSMIC EXPANSION OPTIONS
## UPDATES: 1.2 INCLUDES FITTING OF OBSERVABLES WITH EMCEE
## SCRIPT  : MAIN SCRIPT THAT LOADS INTEGRATIONS ROUTINES, AND WHERET TO DEFINE OBSERVABLES 
## REFERENCES: PLEASE CITE THE ARTICLE WHERE DELOREAN WAS ORIGINALY PUBLISHED: Blaña et al. 2020, MNRAS, 497,3601-3622
## URL: https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.3601B/abstract
## #######################################################################################################################################
## #######################################################################################################################################
from DELOREAN_pylib import *     ## IMPORT STANDARD PYTHON LIBRARIES
from DELOREAN_UC_v1 import *     ## IMPORT CONSTANTS AND UNITS 
from DELOREAN_IC_v1 import *     ## IMPORT INITIAL CONDITIONS GENERATOR FOR DELOREAN INTEGRATOR
from DELOREAN_FUN_v1 import *    ## IMPORT FUNCTIONS 
from DELOREAN_SETUP_v1 import *  ## IMPORT SETUPS FOR EACH ORBITAL INTEGRATION
from DELOREAN_ORBINT_v1 import * ## IMPORT INTEGRATION SUBROUTINES
printdelorean()