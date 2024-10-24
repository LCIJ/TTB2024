#########################################################################################################################################
## ######################################################################################################################################
## CODE: ORBIT CALCULATOR DELOREAN (v1)
## AUTHOR: MATIAS A. BLANA D.
## CHILE, SANTIAGO JULY 2023
## VERSION : MW DWARF SATELLITES INCLUDING MW MULTIPLE POTENTIALS and M31, AND Fornax Cluster potentials, WITH ORBITS WITH COSMIC EXPANSION OPTIONS
## SCRIPT  : LOADS UNITS AND CONSTANTS 
## #######################################################################################################################################
##########################################################################################################################################
from DELOREAN_pylib import *

# import numpy as np
# from astropy import units as u
# import astropy.coordinates as coord
# from astropy.coordinates import Galactocentric

## UNITS and CONSTANTS V1.0 ###########################################################
#######################################################################################
###### SPATIAL ########################################################################
kpc2pc      = 1000.
pc2kpc      = 1./kpc2pc
pc2km       = 3.085677E+13
kpc2km      = pc2km*1.0E+3
km2kpc      = 1./kpc2km
km2m        = 1.E3
m2km        = 1./km2m
m2cm        = 1.E2
cm2m        = 1./m2cm
pc2cm       = 3.0856775814913673e+18
cm2pc       = 1./pc2cm
kpc2cm      = kpc2pc*pc2cm
cm2kpc      = 1./kpc2cm
kpc2m       = kpc2cm*cm2m
m2kpc       = 1./kpc2m
au2km       = 1.496e+8
km2au       = 1./au2km
###### SPATIAL ANGULAR ################################################################
# Dsun_A18    = 420.
# Dsun_C1     = 409. #kpc
deg2as      = 60.*60
as2deg      = 1./deg2as
# as2pc_A18   = (np.pi/180.)/60./60./np.arctan(1./(Dsun_A18*1E3))  # A18 took it from Irwin et al 2007 [I07]
# as2pc_C12   = (np.pi/180.)/60./60./np.arctan(1./(Dsun_C1*1E3))   #Heliocentric (Clementini et al 2012 [C12])
# as2pc     = asec2pc_A18
# as2pc       = as2pc_C12
# pc2as       = 1./as2pc
# as2kpc      = as2pc*pc2kpc
# kpc2as      = 1./as2kpc
# SLcorr      = (Dsun_C1/Dsun_A18)**2
# pc2deg      = pc2as*as2deg
# deg2pc      = 1./pc2deg
###### TIME ##########################################################################
yr2s              = 365.25*24.*60.*60
s2yr              = 1./yr2s
Myr2yr            = 1.E6
yr2Myr            = 1./Myr2yr
Myr2s             = Myr2yr*yr2s
s2Myr             = 1./Myr2s
Gyr2Myr           = 1.E3
Myr2Gyr           = 1./Gyr2Myr
Gyr2s             = Gyr2Myr*Myr2yr
s2Gyr             = 1./Gyr2s
###### VELOCITY #####################################################################
kms2kpcGyr        = 1.02269
kpcGyr2kms        = 1./kms2kpcGyr
kms2pcMyr         = kms2kpcGyr
pcMyr2ms          = 1./kms2pcMyr
kms2kpcMyr        = 1.02269E-3
kpcMyr2kms        = 1./kms2kpcMyr
kms2ms            = 1.E3
ms2kms            = 1./kms2ms
ms2cms            = 1.E2
cms2ms            = 1./ms2cms
kms2cms           = kms2ms*ms2cms
cms2kms           = 1./kms2cms
###### MASS Number den ##############################################################
kg2g              = 1.E3
g2kg              = 1./kg2g
Msun2g            = 1.98847e+33 # Msun in grams
g2Msun            = 1./Msun2g
Msun2kg           = Msun2g*g2kg
kg2Msun           = 1./Msun2kg
H2g               = 1.6735575e-24 ## mass (proton+electron?) HYDROGEN ATOM in grams
mHg               = 1.6735575e-24 ## mass (proton+electron?) HYDROGEN ATOM in grams
mPg               = 1.6726219236951E-24 # mass proton in grams
mN0g              = 1.6749274980495E-24 # mass neutron (N0 or n0 or  n) in grams
mHeg              = 6.6464764e-24 ## MASS (less than 4x Hydrogen)
NHtoMsun          = 8.4152384E-58
MsuntoNH          = 1./NHtoMsun
###### SURFACE DENSITY #####################
Msunkpc2toMsunpc2 = 1./(kpc2pc**2)
Msunpc2toMsunkpc2 = 1./Msunkpc2toMsunpc2
NHcm2toMsunpc2    = NHtoMsun/(cm2pc**2)
Msunpc2toNHcm2    = 1./NHcm2toMsunpc2
NHcm2toMsunkpc2   = NHcm2toMsunpc2*Msunpc2toMsunkpc2
Msunkpc2toNHcm2   = 1./NHcm2toMsunkpc2
###### VOLUME DENSITY ######################
Msunkpc3toMsunpc3 = 1./(kpc2pc**3)
Msunpc3toMsunkpc3 = 1.E9 #1./Mkpc3toMpc3
Msunkpc3togcm3    = Msun2g/(kpc2cm**3)
gcm3toMsunkpc3    = 1./Msunkpc3togcm3
Msunpc3togcm3     = Msunpc3toMsunkpc3*Msunkpc3togcm3
gcm3toMsunpc3     = 1./Msunpc3togcm3

NHcm3toMsunpc3    = NHtoMsun/(cm2pc**3)
NHcm3toMsunkpc3   = NHcm3toMsunpc3*Msunpc3toMsunkpc3
Msunkpc3toNHcm3   = 1./NHcm3toMsunkpc3
Msunpc3toNHcm3    = 1./NHcm3toMsunpc3

###### define hydrogen_massfrac       0.76                mass fraction of hydrogen DICE
MHe2MHrat         = 0.25
nHe2nHrat         = MHe2MHrat*mHg/mHeg
MH2Mgas           = 1./(1-MHe2MHrat)
Mgas2MH           = 1./MH2Mgas
muHHe             = 1./mHg * (mHg + nHe2nHrat * mHeg)/(1. + nHe2nHrat) # mean molecular weight of a H and He composite
###### CONSTANTS #########################################################################
G_AU3yr2Msun      = 39.478           # AU^3/yr^2/Msun
G_AUMsunkm2s2     = G_AU3yr2Msun*au2km**2/(yr2s**2)
G_pcMsunkm2s2     = 4.3009125E-3     # pc/M (km/s)^2
G_SI              = 6.6743015e-11    # m m^2/s^2 /kg !!!!
G_cgs             = 6.6743015E-8     # cm3⋅g−1⋅s−2  # G ≈ 6.674×10−8 cm3⋅g−1⋅s−2.
gamma_mongas      = 5./3.            # = 1.6667
kB_cgs            = 1.380649E-16     # erg/K = 1.380649E-16 #g cm^2/s^2/K= 1.38064852 × 10-23 m2 kg s-2 K-1

G          = 4.3E-3 #pc/M (km/s)^2
G2         = 4.3E-3*(0.001) #kpc/M (km/s)^2
G3         = 4.3E-3*(0.001*pc2km) #km/M (km/s)^2
G4         = G2*kms2kpcMyr**2  # kpc/M (kpc/Myr)^2
G5         = 6.67e-11     # kg m^2/s^2 ???
###########################################################################################
# pc2km      = 3.0857E+13
# kpc2km     = pc2km*1.0E+3
# km2kpc     = 1./kpc2km
# kms2kpcGyr = 1.02269
# kpcGyr2kms = 1./kms2kpcGyr
# kms2kpcMyr = 1.02269E-3
# kpcMyr2kms = 1./kms2kpcMyr
# Myr2s      = 365.*24.*60*60.
# Gpc        = 4.3E-3 #pc/M (km/s)^2
# G          = 4.3E-3 #pc/M (km/s)^2
# G2         = 4.3E-3*(0.001) #kpc/M (km/s)^2
# G3         = 4.3E-3*(0.001*pc2km) #km/M (km/s)^2
# G4         = G2 *kms2kpcMyr**2
# G5         = 39.478 # AU/Msun(AU/yr)^2
# # dsun       = 8.3*u.kpc
# # Zsun       = 27./1000.*u.kpc #kpc

# # dsun       = 8.2*u.kpc  # BG16 + Fritz+18
# # Zsun       = 25./1000.*u.kpc #kpc # BG16 + Fritz+18

# NHcm2toMpc2=8.4152384E-58/(3.2408E-19*3.2408E-19)
# NHcm3toMpc3=8.4152384E-58/(3.2408E-19*3.2408E-19*3.2408E-19)
# Mpc3toNHcm3=1./NHcm3toMpc3
###########################################################################################
## Constants for N-body potential and others
# G          = 4.3E-3 #pc/M (km/s)^2
GSI          = 6.67e-11     # kg m m^2/s^2 !!
umdice2Msun= 1.E10
um2Msun    = 1.E10 # new value # 5.5E10 # Msun : Md mass Miyamoto-Nagai MW disk, value from BO2015 (stellar+gas disk) used for 
Msun2um    = 1./um2Msun
ud2kpc     = 1.0 # new value #2.5  # kpc : Rd disk scale length Miyamoto-Nagai MW disk
kpc2ud     = 1.0/ud2kpc
ud2pc      = ud2kpc*kpc2pc
pc2ud      = 1./ud2pc
uv2kms     = np.sqrt(um2Msun*Msun2kg*kpc2ud*GSI/(kpc2m*1e6))
kms2uv     = 1./uv2kms
# ud2as      = ud2kpc*kpc2as
uac2kpcMyrsq= (uv2kms*kms2kpcMyr)**2/ud2kpc
# print "uv2kms=",uv2kms

##########################################################################################################
### DEFINING OBSERVATIONAL PARAMETERS FOR COORDINATE TRANSFORMATIONS 
##########################################################################################################
U_BG16,V_BG16,W_BG16 = 10.,11.,7.0 #km/s: err+-: 1,2,0.5 km/s from Bland-Hawthorn & Gerhard 2016 (BG16) LSR
u_BG16,v_BG16,w_BG16 = 11.,248.,7.3 #km/s: err+-: 1,3,0.5 km/s from Bland-Hawthorn & Gerhard 2016 (BG16) from Fritz+18
# v_sun_BG16  = coord.CartesianRepresentation([u_BG16,v_BG16,w_BG16]*u.km/u.s)
v_sun_BG16  = coord.CartesianDifferential([u_BG16,v_BG16,w_BG16]*u.km/u.s)

R_BG16      = 8.2*u.kpc #kpc: err+-0.1 from Bland-Hawthorn & Gerhard 2016 (BG16)
galcen_distance_BG16 = R_BG16
z_sun_BG16  = 25*u.pc

dsun       = 8.2*u.kpc  # BG16 + Fritz+18
Zsun       = 25./1000.*u.kpc #kpc # BG16 + Fritz+18

Xsun       = (dsun**2-Zsun**2)**0.5

R_A19    = 8.178 # kpc error +-13pc sta +- 22pc sys  GRAVITY coll. Abuter+19
galcen_distance_R_A19 = R_A19*u.kpc
mu_RB04  = 6.379 #err +-0.024mas/yr Read & Brunthaler 2004 RB04
v_RA     = (np.tan(mu_RB04/1000/60/60/180*np.pi)*R_A19*3.086e+16/(365*24*60*60)) #km/s: err+-: combining RB04 and Abut+19
v_sun_RA = coord.CartesianDifferential([u_BG16,v_RA,w_BG16]*u.km/u.s)

coord.Galactocentric_def  = coord.Galactocentric()
coord.Galactocentric_BG16 = coord.Galactocentric(galcen_distance=galcen_distance_BG16,galcen_v_sun=v_sun_BG16,z_sun=z_sun_BG16)
coord.Galactocentric_RA   = coord.Galactocentric(galcen_distance=galcen_distance_R_A19,galcen_v_sun=v_sun_RA,z_sun=z_sun_BG16)

v_sun                   = coord.Galactocentric_def.galcen_v_sun.to_cartesian()
v_sun_dif_BG16          = coord.Galactocentric_BG16.galcen_v_sun.to_cartesian()
v_sun_dif_RA            = coord.Galactocentric_RA.galcen_v_sun.to_cartesian()

# print(v_sun.dot)
# print(v_sun_dif_BG16.dot)
# print(v_sun_dif_RA)