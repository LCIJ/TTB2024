#########################################################################################################################################
## ######################################################################################################################################
## CODE: ORBIT CALCULATOR DELOREAN (v2)
## AUTHOR: MATIAS A. BLANA D.
## CHILE, SANTIAGO JULY 2023
## VERSION : MW DWARF SATELLITES INCLUDING MW MULTIPLE POTENTIALS and M31, AND Fornax Cluster potentials, WITH ORBITS WITH COSMIC EXPANSION OPTIONS
## SCRIPT  : CAN LOAD AND STORE INITIAL CONDITIONS OF ASTRONOMICAL OBJECTS TO CALCULATE THEIR ORBITS LATER, AND  
## #######################################################################################################################################
##########################################################################################################################################
from DELOREAN_pylib import *
from DELOREAN_UC_v1 import *
from DELOREAN_FUN_v1 import *

##################################################################
##################################################################
######<b> Define instance/class objects generator (satellites, clusters, etc) <b>#########
# class Galaxy:
#     pass
class Emcee:
    pass

class Galaxy:
    def __init__(self):
        self.emcee = Emcee()
        self.emcee.use = False
    pass

class GalaxyCluster:
    def __init__(self):
        self.emcee = Emcee()
        self.emcee.use = False    
    pass

class StarCluster:
    def __init__(self):
        self.emcee = Emcee()
        self.emcee.use = False    
    pass

class Stars:
    def __init__(self):
        self.emcee = Emcee()
        self.emcee.use = False    
    pass


###################################################################################################
### FUNCTIONS TO GENERATE VALUES FOR THE PARAMETER SPACE EXPLORATION ##############################
def ParamSpaceSats(test=False):
    if test:
        vt_i,vt_f,vt_stp   = 0.0,200.,200.
        nvels_out          = np.arange(vt_i,vt_f+vt_stp,vt_stp)
        Ang_ut_utt         = 0. # 55.
        dAng               = 90.
        vAng_out           = np.arange(0,dAng,dAng) + Ang_ut_utt
    elif not test:
        vt_i,vt_f,vt_stp   = 0.0,10.,0.5
        nvels1             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
        vt_i,vt_f,vt_stp   = 10.,25.,2.
        nvels2             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
        vt_i,vt_f,vt_stp   = 28.,170.,2.
        nvels3             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
        vt_i,vt_f,vt_stp   = 180.,250.,10.
        nvels4             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
        vt_i,vt_f,vt_stp   = 300.,350.,50.
        nvels5             = np.arange(vt_i,vt_f+vt_stp,vt_stp)    
        nvels_out          = np.concatenate((nvels1,nvels2,nvels3,nvels4,nvels5),axis=0)
        Ang_ut_utt         = 0. # 55.
        dAng               = 10.
        vAng0              = np.arange(0,360,dAng) + Ang_ut_utt
        vAng_out           = np.zeros((len(vAng0)+1))
        vAng_out[:6]       = vAng0[:6]
        vAng_out[6]        = 55. 
        vAng_out[7:]       = vAng0[6:]

    nvel = len(nvels_out)
    nang = len(vAng_out)
    npartot   =  nvel*nang
    print(f"norbs_vt={nvel}")
    print(f"norbs_ang={nang}")
    print(f"norbs_tot={npartot}")
    npars = 2
    paramvars = []
    cont = 0
    for i in range(0,nvel):
        for j in range(0,nang):
            paramvars.append([nvels_out[i],vAng_out[j]])            
            cont = cont +1
    return paramvars

def SetICveltanLeoT():
    vt_i,vt_f,vt_stp   = 0.0,10.,0.5
    nvels1             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
    vt_i,vt_f,vt_stp   = 10.,25.,2.
    nvels2             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
    vt_i,vt_f,vt_stp   = 28.,170.,2.
    nvels3             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
    vt_i,vt_f,vt_stp   = 180.,250.,10.
    nvels4             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
    vt_i,vt_f,vt_stp   = 300.,350.,50.
    nvels5             = np.arange(vt_i,vt_f+vt_stp,vt_stp)    
    nvels_out          = np.concatenate((nvels1,nvels2,nvels3,nvels4,nvels5),axis=0)
    Ang_ut_utt         = 0. # 55.
    dAng               = 10.
    vAng0              = np.arange(0,360,dAng) + Ang_ut_utt
    vAng_out           = np.zeros((len(vAng0)+1))
    vAng_out[:6]       = vAng0[:6]
    vAng_out[6]        = 55. 
    vAng_out[7:]       = vAng0[6:]
    print("norbs_vt=",len(nvels_out))
    print("norbs_ang=",len(vAng_out))
    print("norbs_tot=",len(vAng_out)*len(nvels_out))
    return nvels_out,vAng_out

def SetICveltanLeoTNbody():
    nvels_out          = np.array([0.,30.,50.,100.,150.,200.,250.])
    vAng_out           = np.array([50.])
    print("norbs_vt=",len(nvels_out))
    print("norbs_ang=",len(vAng_out))
    print("norbs_tot=",len(vAng_out)*len(nvels_out))
    return nvels_out,vAng_out

def SetICveltanLeoTNbody2():
    nvels_out          = np.array([10.,20.,40.,60.,70.,80.,90.])
    vAng_out           = np.array([50.])
    print("norbs_vt=",len(nvels_out))
    print("norbs_ang=",len(vAng_out))
    print("norbs_tot=",len(vAng_out)*len(nvels_out))
    return nvels_out,vAng_out

def SetICveltanSat():
    vt_i,vt_f,vt_stp   = 0.0,10.,0.5
    nvels1             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
    vt_i,vt_f,vt_stp   = 10.,30.,2.
    nvels2             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
    vt_i,vt_f,vt_stp   = 5.,170.,5.
    nvels3             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
    vt_i,vt_f,vt_stp   = 180.,250.,10.
    nvels4             = np.arange(vt_i,vt_f+vt_stp,vt_stp)
    vt_i,vt_f,vt_stp   = 300.,350.,50.
    nvels5             = np.arange(vt_i,vt_f+vt_stp,vt_stp)    
    nvels_out          = np.concatenate((nvels1,nvels2,nvels3,nvels4,nvels5),axis=0)
    Ang_ut_utt         = 0. # 55.
    dAng               = 45.
    vAng0              = np.arange(0,360,dAng) + Ang_ut_utt
    vAng_out           = vAng0
    print("norbs_vt=",len(nvels_out))
    print("norbs_ang=",len(vAng_out))
    print("norbs_tot=",len(vAng_out)*len(nvels_out))
    return nvels_out,vAng_out


def SetICveltanSat_Test(vt,ang):
    nvels_out          = np.array([vt])
    vAng_out           = np.array([ang])
    
    nvel = len(nvels_out)
    nang = len(vAng_out)
    npartot   =  nvel*nang
    print(f"norbs_vt={nvel}")
    print(f"norbs_ang={nang}")
    print(f"norbs_tot={npartot}")
    npars = 2
    # paramvars = np.zeros((npartot,2))
    paramvars = []
    cont = 0
    for i in range(0,nvel):
        for j in range(0,nang):
            # paramvars[cont,:] = nvels_out[i],vAng_out[j]
            paramvars.append([nvels_out[i],vAng_out[j]])            
            cont = cont +1
    return paramvars

## <b>Functions to calculate coordinates <b> ########################################################
def VlosGSR(self): #Function to convert VlosHe to GSR when Vtan=0
    galaxy_in_icrs      = coord.ICRS(ra=self.ra, dec=self.dec, distance=self.DHe)
    galaxy_in_gal       = galaxy_in_icrs.transform_to(coord.Galactic)
    self.vlosGSR      = self.vlosHe +v_sun_dif_BG16.dot(galaxy_in_gal.data.to_cartesian())/galaxy_in_gal.data.to_cartesian().norm()
#     print(str(self.name)+".vlosGSR=",self.vlosGSR)
    return self.vlosGSR

def Vlos2Coord(galaxy): # Function to convert coords to GC,gal,icrs when Vtan=0
    galaxy_icrs         = coord.ICRS(ra=galaxy.ra, dec=galaxy.dec, distance=galaxy.DHe)
    galaxy_gal          = galaxy_icrs.transform_to(coord.Galactic)
    galaxy_GC           = galaxy_gal.transform_to(coord.Galactocentric_BG16)
    lgal                = galaxy_gal.l.radian  #theta=l
    bgal                = galaxy_gal.b.radian  #phi=b
    thetaGSR            = lgal
    gamma_tilt          = 0.13 ## BG2016 Fig 5. Plane GC and Plane l,b is tilted by gamma_tilt
    bci                 = 0.046 #deg
    zsuni               = 0.25 #kpc
    Rsuni               = 8.2 #kpc
#     gamma_tilt          = (180.-(180.-np.arctan(zsuni/(Rsuni**2-zsuni**2)**0.5))-bci)
    phiGSR              = bgal # - gamma_tilt
    galaxy.VthetaGSR    = 0. #km/s
    galaxy.VphiGSR      = 0. #km/s
    galaxy_GC.VXGSR     = (np.cos(phiGSR)*np.cos(thetaGSR)*galaxy.vlosGSR + 0. + 0.)#*u.km/u.s  CARTESIAN AROUND THE SUN!! 
    galaxy_GC.VYGSR     = (np.cos(phiGSR)*np.sin(thetaGSR)*galaxy.vlosGSR + 0. + 0.)#*u.km/u.s
    galaxy_GC.VZGSR     = (np.sin(phiGSR)*galaxy.vlosGSR - 0.) #*u.km/u.s
    galaxy_GC.VXGSR,galaxy_GC.VYGSR,galaxy_GC.VZGSR = Rotysc(galaxy_GC.VXGSR.value,galaxy_GC.VYGSR.value,galaxy_GC.VZGSR.value,-gamma_tilt)
    galaxy_GC.VXGSR *= u.km/u.s
    galaxy_GC.VYGSR *= u.km/u.s
    galaxy_GC.VZGSR *= u.km/u.s
    galaxy.GC   = coord.Galactocentric(x=galaxy_GC.x, y=galaxy_GC.y,z=galaxy_GC.z, v_x=galaxy_GC.VXGSR, v_y=galaxy_GC.VYGSR, \
                                       v_z=galaxy_GC.VZGSR,galcen_distance=galcen_distance_BG16,galcen_v_sun=v_sun_BG16,z_sun=z_sun_BG16)
    galaxy.icrs = galaxy.GC.transform_to(coord.FK5)
    galaxy.gal  = galaxy.icrs.transform_to(coord.Galactic)
    return galaxy

def Vlos2Vrad(galaxy): # Function to convert coords to GC,gal,icrs when Vtan=0
    galaxy_icrs         = coord.ICRS(ra=galaxy.ra, dec=galaxy.dec, distance=galaxy.DHe)
    galaxy_gal          = galaxy_icrs.transform_to(coord.Galactic)
    galaxy_GC           = galaxy_gal.transform_to(coord.Galactocentric_BG16)
    lgal                = galaxy_gal.l.radian  #theta=l
    bgal                = galaxy_gal.b.radian  #phi=b
    thetaGSR            = lgal
    gamma_tilt          = 0.13 ## BG2016 Fig 5. Plane GC and Plane l,b is tilted by gamma_tilt
    bci                 = 0.046 #deg
    zsuni               = 0.25 #kpc
    Rsuni               = 8.2 #kpc
#     gamma_tilt          = (180.-(180.-np.arctan(zsuni/(Rsuni**2-zsuni**2)**0.5))-bci)
    phiGSR              = bgal # - gamma_tilt
    galaxy.VthetaGSR    = 0. #km/s
    galaxy.VphiGSR      = 0. #km/s
    galaxy_GC.VXGSR     = (np.cos(phiGSR)*np.cos(thetaGSR)*galaxy.vlosGSR + 0. + 0.)#*u.km/u.s  CARTESIAN AROUND THE SUN!! 
    galaxy_GC.VYGSR     = (np.cos(phiGSR)*np.sin(thetaGSR)*galaxy.vlosGSR + 0. + 0.)#*u.km/u.s
    galaxy_GC.VZGSR     = (np.sin(phiGSR)*galaxy.vlosGSR - 0.) #*u.km/u.s
    galaxy_GC.VXGSR,galaxy_GC.VYGSR,galaxy_GC.VZGSR = Rotysc(galaxy_GC.VXGSR.value,galaxy_GC.VYGSR.value,galaxy_GC.VZGSR.value,-gamma_tilt)
    galaxy_GC.VXGSR *= u.km/u.s
    galaxy_GC.VYGSR *= u.km/u.s
    galaxy_GC.VZGSR *= u.km/u.s
    galaxy.GC   = coord.Galactocentric(x=galaxy_GC.x, y=galaxy_GC.y,z=galaxy_GC.z, v_x=galaxy_GC.VXGSR,\
                                       v_y=galaxy_GC.VYGSR,v_z=galaxy_GC.VZGSR,galcen_distance=galcen_distance_BG16,galcen_v_sun=v_sun_BG16,z_sun=z_sun_BG16)
    galaxy.icrs = galaxy.GC.transform_to(coord.FK5)
    galaxy.gal  = galaxy.icrs.transform_to(coord.Galactic)
    
    rad   = (galaxy_GC.x.value**2+galaxy_GC.y.value**2+galaxy_GC.z.value**2)**0.5
    print("r=",rad)
    theta   = np.arctan(galaxy_GC.y.value/galaxy_GC.x.value)
    phi = np.arcsin(galaxy_GC.z.value/rad)
    
    Vrad = (galaxy.GC.x.value*galaxy.GC.v_x.value + galaxy.GC.y.value*galaxy.GC.v_y.value + galaxy.GC.z.value*galaxy.GC.v_z.value)/rad
    
    VXrad = galaxy.GC.x.value/rad*Vrad
    VYrad = galaxy.GC.y.value/rad*Vrad
    VZrad = galaxy.GC.z.value/rad*Vrad
    
    VXtan = -(galaxy.GC.v_x.value-VXrad)
    VYtan = -(galaxy.GC.v_y.value-VYrad)
    VZtan = -(galaxy.GC.v_z.value-VZrad)
    
    VXtan *=u.km/u.s
    VYtan *=u.km/u.s
    VZtan *=u.km/u.s
    
    galaxy.Vtan   = coord.Galactocentric(x=galaxy.GC.x, y=galaxy.GC.y,z=galaxy.GC.z, v_x=VXtan,\
                                         v_y=VYtan, v_z=VZtan,galcen_distance=galcen_distance_BG16,
                                         galcen_v_sun=v_sun_BG16,z_sun=z_sun_BG16)
    galaxy.Vtan.icrs = galaxy.Vtan.transform_to(coord.FK5)
    galaxy.Vtan.gal  = galaxy.Vtan.transform_to(coord.Galactic)    
    
    galaxy.Vrad   = coord.Galactocentric(x=galaxy.GC.x, y=galaxy.GC.y,z=galaxy.GC.z, v_x=VXrad*u.km/u.s, v_y=VYrad*u.km/u.s,\
                                         v_z=VZrad*u.km/u.s,galcen_distance=galcen_distance_BG16,galcen_v_sun=v_sun_BG16,z_sun=z_sun_BG16)
    galaxy.Vrad.icrs = galaxy.Vrad.transform_to(coord.FK5)
    galaxy.Vrad.gal  = galaxy.Vrad.transform_to(coord.Galactic)  
    
    return galaxy

def Vtan2GC(galaxy): # Function to obtain GC vectors perp to Rlos (or Vlos)
    # Calculating L vector and unitary velocity vectors nu_ut nu_utt tangential to vradial GC
    x_GC     = galaxy.GC.x.value #kpc
    y_GC     = galaxy.GC.y.value #kpc
    z_GC     = galaxy.GC.z.value #kpc
    r_GC     = (x_GC**2+y_GC**2+z_GC**2)**0.5 #kpc
    nrv_GC   = np.array([x_GC,y_GC,z_GC])/r_GC
    x_GSR     = x_GC + (dsun.value**2-Zsun.value**2)**0.5 #kpc
    y_GSR     = y_GC #kpc
    z_GSR     = z_GC -Zsun.value #kpc
    r_GSR     = (x_GSR**2+y_GSR**2+z_GSR**2)**0.5 #kpc
#     print("r_GC=",r_GC,"[kpc]"
#     print("r_GC=",r_GC/288.,"[rvir] with Rvir=",288,"[kpc]"
#     print("r_GSR=",r_GSR,"[kpc]"
    nrv_GSR   = np.array([x_GSR,y_GSR,z_GSR])/r_GSR
    galaxy.rGSR = np.array([x_GSR,y_GSR,z_GSR])*u.kpc
    vx    = galaxy.GC.v_x.value  #km/s
    vy    = galaxy.GC.v_y.value  #km/s
    vz    = galaxy.GC.v_z.value  #km/s
    vtot  = (vx**2+vy**2+vz**2)**0.5
    nv    = np.array([vx,vy,vz])/vtot
#     print('nrv_GSR=',nrv_GSR
#     print('nv_GSR_GC=',nv
    L    = np.cross(nrv_GC,nv)
    nL   = L/(L[0]**2+L[1]**2+L[2]**2)**0.5
    nu_tt= nL.copy()
    u_t   = np.cross(nL,nrv_GSR)
    nu_t    = u_t/(u_t[0]**2+u_t[1]**2+u_t[2]**2)**0.5
    galaxy.nut  =  nu_t
    galaxy.nutt =  nu_tt
#     print("nu_t=",nu_t
#     print("nu_tt=",nu_tt,"=nL" 
    return galaxy

def GetCoord(self,*argv):
    if self.usepropmot:
#         galaxy.pm_ra_cosdec = argv[0]
#         galaxy.pm_dec       = argv[1]
        self.icrs         = coord.ICRS(ra=self.ra, dec=self.dec, distance=self.DHe, radial_velocity=self.vlosHe, pm_ra_cosdec=self.pm_ra_cosdec, pm_dec=self.pm_dec)
        self.gal          = self.icrs.transform_to(coord.Galactic)    
        self.GC           = self.gal.transform_to(coord.Galactocentric_BG16)
        self.vlosGSR      = VlosGSR(self)
        self              = Vtan2GC(self)
    else:
        self.vlosGSR      = VlosGSR(self)
        self              = Vlos2Coord(self)
        self              = Vtan2GC(self)
#     print(str(self.name)+".vlosGSR=",self.vlosGSR
#     print(str(self.name)+".GC=",self.GC
    return self

def GetIC(galaxy):
    galpos    = np.array([galaxy.GC.x.value,galaxy.GC.y.value,galaxy.GC.z.value])
    galvelos  = np.array([galaxy.GC.v_x.value,galaxy.GC.v_y.value,galaxy.GC.v_z.value])
    gal_nu_t  = galaxy.nut
    gal_nu_tt = galaxy.nutt
    return galpos,galvelos,gal_nu_t,gal_nu_tt

# def GetCoord(galaxy,*argv):
#     if galaxy.usepropmot:
# #         galaxy.pm_ra_cosdec = argv[0]
# #         galaxy.pm_dec       = argv[1]
#         galaxy.icrs         = coord.ICRS(ra=galaxy.ra, dec=galaxy.dec, distance=galaxy.DHe, radial_velocity=galaxy.vlosHe, pm_ra_cosdec=galaxy.pm_ra_cosdec, pm_dec=galaxy.pm_dec)
#         galaxy.gal          = galaxy.icrs.transform_to(coord.Galactic)    
#         galaxy.GC           = galaxy.gal.transform_to(coord.Galactocentric_BG16)
#         galaxy.vlosGSR      = VlosGSR(galaxy)
#         galaxy              = Vtan2GC(galaxy)
#     else:
#         galaxy.vlosGSR      = VlosGSR(galaxy)
#         galaxy              = Vlos2Coord(galaxy)
#         galaxy              = Vtan2GC(galaxy)
#     print(str(galaxy.name)+".vlosGSR=",galaxy.vlosGSR
#     print(str(galaxy.name)+".GC=",galaxy.GC
#     return galaxy

def vtangsrerr(self):
    self.usepropmot = True
    self            = GetCoord(self)
    vtotXYZ0         = np.array([self.GC.v_x.value, self.GC.v_y.value, self.GC.v_z.value])
    pmRAs     = self.pm_ra_cosdec.value
    pmDEC     = self.pm_dec.value
    print("vtotXYZ0=",vtotXYZ0)
    self.usepropmot = False
    self            = GetCoord(self)
    vlosXYZ0        = np.array([self.GC.v_x.value, self.GC.v_y.value, self.GC.v_z.value])
    print("vlosXYZ0=",vlosXYZ0)
    vtanXYZ0         = vtotXYZ0-vlosXYZ0
    NvtanXYZ0        = np.linalg.norm(vtanXYZ0)
    print("vtanXYZ0=",NvtanXYZ0)
    pmRAsE    = self.pm_ra_cosdec_err.value
    pmDECE    = self.pm_dec_err.value
    nsteps    = 100
    pmRA_arr  = np.arange(pmRAs-pmRAsE,pmRAs+pmRAsE+2.*pmRAsE/nsteps,2.*pmRAsE/nsteps)
    pmDEC_arr = np.arange(pmDEC-pmDECE,pmDEC+pmDECE+2.*pmDECE/nsteps,2.*pmDECE/nsteps)
    NvtanXYZmin = 1.E20
    NvtanXYZmax = 0.
    self.usepropmot = True
    for i in range(0,len(pmRA_arr)):
        for j in range(0,len(pmDEC_arr)):
            self.pm_ra_cosdec = pmRA_arr[i]*u.mas/u.yr
            self.pm_dec       = pmDEC_arr[j]*u.mas/u.yr
            self            = GetCoord(self)
            vtotXYZ         = np.array([self.GC.v_x.value, self.GC.v_y.value, self.GC.v_z.value])    
            vtanXYZ         = vtotXYZ-vlosXYZ0
            NvtanXYZ        = np.linalg.norm(vtanXYZ)
            NvtanXYZmin     = np.amin(np.array([NvtanXYZ,NvtanXYZmin]))
            NvtanXYZmax     = np.amax(np.array([NvtanXYZ,NvtanXYZmax]))
        print("N=",(i+1)*len(pmDEC_arr) ,"of",len(pmRA_arr)*len(pmDEC_arr))
    print("vtanXYZ0=",NvtanXYZ0,"+",NvtanXYZmax-NvtanXYZ0,"-",NvtanXYZ0-NvtanXYZmin)
    self.pm_ra_cosdec = pmRAs*u.mas/u.yr
    self.pm_dec       = pmDEC*u.mas/u.yr
    self              = GetCoord(self)
    self.vtanGSR      = NvtanXYZ0*u.km/u.s
    self.vtanGSRmax   = (NvtanXYZmax)*u.km/u.s
    self.vtanGSRmin   = (NvtanXYZmin)*u.km/u.s
    return

def vtangsrerr_v2(self): #this version uses the correlation between mu_alpha* and mu_dec
    self.usepropmot = True
    self            = GetCoord(self)
    vtotXYZ0         = np.array([self.GC.v_x.value, self.GC.v_y.value, self.GC.v_z.value])
    pmRAs     = self.pm_ra_cosdec.value
    pmDEC     = self.pm_dec.value
    print("vtotXYZ0=",vtotXYZ0)
    self.usepropmot = False
    self            = GetCoord(self)
    vlosXYZ0        = np.array([self.GC.v_x.value, self.GC.v_y.value, self.GC.v_z.value])
    print("vlosXYZ0=",vlosXYZ0)
    vtanXYZ0         = vtotXYZ0-vlosXYZ0
    NvtanXYZ0        = np.linalg.norm(vtanXYZ0)
    print("vtanXYZ0=",NvtanXYZ0)
    pmRAsE    = self.pm_ra_cosdec_err.value
    pmDECE    = self.pm_dec_err.value
    pmCORR    = self.pm_corr.value
    
    sigma2_muas=pmRAsE**2
    sigma2_mud =pmDECE**2
    
    COV      = pmCORR*pmRAsE*pmDECE
    sigma2_U = 0.5*(sigma2_muas+sigma2_mud+((sigma2_muas+sigma2_mud)**2-4.*COV**2)**0.5)
    sigma2_V = 0.5*(sigma2_muas+sigma2_mud-((sigma2_muas+sigma2_mud)**2-4.*COV**2)**0.5)
    tantheta    = (sigma2_muas-sigma2_mud)/COV
    sigma_U = sigma2_U**0.5
    sigma_V = sigma2_V**0.5
    nsteps    = 100
    U_arr  = np.arange(-sigma_U, sigma_U+2.*sigma_U/nsteps,2.*sigma_U/nsteps)
    V_arr = np.arange(-sigma_V, sigma_V+2.*sigma_V/nsteps,2.*sigma_V/nsteps)
    
    print("problem, sign error RA, DEC, two solutions??")
    
    NvtanXYZmin = 1.E20
    NvtanXYZmax = 0.
    self.usepropmot = True
    for i in range(0,len(U_arr)):
        for j in range(0,len(V_arr)):
            pmRA_arr  = pmRAsE + U_arr[i]
            pmDEC_arr = pmDECE + V_arr[j]
            self.pm_ra_cosdec = pmRA_arr*u.mas/u.yr
            self.pm_dec       = pmDEC_arr*u.mas/u.yr
            self            = GetCoord(self)
            vtotXYZ         = np.array([self.GC.v_x.value, self.GC.v_y.value, self.GC.v_z.value])    
            vtanXYZ         = vtotXYZ-vlosXYZ0
            NvtanXYZ        = np.linalg.norm(vtanXYZ)
            NvtanXYZmin     = np.amin(np.array([NvtanXYZ,NvtanXYZmin]))
            NvtanXYZmax     = np.amax(np.array([NvtanXYZ,NvtanXYZmax]))
        print("N=",(i+1)*len(pmDEC_arr) ,"of",len(pmRA_arr)*len(pmDEC_arr))
    print("vtanXYZ0=",NvtanXYZ0,"+",NvtanXYZmax-NvtanXYZ0,"-",NvtanXYZ0-NvtanXYZmin)
    self.pm_ra_cosdec = pmRAs*u.mas/u.yr
    self.pm_dec       = pmDEC*u.mas/u.yr
    self              = GetCoord(self)
    self.vtanGSR      = NvtanXYZ0*u.km/u.s
    self.vtanGSRmax   = (NvtanXYZmax)*u.km/u.s
    self.vtanGSRmin   = (NvtanXYZmin)*u.km/u.s
    return


##################################################################
######<b>Loading and defininf Satellite and properties<b>#########
# List satellites for the paper: 
# LeoT        DGC=414kpc
# Cetus       DGC=756kpc (far away gas-less backsplash?) 
# Eridanus II DGC=365kpc
# Phoenix I   DGC=419kpc
# Leo I ?     DGC=273kpc
# Leo II ?    DGC=227kpc
# CVen I ?    DGC=211kpc
# Piscis II ? DGC=182kpc
# Leo V
# Leo IV

# List of satellites to set!
# setLeoT   = True
setLeoT   = False
setLeoV   = False
setLeoIV  = False
setHer    = False
setCet    = False
setMW     = False
setIC10   = False
setIC1613 = False
setEriII  = False
setPhxI   = False
setNGC6569b=False

setM31 = False
setMW  = False

# import pickle
# saving data as class
# fnameIC='Data/Data_IC/satellites_IC_v2.pkl'
fnameIC='Data/Data_Inputs/Data_IC/satellites_IC_v3.pkl'
cwd = os.getcwd()
print("directory ",cwd)
# os.system("rm "+fnameIC)

if setLeoT:
    LeoT              = Galaxy()
    LeoT.name         = "LeoT"
    LeoT.coordtype    = "Blana2020"
    print(f"Setting: {LeoT.name}")
    LeoT.Mst          = 2.E+5 # Msun assumed in Adams+18
    LeoT.Lv           = 1.41E+05 #Lsun from de Jong 2008 Mv,tot=-8mag
    LeoT.MHI          = 4.1E+5 #Msun # Adams et al 2018 but! DHe=420kpc! Need Correction! fc=409/420=0.9738095238095238 ?
    LeoT.MHe          = 1.35E+5 #Msun # Adams et al 2018 but! distance=420kpc!!
    LeoT.Mgas         = 5.45E5 #=LeoT.MHI+LeoT.MHe
    ## Measured quantities ***************************************
    ## Equatorial coordinates SPATIAL, Distance and vlos *****************
    LeoT.ra           = coord.Angle((9., 34., 53.5), unit='hourangle') #C12
    LeoT.dec          = coord.Angle((17.,3.,4.), unit=u.deg)  #C12
    LeoT.DHeI07       = 420.*u.kpc #kpc Heliocentric (Irwin et al 2007 [I07])
    LeoT.DHeSG7       = 417.*u.kpc #kpc Heliocentric (Simon & Geha 2007 [I07])
    LeoT.DHeC12       = (409.)*u.kpc #+29 -27kpc Heliocentric (Clementini et al 2012 [C12])
    LeoT.vlosHe       = 38.1*u.km/u.s  # Simon & Geha 2007 average line-of-sight stellar velocity
    LeoT.usepropmot   = False
    ## ************************************************************
    LeoT.DHe          = LeoT.DHeC12
    LeoT              = GetCoord(LeoT)
    print(f'test={LeoT.GC}')
    # LeoT.vlosGSR        = -59.4977112195 #km/s # old with vsun wrong 237km/s
    # kpc2a             = np.arctan(1./LeoT.DHe.value)*180./np.pi*3600.
    # pc2a              = kpc2a/1000.
    # a2pc              = 1./pc2a
    # f_area1to2        = (np.arctan(1./LeoT.DHe.value)/np.arctan(1./LeoT.DHeI07.value))**2.
    # a2kpc             = 1./kpc2a
    # kpc2a             = 1./pc2a
    # kpc2aA18          = np.arctan(1./LeoT.DHeI07.value)*180./np.pi*3600.
    # pc2aA18           = kpc2aA18/1000.
    # a2kpcA18          = 1./kpc2aA18
    # kpc2aA18          = 1./pc2aA18
    # print(LeoT)
    with open(fnameIC, 'wb') as output:
          pickle.dump(LeoT, output, pickle.HIGHEST_PROTOCOL)

# LeoT.vlosGSR= -65.9193721543 km / s
# r_GC= 413.920671763 [kpc]
# r_GC= 1.43722455473 [rvir] with Rvir= 288 [kpc]
# r_GSR= 409.0 [kpc]
# nu_t= [-0.80419409  0.36399235 -0.46987385]
# nu_tt= [-0.0570412  -0.83416057 -0.54856398] =nL
# LeoT.vlosGSR= -65.9193721543 km / s
# LeoT.GC= <Galactocentric Coordinate (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
#     (266.4051, -28.936175)>, galcen_distance=8.2 kpc, galcen_v_sun=(11., 248., 7.3) km / s, z_sun=25.0 pc, roll=0.0 deg): (x, y, z) in kpc
#     (-250.14785843, -169.09021362, 283.13401614)
#  (v_x, v_y, v_z) in km / s
#     (39.13408255, 27.25258745, -45.51025912)>
#     <FK5 Coordinate (equinox=J2000.000): (ra, dec, distance) in (deg, deg, kpc)
#     (143.72292212, 17.0511058, 409.)
#  (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
#     (-0.01509574, -0.11533449, 38.0640209)>
##########################################################################################
# LeoT.vlosGSR= -65.9193721543 km / s
# r_GC= 413.920671763 [kpc]
# r_GC= 1.43722455473 [rvir] with Rvir= 288 [kpc]
# r_GSR= 409.059164189 [kpc]
# nu_t= [-0.80407776  0.3640738  -0.47000981]
# nu_tt= [-0.0570412  -0.83416057 -0.54856398] =nL
# LeoT.vlosGSR= -65.9193721543 km / s
# LeoT.GC= <Galactocentric Coordinate (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
#     (266.4051, -28.936175)>, galcen_distance=8.2 kpc, galcen_v_sun=(11., 248., 7.3) km / s, z_sun=25.0 pc, roll=0.0 deg): (x, y, z) in kpc
#     (-250.14785843, -169.09021362, 283.13401614)
#  (v_x, v_y, v_z) in km / s
#     (39.13408255, 27.25258745, -45.51025912)>
#########################################################################################
# val = 38.1*u.km/u.s
# val = 39.6*u.km/u.s

# LeoT.vlosHe=val
# VlosGSR(LeoT)
# err = 0.1*u.km/u.s
# LeoT.vlosHe=val + err
# VlosGSR(LeoT)
# err = -0.1*u.km/u.s
# LeoT.vlosHe=val + err
# VlosGSR(LeoT)

if setHer:
    Her              = Galaxy()
    Her.name         = "Hercules"
    print(f"Setting: {Her.name}")    
    Her.DHe          = 138.*u.kpc
    Her.ra           = coord.Angle((16., 31., 2.), unit='hourangle')
    Her.dec          = coord.Angle((12.,47.,13.83), unit=u.deg)
    Her.vrHe         = 45.*u.km/u.s
    Her.pm_ra_cosdec = 0.*u.mas/u.yr
    Her.pm_dec       = 0.*u.mas/u.yr
    Her_icrs         = coord.ICRS(ra=Her.ra, dec=Her.dec, distance=Her.DHe, radial_velocity=Her.vrHe, pm_ra_cosdec=Her.pm_ra_cosdec, pm_dec=Her.pm_dec)
    Her_gal          = Her_icrs.transform_to(coord.Galactic)
    v_sun            = coord.Galactocentric.galcen_v_sun.to_cartesian()
    cart_data        = Her_gal.data.to_cartesian()
    v_proj           = v_sun.dot(cart_data)/cart_data.norm()
    rv_gsr           = Her_icrs.radial_velocity + v_proj
    print('Vr_GSR=',rv_gsr)
    with open(fnameIC, 'wb') as output:
        pickle.dump(Her, output, pickle.HIGHEST_PROTOCOL)
        
        
if setCet:
    Cet              = Galaxy()
    Cet.name         = "Cetus"
    Cet.usepropmot   = False
    Cet.DHe          = 755*u.kpc # 755+-23 McConnachie et al. 2005
    Cet.ra           = coord.Angle((00., 26., 10.5), unit='hourangle')
    Cet.dec          = coord.Angle((-11.,2.,32.), unit=u.deg)
    # Cet.vlosHe         = -83.9*u.km/u.s # Kirby+14
    Cet.vlosHe       = -78.9*u.km/u.s # -78.9+1.7-1.6 kms-1 Taibi+18
    Cet              = GetCoord(Cet)
    Cet.icrs
    print(Cet.vlosGSR)
    with open(fnameIC, 'wb') as output:
        pickle.dump(Cet, output, pickle.HIGHEST_PROTOCOL)
    
# Cetus.vlosGSR= -14.8533463202 km / s
# r_GC= 755.524316756 [kpc]
# r_GC= 2.62334832207 [rvir] with Rvir= 288 [kpc]
# r_GSR= 755.0 [kpc]
# nu_t= [-0.9980996  -0.01333344  0.06016144]
# nu_tt= [0.00464452 0.9572549  0.28920838] =nL
# Cetus.vlosGSR= -14.8533463202 km / s
# Cetus.GC= <Galactocentric Coordinate (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
#     (266.4051, -28.936175)>, galcen_distance=8.2 kpc, galcen_v_sun=(11., 248., 7.3) km / s, z_sun=25.0 pc, roll=0.0 deg): (x, y, z) in kpc
#     (-54.5915266, 218.17787326, -721.27330053)
#  (v_x, v_y, v_z) in km / s
#     (0.86943235, -4.29225849, 14.19304412)>
# <FK5 Coordinate (equinox=J2000.000): (ra, dec, distance) in (deg, deg, kpc)
#     (6.54375523, -11.04222034, 755.)
#  (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
#     (0.03737741, -0.05566031, -78.86961128)>

# arr1=np.array([-54.5915266, 218.17787326, -721.27330053])
# arr=np.array([-379.71375482, 612.66794081, -281.95992153])-arr1
# sum(arr**2)**0.5


if setPhxI:
    PhxI              = Galaxy()
    PhxI.name         = "PhoenixI"
    PhxI.ra           = coord.Angle((01., 51., 6.3), unit='hourangle')
    PhxI.dec          = coord.Angle((-44.,26.,41.), unit=u.deg)
    PhxI.DHe          = 409.*u.kpc #+-23Battaglia2012  #415.*u.kpc #+-19 McConnachie2005   #419*u.kpc # McConnachie2012
    PhxI.vlosHe       = -21.2*u.km/u.s# +- 1.0Kacharov2017
    PhxI.usepropmot   = False
    # PhxI.usepropmot   = True
    PhxI.pm_ra_cosdec = 0.079*u.mas/u.yr # Fritz+18 from Gaia
    PhxI.pm_dec       = -0.049*u.mas/u.yr # Fritz+18 from Gaia    
    PhxI.pm_ra_cosdec_err = (0.099+0.04)*u.mas/u.yr # Fritz+18 from Gaia
    PhxI.pm_dec_err       = (0.12+0.04)*u.mas/u.yr # Fritz+18 from Gaia
    # PhxI.pm_ra_cosdec_err = (0.099)*u.mas/u.yr # Fritz+18 from Gaia
    # PhxI.pm_dec_err       = (0.12)*u.mas/u.yr # Fritz+18 from Gaia
    # PhxI.vtanGSR
    PhxI              = GetCoord(PhxI)
    # PhxI.vtanGSR    = 67.993441*u.km/u.s
    # PhxI.vtanGSRmax = 473.96937*u.km/u.s
    # PhxI.vtanGSRmin = 2.4158523*u.km/u.s
    PhxI.icrs
    # vlosXYZ = np.array([-1.82917991, 41.9529042, 109.07651144])
    # vtotXYZ = np.array([-40.36683888, 94.0996993, 88.615118])
    # vtanXYZ = vtotXYZ-vlosXYZ
    # print(" ")
    # print(" ")
    # print("vlosGSR=",PhxI.vlosGSR)
    # print("vlosXYZ=",np.linalg.norm(vlosXYZ))
    # print("vtotXYZ=",np.linalg.norm(vtotXYZ))
    # print("vtanXYZ=",np.linalg.norm(vtanXYZ))
    print(PhxI.vlosGSR)
    print(PhxI.GC)
    with open(fnameIC, 'wb') as output:
        pickle.dump(PhxI, output, pickle.HIGHEST_PROTOCOL)
# PhxI.vtanGSRmin
# vtangsrerr(PhxI) # vtanXYZ0= 67.9934409984874 + 405.97593182329734 - 65.57758870895746

if setEriII:
    EriII              = Galaxy()
    EriII.name         = "EridanusII"
    EriII.DHe          = 366*u.kpc #  /366kpc Crnojevic+16  // 365 DGC !! Fritz+18
    EriII.ra           = coord.Angle((03., 44., 20.1), unit='hourangle') # Crnojevic+16
    EriII.dec          = coord.Angle((-43.,32.,01.7), unit=u.deg) # Crnojevic+16
    EriII.vlosHe       = 75.6*u.km/u.s # +- 1.3 +-2.0 Li+16
    EriII.pm_ra_cosdec       = 0.159*u.mas/u.yr # Fritz+18 from Gaia
    EriII.pm_dec             = 0.372*u.mas/u.yr # Fritz+18 from Gaia
    EriII.pm_ra_cosdec_err   = (0.292+0.053)*u.mas/u.yr # Fritz+18 from Gaia
    EriII.pm_dec_err         = (0.34+0.053)*u.mas/u.yr # Fritz+18 from Gaia
    EriII.usepropmot   = False
    # EriII.usepropmot   = True
    EriII              = GetCoord(EriII)
    print(EriII.GC)
    # EriII.icrs
    # EriII.vtanGSR      = 761.3711810814606*u.km/u.s
    # EriII.vtanGSRmax   = (1601.263)*u.km/u.s
    # EriII.vtanGSRmin   = (71.351417)*u.km/u.s
    # vlosXYZ = np.array([16.34776326, 44.76409644, 60.32960066])
    # vtotXYZ = np.array([-689.21797584, 330.37701854, 43.34955594])
    # vtanXYZ = vtotXYZ-vlosXYZ
    # print(" ")
    # print(" ")
    # print("vlosGSR=",EriII.vlosGSR)
    # print("vlosXYZ=",np.linalg.norm(vlosXYZ))
    # print("vtotXYZ=",np.linalg.norm(vtotXYZ))
    # print("vtanXYZ=",np.linalg.norm(vtanXYZ))
    # print(((EriII.GC.v_x**2 +EriII.GC.v_y**2 +EriII.GC.v_z**2)-(76.8*u.km/u.s)**2)**0.5)
    print(EriII.vlosGSR)
    # vtangsrerr(EriII)
    with open(fnameIC, 'wb') as output:
        pickle.dump(EriII, output, pickle.HIGHEST_PROTOCOL)
        
if setIC1613:
    IC1613              = Galaxy()
    IC1613.name         = "IC1613"
    IC1613.DHe          = 755.*u.kpc  # +-42 Mcchonacchie 2012
    IC1613.ra           = coord.Angle((1.,4.,47.79), unit='hourangle') # ned.ipac
    IC1613.dec          = coord.Angle((2.,7.,4.), unit=u.deg)  # ned.ipac
    IC1613.vlosHe       = -234.137943*u.km/u.s  # +/- 0.899378 ned.ipac
    IC1613.usepropmot   = False
    IC1613              = GetCoord(IC1613)
    print(IC1613.icrs)
    print(IC1613.vlosGSR)
    with open(fnameIC, 'wb') as output:
        pickle.dump(IC1613, output, pickle.HIGHEST_PROTOCOL)
        
if setIC10:
    IC10              = Galaxy()
    IC10.name         = "IC10"
    IC10.DHe          = 794.*u.kpc  # +-44 Mcchonacchie 2012
    IC10.ra           = coord.Angle((0.,20.,17.34), unit='hourangle') # ned.ipac
    IC10.dec          = coord.Angle((59,18.,13.6), unit=u.deg)  # ned.ipac
    IC10.vlosHe       = -348.059092*u.km/u.s  # +/- 0.899378ned.ipac
    IC10.usepropmot   = False
    IC10              = GetCoord(IC10)
    print(IC10.icrs)
    print(IC10.vlosGSR)
    with open(fnameIC, 'wb') as output:
        pickle.dump(IC10, output, pickle.HIGHEST_PROTOCOL)

if setLeoV:
    LeoV              = Galaxy()
    LeoV.name         = "LeoV"
    LeoV.DHe          = 173.*u.kpc  # Riley+18
    LeoV.ra           = coord.Angle((172.784,0.,0.), unit='hourangle') #  Riley+18
    LeoV.dec          = coord.Angle((2.222,0.,0.), unit=u.deg)  # Riley+18
    LeoV.vlosHe       = 172.1*u.km/u.s  # Riley+18
    LeoV.pm_ra_cosdec = -0.097*u.mas/u.yr # Fritz+18
    LeoV.pm_dec       = -0.628*u.mas/u.yr    # Fritz+18
    LeoV.usepropmot   = False
    LeoV              = GetCoord(LeoV)
    # # LeoV.DHe          = 175.*u.kpc  #deJong+10  
    # LeoV.ra           = coord.Angle((11.,31.,8.4), unit='hourangle') #deJong+10 
    # LeoV.dec          = coord.Angle((2.,12.,57.), unit=u.deg)  #deJong+10

if setLeoIV:
    LeoIV              = Galaxy()
    LeoIV.name         = "LeoIV"
    LeoIV.DHe          = 154.*u.kpc  # Riley+18
    LeoIV.ra           = coord.Angle((173.233), unit='hourangle')  # Riley+18
    LeoIV.dec          = coord.Angle((-0.540), unit=u.deg)  # Riley+18
    LeoIV.vlosHe       = 132.3*u.km/u.s  # Riley+18
    pm_ra_cosdec       = -0.59*u.mas/u.yr  # Fritz+18
    pm_dec             = -0.449*u.mas/u.yr # Fritz+18
    LeoIV              = GetCoord(LeoIV)

## LeoIV.DHe          = 154.*u.kpc  #deJong+10
## LeoIV.ra           = coord.Angle((11., 32., 58.6), unit='hourangle')  #deJong+10
## LeoIV.dec          = coord.Angle((0.,33.,6.), unit=u.deg)  #deJong+10


if setMW:
    MW           =  Galaxy()
    MW.name      = "MilkyWay"
    ### FOR DICE_v3 #######
    ## now with MW.M200=1.1E12
    ## got MW.Mvir=1.2784124999140897 E12
    MW.M200      = 1.1E+12 ## ??
    MW.Msum      = 1.2784134000539514E12 ## resulting from DICE when MW.MDM=1.3E12
    MW.M200b     = 1.3E12
    MW.Msumb     = 1.469387E+12
    MW.M200c     = 1.14E12
    MW.MDM       = MW.M200c
    # MW.f_MDMvir= MW.MDM200/MW.MDMsum
    # MW.MDMvir  = MW.f_MDMvir*1.E12
    MW.disk      = 6.6E+10 #Ms
    MW.bulge     = 1.12E+10
    MW.dyn       = MW.MDM + MW.disk + MW.bulge
    print("MW_Mdyn=",MW.dyn)
    # print("MW_fMDMvir=",MW.f_MDMvir)
    # print("MW_MDMvir=",MW.MDMvir)
    print("f_MW.MDM=",MW.MDM/MW.dyn)
    print("f_MW.disk=",MW.disk/MW.dyn)
    print("f_MW.bulge=",MW.bulge/MW.dyn)
    print("z/R=",0.3/2.5) #kpc)
    
    ## DICE_MW_v4 #########
    MW.M200d = 1.1350E12
    MW.MDM   = MW.M200d
    # MW.f_MDMvir= MW.MDM200/MW.MDMsum
    # MW.MDMvir  = MW.f_MDMvir*1.E12
    MW.disk    = 5.50E+10 #Ms
    MW.bulge   = 0.50E+10
    MW.dyn     = MW.MDM + MW.disk + MW.bulge
    print("MW_Mdyn=",MW.dyn)
    # print("MW_fMDMvir=",MW.f_MDMvir)
    # print("MW_MDMvir=",MW.MDMvir)
    print("f_MW.MDM=",MW.MDM/MW.dyn)
    print("f_MW.disk=",MW.disk/MW.dyn)
    print("f_MW.bulge=",MW.bulge/MW.dyn)
    print("zh/Rd=",0.3/2.5) #kpc/kpc)
    print("rh_bulge=",0.3,"kpc")
    
if setM31:
    M31a              = Galaxy()
    M31a.name         = "M31a"
    M31a.Mvir         = 1.1E12*u.Msun
    M31a.rvir         = 258.*(340.*0.3/102)**(-1./3)*(M31a.Mvir.value/1.E12)**(1./3)*u.kpc
    M31a.ra           = coord.Angle((10.68333), unit=u.deg)
    M31a.dec          = coord.Angle((41.26917), unit=u.deg)
    M31a.DHe          = 770.*u.kpc  # vdM12a
    M31a.vlosHe       = -301*u.km/u.s #km/s Courteau & van den Bergh 1999
    M31a.pm_ra_cosdec      = 65.*1E-3*u.mas/u.yr  # +-18.E-3mas/yr van der Marel+2019 Gaia
    M31a.pm_dec            = -57.*1E-3*u.mas/u.yr # +-15.E-3mas/yr van der Marel+2019 Gaia
    M31a.usepropmot   = True
    M31a              = GetCoord(M31a)
    print("|V_M31MW|=",(M31a.GC.v_x.value**2+M31a.GC.v_y.value**2+M31a.GC.v_z.value**2)**0.5)
    print("VlosCGR=",M31a.vlosGSR)
    ##M31.ra           = coord.Angle((0., 42., 44.3), unit='hourangle')
    ##M31.dec          = coord.Angle((41.,16.,9.), unit=u.deg)
    ##M31.DHe          = 778.*u.kpc  #
    ##M31_GC
    ##outputM31 = commah.run('WMAP1',zi=0.,Mi=M31.Mvir.value,z=[0.0])
    ##outputM31['c']
    ###########################################################
    M31b              = Galaxy()
    M31b.name         = "M31"
    M31b.Mvir         = 1.1E12*u.Msun
    M31b.rvir         = 258.*(340.*0.3/102)**(-1./3)*(M31b.Mvir.value/1.E12)**(1./3)*u.kpc
    M31b.ra           = coord.Angle((10.68333), unit=u.deg)
    M31b.dec          = coord.Angle((41.26917), unit=u.deg)
    M31b.DHe          = 770.*u.kpc  # vdM12a
    M31b.vlosHe       = -301*u.km/u.s #km/s
    M31b.pm_ra_cosdec      = 65.*1E-3*u.mas/u.yr  # +-18.E-3mas/yr van der Marel+2019 Gaia
    M31b.pm_dec            = -57.*1E-3*u.mas/u.yr # +-15.E-3mas/yr van der Marel+2019 Gaia
    M31b.usepropmot   = False
    M31b              = GetCoord(M31b)
    # fnameIC='Data/Data_IC/M31_IC_v2.pkl'
    fnameIC='Data/Data_IC/M31_IC_v3.pkl'
    os.system("rm "+fnameIC)
    with open(fnameIC, 'wb') as output:
        pickle.dump(M31a, output, pickle.HIGHEST_PROTOCOL)
        pickle.dump(M31b, output, pickle.HIGHEST_PROTOCOL)

# Adding tangential velocity component in the plane of Leo T - MW - M31 <b>
# posM31  = np.array([M31_galcen.x.value,M31_galcen.y.value,M31_galcen.z.value])
# posLeoT = np.array([LeoT_galcen.x.value,LeoT_galcen.y.value,LeoT_galcen.z.value])
# print("posM31=",posM31)
# print("posLeoT=",posLeoT)

# from numpy import linalg as LA
# # numpy.linalg.norm(x, ord=None, axis=None, keepdims=False)
# r1 = posM31
# r2 = posLeoT

# r3 = np.cross(r1,r2)
# r4 = np.cross(r3,r1)
# n4 = r4/LA.linalg.norm(r4)
# print("n4=",n4)

# Distance_LeoT_M31 = np.linalg.norm(np.array([-250.1900275, -169.09021362, 283.18542387])-np.array([-383.74265383, 619.02654549, -284.81140451]))
# print("Distance_LeoT_M31=",Distance_LeoT_M31,"kpc")
# "F_MW/F_M31=",(1./419**2)/(1./Distance_LeoT_M31**2)

############################################################################
## STAR CLUSTERS ###########################################################
if setNGC6569b:
    NGC6569b              = StarCluster()  # http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=NGC+6569
    NGC6569b.name         = "NGC6569b"
    NGC6569b.ra           = coord.Angle((18., 13., 38.88), unit='hourangle') #2006MNRAS.365.1357D
    NGC6569b.dec          = coord.Angle((-31.,49.,35.2), unit=u.deg)
    NGC6569b.DHe          = 10.9 *u.kpc # +- Harris et al 1996 http://physwww.physics.mcmaster.ca/~harris/mwgc.dat
    # NGC6569.vlosHe       = -28.1*u.km/u.s #+- 5.6 1996AJ....112.1487H
    NGC6569b.vlosHe       = -48.8*u.km/u.s #+- 5.6  Johnson et al 2018
    NGC6569b.usepropmot   = True
    NGC6569b.pm_ra_cosdec = -3.05*u.mas/u.yr # +- 0.56 2013A&A...558A..53K
    NGC6569b.pm_dec       = -1.85*u.mas/u.yr # +- 0.56 2013A&A...558A..53K
    NGC6569b              = GetCoord(NGC6569b)
    NGC6569b.icrs
    ## E. Vasiliev 2019 PAPER!! 
    NGC6569              = StarCluster()  
    NGC6569.name         = "NGC6569"
    NGC6569.ra           = coord.Angle((273.412,0.,0.), unit=u.deg) 
    NGC6569.dec          = coord.Angle((-31.827,0.,0.), unit=u.deg)
    NGC6569.DHe          = 10.9 *u.kpc
    NGC6569.vlosHe       = -49.83*u.km/u.s 
    NGC6569.usepropmot   = True
    NGC6569.pm_ra_cosdec = -4.109*u.mas/u.yr 
    NGC6569.pm_dec       = -7.267*u.mas/u.yr
    NGC6569              = GetCoord(NGC6569)
    NGC6569.icrs
    # saving data as class
    fnameIC='Data/Data_IC/starclusters_IC_v1.pkl'
    os.system("rm "+fnameIC)
    with open(fnameIC, 'wb') as output:
        pickle.dump(NGC6569, output, pickle.HIGHEST_PROTOCOL)
        
        
#########################################################################
### FORNAX CLUSTER ######################################################
