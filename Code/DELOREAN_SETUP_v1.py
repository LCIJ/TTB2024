#########################################################################################################################################
## ######################################################################################################################################
## CODE: ORBIT CALCULATOR DELOREAN (v2)
## AUTHOR: MATIAS A. BLANA D.
## CHILE, SANTIAGO JULY 2023
## VERSION : MW DWARF SATELLITES INCLUDING MW MULTIPLE POTENTIALS and M31, AND Fornax Cluster potentials, WITH ORBITS WITH COSMIC EXPANSION OPTIONS
## SCRIPT  : SETUP FUNCTIONS TO SETUP MAIN POTENTIALS, MASSES, LOAD OR PRE-COMPUTE MW AND M31 ORBITS, PRE-COMPUTED COSMOLOGICAL TABLES a(t), H(t), etc
## #######################################################################################################################################
##########################################################################################################################################
from DELOREAN_pylib import *
from DELOREAN_UC_v1 import *    ## IMPORT CONSTANTS AND UNITS 
from DELOREAN_IC_v1 import *    ## IMPORT INITIAL CONDITIONS GENERATOR FOR DELOREAN INTEGRATOR
from DELOREAN_FUN_v1 import *
# from DELOREAN_FUN_v1 import rvirf
path_inputs = 'Data/Data_Inputs/'

####################################################################################
#### SETUP BUILDER AND GATHERER ####################################################
####################################################################################
class MySetup:
    SETUPS = {}
    def __init__(self, setup_name=None, case=None):
        self.register_setups()
        if setup_name is not None:
            self.set_setup(setup_name, case)
        else:
            self.list_setups()

    @classmethod
    def register_setups(cls):
        # Register setups if not already registered
        for setup_class in SetupBase.__subclasses__():
            setup_name = getattr(setup_class, 'setup_name', None)
            if setup_name:
                cls.SETUPS[setup_name] = setup_class

    def set_setup(self, setup_name, case=None):
        setup_class = self.SETUPS.get(setup_name)
        if setup_class:
            self.setup = setup_class(case=case)
        else:
            raise ValueError(f"Invalid setup name: {setup_name}")

    def list_setups(self):
        print("Make your setup, or choose among the available Setups:", list(self.SETUPS.keys()))

    def call_setup_function(self, function_name, x):
        if hasattr(self, 'setup') and self.setup is not None:
            setup_function = getattr(self.setup, function_name, None)
            if callable(setup_function):
                setup_function(x)
            else:
                print(f"{function_name} not found or not callable in the current setup.")
        else:
            print("No setup selected.")

class SetupBase:
    setup_name = None
    @classmethod
    def register_setup(cls):
        MySetup.register_setups()  # Ensure setups are registered
        MySetup.register_setup(cls.setup_name, cls)

################################################################################
### DEFINE SETUPS BELOW ########################################################
################################################################################

################################################################################
#### SETUP MILKY WAY SATELLITES ################################################
################################################################################
class Setup1(SetupBase):
    setup_name = 'Setup:Blana2020'
    def __init__(self, case=None):
        super().__init__()
        self.name      = self.setup_name
        self.case      = case
        self.MRcfun_v2 = self.MRcfun_v2
        self.accase    = self.accase

    def MRcfun_v2(self,objclass,t_in):
        case = objclass.setup.case
        global GLfrac_Mdisk, GLfrac_Mbulge,  GLMdisk, GLMbulge
        if case=="case2" or case=="case4" or case=="case6" or case=="case2_cosmo" or case=="case4_cosmo":
            Mviri  = f_Mvir_t(t_in)
            Rviri  = f_Rvir_t(t_in)
            ci     = f_c_t(t_in)
            GLfrac_Mdisk  = GLMdisk/Mviri[0]
            GLfrac_Mbulge = GLMbulge/Mviri[0]           
        elif case=="case2b" or case=="case4b":  
            Mviri  = f_Mvir_t2b(t_in)
            Rviri  = f_Rvir_t2b(t_in)
            ci     = f_c_t2b(t_in)
            GLfrac_Mdisk  = GLMdisk/Mviri[0]
            GLfrac_Mbulge = GLMbulge/Mviri[0]               
        elif case=="case2a" or case=="case4a":    
            Mviri  = f_Mvir_t2a(t_in)
            Rviri  = f_Rvir_t2a(t_in)
            ci     = f_c_t2a(t_in)
            GLfrac_Mdisk  = GLMdisk/Mviri[0]
            GLfrac_Mbulge = GLMbulge/Mviri[0]               
        elif case=="case1" or case=="case3" or case=="case3fut" or case=="case5" or case=="case1nbody" or case=="case1_cosmo" or case=="case3_cosmo" or case=="case1_kepler":
            Mviri  = 1.3E+12*np.ones(len(t_in)) #Ms
            Rviri  = ((Mviri)**(1./3)/(4.*np.pi/3.*Planck15.critical_density(0).to(u.solMass/u.kpc**3).value*Planck15.Om(0)*f_deltavir_t(0))**(1./3)/f_a_t(0))*np.ones(len(t_in)) #kpc
            ci     = 8.63263214*np.ones(len(t_in))
            GLfrac_Mdisk  = GLMdisk/Mviri[0]
            GLfrac_Mbulge = GLMbulge/Mviri[0]               
        elif case=="case1a" or case=="case3a":
            Mviri  = 1.0E+12*np.ones(len(t_in)) #Ms
            Rviri  = ((Mviri)**(1./3)/(4.*np.pi/3.*Planck15.critical_density(0).to(u.solMass/u.kpc**3).value*Planck15.Om(0)*f_deltavir_t(0))**(1./3)/f_a_t(0))*np.ones(len(t_in))*np.ones(len(t_in)) #kpc
            ci     = 8.84952006*np.ones(len(t_in))     
            GLfrac_Mdisk  = GLMdisk/Mviri[0]
            GLfrac_Mbulge = GLMbulge/Mviri[0]               
        elif case=="case1b" or case=="case3b":
            Mviri  = 1.6E+12*np.ones(len(t_in)) #Ms
            Rviri  = ((Mviri)**(1./3)/(4.*np.pi/3.*Planck15.critical_density(0).to(u.solMass/u.kpc**3).value*Planck15.Om(0)*f_deltavir_t(0))**(1./3)/f_a_t(0))*np.ones(len(t_in))*np.ones(len(t_in)) #kpc
            ci     = 8.43986264*np.ones(len(t_in))  
            GLfrac_Mdisk  = GLMdisk/Mviri[0]
            GLfrac_Mbulge = GLMbulge/Mviri[0]               
        elif case=="case1old":
            Mviri  = 1.3E+12*np.ones(len(t_in)) #Ms
            Rviri  = 282.*np.ones(len(t_in)) #kpc
            ci     = 11.5*np.ones(len(t_in))     
            GLfrac_Mdisk  = GLMdisk/Mviri[0]
            GLfrac_Mbulge = GLMbulge/Mviri[0]
        else:
            Mviri  = 1.3E+12*np.ones(len(t_in)) #Ms
            Rviri  = 288.*np.ones(len(t_in)) #kpc
            ci     = 8.43986264*np.ones(len(t_in))
            GLfrac_Mdisk  = GLMdisk/Mviri[0]
            GLfrac_Mbulge = GLMbulge/Mviri[0]   
        return Mviri,Rviri,ci
        
    def accase(self,objclass,case,ti,dto,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci):
        case = objclass.setup.case
        if case=="case1" or case=="case1a" or case=="case1b" or case=="case2" or case=="case2a" or case=="case2b":
            acc = acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
        elif case=="case3" or case=="case3a" or case=="case3b" or case=="case4" or case=="case4a" or case=="case4b" or case=="case3_cosmo" or case=="case4_cosmo":
            acc = acchalonfw_M31(ri-riPOT2)+acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
        elif case=="case3_cosmo" or case=="case4_cosmo":
            acc = acchalonfwnostars_MW_ac(ri-riPOT2,Mviri,Rviri,ci)+acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
        elif case=="case3" or case=="case3a" or case=="case3b" or case=="case4" or case=="case4a" or case=="case4b":
            acc = acchalonfw_M31(ri-riPOT2)+acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
        elif case=="case5" or case=="case6" or case=="case5fut":
            acc = acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri) + accdynfr(dto,vi-viPOT1,ri-riPOT1,Mviri,Rviri,ci)
        elif case=="case1fut":
            acc = acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
        elif case=="caseCOtest":
            acc = acchalonfwnostars_MW_ac(ri-riPOT1)
        elif case=="case1nbody":
            acc = accnbody(ri-riPOT1)
        elif case=="case1_cosmo" or case=="case2_cosmo": 
            acc = acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
        elif case=="case1_kepler":
            acc = acchalokeplernostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)
        elif case=="case1_burk":
            acc = acc_burkert(ri-riPOT1)
        elif case=="case2_burk":
            acc = acc_burkert2(ri-riPOT1)    
        elif case=="case1_nfw":
            acc = acc_nfw(ri-riPOT1)             
        elif case=="case_bar1":
            # acc = accbar(ri,ti)   +acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)  #  rvin-riPOT1
            acc = accbarLM(ri,ti) +accdisc_MW_ac(ri-riPOT1,Mviri) +acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)  #  rvin-riPOT1   
        elif case=="case_bar2":
            acc = accbarLM(ri,ti) +acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)  #  rvin-riPOT1               
        return acc

################################################################################
#### SETUP FORNAX CLUSTER ######################################################
################################################################################
class Setup2(SetupBase):
    setup_name = 'Setup:Fornax'
    def __init__(self, case=None):
        super().__init__()
        self.name      = self.setup_name
        self.case      = case
        self.MRcfun_v2 = self.MRcfun_v2
        self.accase    = self.accase
        
    def MRcfun_v2(self,objclass,t_in):
        case = objclass.setup.case
        if case=="case0":
            Mviri  = 1.0E+14*np.ones(len(t_in)) #Ms # Mastroprieto 2021* and ref there
            Rviri  = 978.0*np.ones(len(t_in)) #kpc
            ci     = 8.15*np.ones(len(t_in))      
        elif case=="case1":
            Mviri  = 7.0E+13*np.ones(len(t_in)) #Ms # Mastroprieto 2021* and ref there
            Rviri  = 1063.29160735*np.ones(len(t_in)) #kpc
            ci     = 5.78330723*np.ones(len(t_in))              
        elif case=="case2":
            Mviri  = f_Mvir_tF2(t_in)
            Rviri  = f_Rvir_tF2(t_in)
            ci     = f_c_tF2(t_in)
        return Mviri,Rviri,ci
    
    def accase(self,objclass,case,ti,dto,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci):      
        if case=="case1":
            acc = acc_Fornax_nfw(ri-riPOT1)  
        elif case=="case2":
            acc = acc_Fornax_nfw_ac(ri-riPOT1,Mviri,Rviri,ci)
        return acc
    
################################################################################
#### SETUP MALIN-1 GALAXY ######################################################
################################################################################
class Setup3(SetupBase):
    setup_name = 'Setup:Malin1'
    def __init__(self, case=None):
        super().__init__()
        self.name      = self.setup_name
        self.case      = case
        self.MRcfun_v2 = self.MRcfun_v2
        self.accase    = self.accase

    def MRcfun_v2(self,objclass,t_in):
        case = objclass.setup.case
        if case=="case1" or case=="case1i":
            Mviri  = 2.0E+12*np.ones(len(t_in)) #Ms # Mastroprieto 2021* and ref there
            Rviri  = rvirf(Mviri,0.) #kpc
            ci     = 3.*np.ones(len(t_in))
        elif case=="case2":
            Mviri  = 2.0E+12*np.ones(len(t_in)) #Ms # Mastroprieto 2021* and ref there
            Rviri  = rvirf(Mviri,0.) #kpc
            ci     = 20.*np.ones(len(t_in))    
        elif case=="case3":
            Mviri  = 4.0E+12*np.ones(len(t_in)) #Ms # Mastroprieto 2021* and ref there
            Rviri  = rvirf(Mviri,0.) #kpc
            ci     = 10.*np.ones(len(t_in))
        elif case=="case4":
            Mviri  = 4.0E+12*np.ones(len(t_in)) #Ms # Mastroprieto 2021* and ref there
            Rviri  = rvirf(Mviri,0.) #kpc
            ci     = 10.*np.ones(len(t_in))            
        return Mviri,Rviri,ci

    def accase(self,objclass,case,ti,dto,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci):      
        if case=="case1":
            acc = acc_Malin1(ri-riPOT1,Mviri,Rviri,ci)
        elif case=="case1i":
            acc = acc_Malin1(ri-riPOT1,Mviri,Rviri,ci,i=35.,PA=90.)
        elif case=="case2":
            acc = acc_Malin1(ri-riPOT1,Mviri,Rviri,ci)
        elif case=="case3":
            acc = acc_Malin1(ri-riPOT1,Mviri,Rviri,ci)     #  rvin-riPOT1
        elif case=="case4":
            acc = acc_nfw1(ri,Mviri,Rviri,ci) + acc_plumer1(ri,Mpl=1e10,rpl=1e0) + acc_miyamotonagai1(ri,Md=1e11,Rd=1e1,zh=1e0,incl=0.,PA=90.)
        return acc #kpc/Myr^2


################################################################################
#### EXAMPLE SETUP #############################################################
################################################################################
class Setup_Example(SetupBase):
    setup_name = 'Setup:Example'
    def __init__(self, case=None):
        super().__init__()
        self.name      = self.setup_name
        self.case      = case
        self.MRcfun_v2 = self.MRcfun_v2
        self.accase    = self.accase

    def MRcfun_v2(self,objclass,t_in):
        case = objclass.setup.case
        if case=="case1":
            Mviri  = 1.0E+13*np.ones(len(t_in)) #Ms 
            Rviri  = 700*np.ones(len(t_in)) #kpc
            ci     = 10*np.ones(len(t_in))              
        elif case=="case2":
            Mviri  = f_Mvir_tF2(t_in)
            Rviri  = f_Rvir_tF2(t_in)
            ci     = f_c_tF2(t_in)
        # elif case=='case3':
        # ...
        return Mviri,Rviri,ci
    
    def accase(self,objclass,case,ti,dto,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci):      
        if case=="case1":
            acc = acc_nfw1(ri,Mviri,Rviri,ci)  
        elif case=="case2":
            acc = acc_nfw1(ri,Mviri,Rviri,ci)
        # elif case=='case3':
        # ...
        return acc

##############################################################################
##############################################################################

#################################################################################################################################
## FUNCTION TO CONVERT INPUT USER COORDINATES (AS VARIABLE OR FIXED PARAMETERS) TO CARTESIAN COORDINATES FOR THE ORBIT INTEGRATOR
def CoordIN2CoordORB_v2(self,params):
    if self.coordtype == 'Blana2020': #(uses Vtan and Angle to Supergalactic) # Vtan, Angle to x,y,z,vx,vy,vz,time (paper)    
        nvels_in,vAng_in = params # 2 parameters: |Vlos|, and angle (direction) of satellite.
        objpos,objvels,obj_nut,obj_nutt = GetIC(self)
        ut      = np.array([-nvels_in * np.cos(vAng_in/180.*np.pi)])  ## perp to L
        utt     = np.array([nvels_in* np.sin(vAng_in/180.*np.pi)])    ## parall to L
        objvels = objvels + utt*obj_nutt + ut*obj_nut
        objcoords= [objpos[:],objvels[:]]
        
    if self.coordtype == 'Equatorial': #(uses RA,DEC, etc from astropy to Supergalactic)
        RA            = params[0]*u.deg
        DEC           = params[1]*u.deg
        dist          = params[2]*u.kpc
        vlos          = params[3]*u.km/u.s
        pmra,pmcosdec = params[4]*u.mas/u.yr,params[5]*u.mas/u.yr
        # obj_icrs      = coord.ICRS(ra=RA, dec=DEC, distance=dist, radial_velocity=vlos, pm_ra_cosdec=pmra, pm_dec=pmcosdec)
        obj_icrs      = coord.SkyCoord(ra=RA, dec=DEC, distance=dist, radial_velocity=vlos, pm_ra_cosdec=pmra, pm_dec=pmcosdec, frame='icrs')        
        # obs parameters are defined in DELOREAN_UC_v1.py, astropy conv.: https://docs.astropy.org/en/stable/generated/examples/coordinates/plot_galactocentric-frame.html#sphx-glr-generated-examples-coordinates-plot-galactocentric-frame-py
        obj_GC        = obj_icrs.transform_to(coord.Galactocentric_BG16) 
        # obj_GC        = obj_icrs.transform_to(coord.Galactocentric_def)
        try:
            if self.usebaumgardtcoord:
                fact = -1.
            else:
                fact = 1.
        except:
            fact=1    
        objcoords     = [np.array([fact*obj_GC.x.value,obj_GC.y.value,obj_GC.z.value]),np.array([fact*obj_GC.v_x.value,obj_GC.v_y.value,obj_GC.v_z.value])]
        
    if self.coordtype == 'Extragal': #(converts RA,DEC to ANGLES in respect to a reference frame (host) to Cartesian on the sky (FLAT SKY APPROXIMATION))    
        Host          = CoordHost(self.coordtype.host)
        RA            = params[0]
        DEC           = params[1]
        zik           = params[2]
        vx,vy,vz      = params[3],params[4],params[5]
        c1            = coord.SkyCoord(RA, DEC, frame='icrs', unit=u.deg)
        c1_to_host    = c1.transform_to(Host.skyoffset_frame())
        xi, eta       = -c1_to_host.lon.deg, c1_to_host.lat.deg
        try:
            redshift  = self.host.redshift
            # asec2kpc_host = astropy.cosmology.funcs.angular_diameter_distance(redshift, cosmo=None)
            asec2kpc_host = (Planck15.angular_diameter_distance(0.1).to(u.kpc)/(u.rad.to(u.arcsec))).value
        except:
            distance      = (host.distance.to(u.kpc)).value
            asec2kpc_host = 1./(distance) /np.pi*648000. # 1kpc = asec2kpc[arcsec]: = 1kpc/D/pi*180*60*60 # VALID ONLY FOR LOW REDSHIFT!! z<<1, see https://en.wikipedia.org/wiki/Angular_diameter_distance
        xik, etak     = xi*deg2as*asec2kpc_host, eta*deg2as*asec2kpc_host
        objcoords= [np.array([x,y,z]),np.array([vx,vy,vz])]  
        
    if self.coordtype == 'Cartesian': #(direct conversion)
        x,y,z,vx,vy,vz = params
        objcoords      = [np.array([x,y,z]),np.array([vx,vy,vz])]
        
    if self.coordtype == 'Malin1': #(direct conversion)
        x,y,vz     = self.paramfix
        z,vt,vtang = params
        vy,vx      = vt*np.sin(vtang*np.pi/180.), vt*np.cos(vtang*np.pi/180.)
        objcoords  = [np.array([x,y,z]),np.array([vx,vy,vz])]        
        
    return objcoords

###################################################################################################
### Coordinates for host systems: #################################################################
def CoordHost(hostname): # Malin 1 Optical Source: http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=malin+1&submit=SIMBAD+search
    if hostname=='Malin1':
        ra, dec, distance = '12h36m59.3463193464s', '+14d19m49.164305988s', 366. * u.Mpc
        coords  = SkyCoord(ra, dec, frame='icrs', distance = distance)
    return coords
##################################################################################################


##########################################################
##### CIRCULAR VELOCITIES FROM POTENTIAL ACCELERATION ####
##########################################################
def Vcirc(rvin,objclass,ti):
    case = objclass.setup.case
    # if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
    #     Mviri,Rviri,ci,ai,Hi = objclass.setup.MRcfun(objclass,ti)
    # else:
    Mviri,Rviri,ci = objclass.setup.MRcfun_v2(objclass,ti)
    # accase(case,ti,dto,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci)    
    acc = objclass.setup.accase(objclass,case,ti,0,rvin,np.zeros((rvin.shape)),0,0,0,Mviri,Rviri,ci)
    accn = (np.sum(acc**2,axis=0))**0.5
    #using ac = Vc^2/r =>  Vc = (ac*r)**0.5
    rvec  = rvin
    y,x   = rvec[1], rvec[0]
    # theta = np.degrees(np.arctan2(y,x))
    R     = (x**2 + y**2)**0.5 #(np.sum(rvec**2,axis=0))**0.5
    vc    = (abs(accn)*R)**0.5*kpcMyr2kms
    vx,vy = -vc*y/R , vc*x/R
    vz    = 0.
    return vx,vy,vz

##############################################################################
########## FUNCTIONS: ACCELERATIONS FOR DIFFERENT POTENTIAL PROFILES #########
##############################################################################
def acc_logpot1(rvin,rvir,vo=150,d=0.5):
    vo    = vo*kms2kpcMyr # km/s*((kpc/Myr)/(km/s)) kpc/Myr 186.
    rin   = (np.dot(rvin,rvin))**0.5
    acc   = -vo**2*(1./(rin**2+d**2))*rvin #(kpc/Myr)^2 kpc^-1
    accVir= -(vo*kpcMyr2kms)**2*(1./(rvir**2+d**2))*rvir #(kpc/Myr)^2 kpc^-1
    Mvir  = -accVir*rvir**2/G2
    if (rin>=rvir): acc=-G2*Mvir/(rin*rin)*(rvin/rin)*kms2kpcMyr**2 #(kpc/Myr)^2/kpc
    return acc #kpc/Myr^2

def acc_kepler1(rvin,Mvir,eps=1.):
    rin  = (np.dot(rvin,rvin))**0.5
    acc  = -G4*Mvir/(rin**2+eps**2)**(3./2)*rvin #(kpc/Myr)^2/kpc
    return acc #kpc/Myr^2

def acc_nfw1(rvin,Mvir,rvir,cin):
    rs   = rvir/cin
    A    = -G4*(Mvir)/(np.log(1.+cin)-cin/(1.+cin))
    rin  = (np.dot(rvin,rvin))**0.5
    acc  = -A*(rvin/rin*(1./(rin*(rin+rs)))-np.log(1.+rin/rs)/rin**2*rvin/rin)
    if (rin>=rvir): acc =-G4*Mvir/(rin*rin)*(rvin/rin) #(kpc/Myr)^2/kpc
    return acc #kpc/Myr^2

def acc_burkert1(rvin,rho_in,rs_in):
    rinf  = 0.0001 # =0.1pc # Mass undefined/diverges for r=0
    rin  = (np.dot(rvin,rvin))**0.5
    if rin<=rinf: rin = rinf
    MbR   = 4.*np.pi*rho_in*rs_in**3.*1./2*(1./2.*np.log(1.+rin**2/rs_in**2) + np.log(1.+rin/rs_in) - np.arctan(rin/rs_in)) #- Cte
    acc   = -G4*MbR/(rin**3)*rvin + plum_LeoT_ac(rvin)
    return acc #kpc/Myr^2

def acc_plumer1(rvin,Mpl=1e10,rpl=1e0):
    rin  = (np.dot(rvin,rvin))**0.5
    acc   = -G4*Mpl/(rin**2+rpl**2)**1.5*rvin
    return acc #kpc/Myr^2

def acc_miyamotonagai1(rvin,Md=1e11,Rd=1e1,zh=1e0,incl=0.,PA=90.):
    PAang   = PA-90.
    inclang = incl
    xi   = rvin[0]
    yi   = rvin[1]
    zi   = rvin[2]
    if (inclang!=0.) or (PAang!=0.):
        xir,yir,zir    = Rotz_v2(xi,yi,zi,PAang)
        xirr,yirr,zirr = Rotx_v2(xir,yir,zir,inclang)
        xi,yi,zi       = xirr,yirr,zirr
        ri             = (xi**2 + yi**2 + zi**2)**0.5
        Ri             = (xi**2 + yi**2)**0.5
        A              = Ri**2 + (Rd + (zi**2 + zh**2)**0.5)**2
        B              = Rd + (zi**2 + zh**2)**0.5        
        acc            = -G4*Md/A**(3./2)*np.array([xi,yi,B*zi/(zi**2 + zh**2)**0.5])
        accx_r,accy_r,accz_r    = Rotx_v2(acc[0],acc[1],acc[2],-inclang)
        accx_rr,accy_rr,accz_rr = Rotz_v2(accx_r,accy_r,accz_r,-PAang)
        acc = np.array([accx_rr,accy_rr,accz_rr])
    else: # avoids rotations in case PA,incl = 0 to make it faster
        ri             = (xi**2 + yi**2 + zi**2)**0.5
        Ri             = (xi**2 + yi**2)**0.5
        A              = Ri**2 + (Rd + (zi**2 + zh**2)**0.5)**2
        B              = Rd + (zi**2 + zh**2)**0.5        
        acc            = -G4*Md/A**(3./2)*np.array([xi,yi,B*zi/(zi**2 + zh**2)**0.5])        
    return acc #kpc/Myr^2

####################################################################
### TESTING BARRED POTENTIAL BY MONARI ET AL 2016 ##################
def accbarLM(rvin,ti,Mbar=1.0E10,a=3.,b=1.5,c=2.,Omegab=40.,barangle0=27.): # Long & Murali 1992 ApJ 397,44
    # a,b,c = 3.,1.5, 2.  # bar length (bar major axis) kpc, (bar major axis disky scale length // bar minor axis)kpc, (scale height) kpc 
    xi    = rvin[0]
    yi    = rvin[1]
    zi    = rvin[2]
    if (Omegab!=0.) or (barangle0!=0.):
        Angletime = -barangle0 # (by definition +27deg bar angle MW is clock-wise)
        if Omegab!=0.:
            Periodbar = 1./(Omegab*kms2kpcMyr)*2.*np.pi
            Angletime = Angletime + Periodbar*ti # degree    
        xi,yi,zi = Rotz_v2(xi,yi,zi,Angletime)
        
    cf    = (c**2 + zi**2)**0.5
    bf    = (b + cf)
    y2f   = yi**2 + bf**2
    
    Tp    = ((a+xi)**2 + y2f)**0.5
    Tm    = ((a-xi)**2 + y2f)**0.5
    
    ax    = -2.*G4*Mbar*xi/(Tp*Tm)*(1./(Tp+Tm))
    ay    = -G4*Mbar*yi/2./(Tp*Tm)*(1./y2f)*(Tp+Tm - 4.*xi**2/(Tp+Tm))
    az    = -G4*Mbar*zi/2./(Tp*Tm)*(1./y2f)*(Tp+Tm - 4.*xi**2/(Tp+Tm))*bf/cf

    if (Omegab!=0.) or (barangle0!=0.):
        ax,ay,az = Rotz_v2(ax,ay,az,-Angletime) # restore to non-rotated frame

    acc  = np.array([ax,ay,az])
    return acc #kpc/Myr^2


####################################################################
### TESTING BARRED POTENTIAL BY DEHNEN 2000 (not working yet) ######
def potbarv1(xi,yi,zi,ti,Af,Rb,Omegabi):
    Ri2  = xi**2 + yi**2
    Ri   = Ri2**0.5
    ri2  = Ri2 + zi**2
    ri   = ri2**0.5
    phi  = np.arctan2(yi, xi) # radians
    
    if Omegabi==0.:Tb   = 1E40
    if Omegabi!=0.:Tb   = 2.*np.pi/Omegabi
    # Af  = 1.0
    Ab  = Af
    if ri>=Rb: rf = -(Rb/ri)**3.
    if ri<Rb:  rf = (ri/Rb)**3. - 2.
    # pot = Ab*np.cos(2.*(phi-Omegabi*ti))*(Ri2/ri2)*rf
    pot = Ab*np.cos(2.*(phi))*(Ri2/ri2)*rf    
    return pot #

def potbar(xyzi,ti,Af,Rb,Omegabi):
    xi,yi,zi = xyzi
    Ri2  = xi**2 + yi**2
    Ri   = Ri2**0.5
    ri2  = Ri2 + zi**2
    ri   = ri2**0.5
    phi  = np.arctan2(yi, xi) # radians
    
    if Omegabi==0.:Tb   = 1E40
    if Omegabi!=0.:Tb   = 2.*np.pi/Omegabi
    # Af  = 1.0
    Ab  = Af
    if ri>=Rb: rf = -(Rb/ri)**3.
    if ri<Rb:  rf = (ri/Rb)**3. - 2.
    # pot = Ab*np.cos(2.*(phi-Omegabi*ti))*(Ri2/ri2)*rf
    pot = Ab*np.cos(2.*(phi))*(Ri2/ri2)*rf    
    return pot #

def accbarv0(rvin,ti):
    # xi, yi, zi = sp.symbols('xi yi zi')
    # ti, Afi, Rb, Omega = sp.symbols('ti, Afi, Rb, Omega', real=True, positive=False)
    
    Afi     = 1.0
    Rb      = 5. #kpc
    Omegabi = 0.
    xi   = rvin[0]
    yi   = rvin[1]
    zi   = rvin[2]
    
    # pot = potbar(xi,yi,zi,ti,Afi,Rb,Omegab)

    # Define the scalar function f with the constant parameters
    # f = x**2 + y**3 * z**2 + sp.sin(omega * time)
    
    # Compute the gradient vector
    # gradpotbar = [sp.diff(pot, xi), sp.diff(pot, yi), sp.diff(pot, zi)]
    
    # pot  = potbar(xi,yi,zi,ti,Afi,Rb,Omegab)
    # gradpotbar = np.gradient(pot,xi,yi,zi)
    xyzi = [xi,yi,zi]
    gradpotbar = partial_derivatives_bar(potbar, xyzi, ti, Afi, Rb, Omegabi)
    acc  = -1.*gradpotbar
    # Ri2  = xi**2 + yi**2
    # Ri   = Ri2**0.5
    # ri2  = Ri2 + zi**2
    # ri   = ri**0.5
    
    return acc #kpc/Myr^2

def U(r,Rb):
    if r>=Rb: Ur = -(r/Rb)**-3
    if r<Rb:  Ur =  (r/Rb)**3 - 2.
    return Ur

def factgradU(r,Rb):
    if r>=Rb: Ur = 3.*(r/Rb)**-4/Rb
    if r<Rb:  Ur = 3.*(r/Rb)**2/Rb
    return Ur

def gradphi0(xi,yi):
    i     = np.array([1.,0.,0.])
    j     = np.array([0.,1.,0.])    
    if xi==0.: xi=yi/1.E20
    gradphi = 1./(1.+(yi/xi)**2)*(j/xi - yi/xi**2*i)                            
    return gradphi

def gradphi(xi,yi):
    Ri2   = xi**2 + yi**2 
    i     = np.array([1.,0.,0.])
    j     = np.array([0.,1.,0.])    
    gradphi = xi/Ri2*j - yi/Ri2*i                            
    return gradphi

def accbar(rvin,ti):
    alpha,Vo,Ro,Rb = 0.1,210.*kms2kpcMyr,8.,4.
    Ao      = alpha*Vo**2/3.*(Ro/Rb)**3
    Afi     = 1.0
    Omegabi = 1./100000. #1./200. #1/Myr
    phib    = 0.
    xi      = rvin[0]
    yi      = rvin[1]
    zi      = rvin[2]
    Ri2     = xi**2 + yi**2
    Ri      = Ri2**0.5
    ri2     = Ri2 + zi**2
    ri      = ri2**0.5 
    phi     = np.arctan2(yi,xi) #np.pi
    if phi<0:phi=2.*np.pi+phi
    # print('phi=',phi)
    gammab  = 2.*(phi - phib - Omegabi*ti)
    cosb    = np.cos(gammab)
    gradr   = np.array([xi,yi,zi])/ri
    gradR   = np.array([xi,yi,0.])/Ri
    gradRr  = gradR/ri - Ri/ri**2*gradr
    gradU   = factgradU(ri,Rb) * gradr
    gradcosb = -np.sin(gammab) * gradphi(xi,yi)*2.
    gradpotbar = Ao * ( (cosb*U(ri,Rb)*2.*(Ri/ri)*(gradRr) ) + ( cosb*(Ri/ri)**2*(gradU) )+ U(ri,Rb)*(Ri/ri)**2*(gradcosb) )
    acc  = -1.*gradpotbar
    
    return acc #kpc/Myr^2

# import numpy as np
from scipy.optimize import approx_fprime

def scalar_function(xyz, omega, time):
    x, y, z = xyz
    return x**2 + omega*y + np.sin(time)*z

def partial_derivatives_bar(scalar_fun, xyz, ti, Af, Rb, Omegabi):
    """
    Calculate the partial derivatives of the scalar function with respect to x, y, and z.

    Parameters:
        scalar_function (callable): The scalar function of x, y, z and other parameters.
        xyz (array-like): The values of x, y, and z at which to compute the derivatives.
        omega (float): A constant parameter.
        time (float): A constant parameter.

    Returns:
        array-like: The partial derivatives with respect to x, y, and z at the specified point.
    """
    grad_func = lambda xyz: scalar_fun(xyz, ti, Af, Rb, Omegabi)
    partials = approx_fprime(xyz, grad_func, epsilon=1e-6)
    return partials


######################################################################
########## POTENTIALS FOR THE MILKY WAY ##############################
######################################################################
# Define some global parameters
# GLMWMvir = 1.3E+12

#GLMbulge   = 1.12E+10  # Paczy´nski (1990) used in Blana15
global GLMbulge
GLMbulge = 0.5E+10  # BO2015, Bovy2015 (0-25%CB)
# GLfrac_Mbulge = GLMbulge/GLMWMvir
# GLfrac_Mbulge = 0.

#GLMd   = 6.6E+10 #Ms # from Bovy2015??
global GLMdisk
GLMdisk  = 5.5E+10  #Ms BO15 = stellar disk 5E10Msun + gas disk 5E9Msun
# GLfrac_Mdisk = GLMdisk/GLMWMvir
# GLfrac_Mdisk = 0.

def acchalolog_MW(rvin):
    rvir  = 200.
    d     = 0.5
    vo    = 150*kms2kpcMyr # km/s*((kpc/Myr)/(km/s)) kpc/Myr 186.
    rin   = (np.dot(rvin,rvin))**0.5
    acc   = -vo**2*(1./(rin**2+d**2))*rvin #(kpc/Myr)^2 kpc^-1
    accVir= -(vo*kpcMyr2kms)**2*(1./(rvir**2+d**2))*rvir #(kpc/Myr)^2 kpc^-1
    Mvir  = -accVir*rvir**2/G2
    if (rin>=rvir): acc=-G2*Mvir/(rin*rin)*(rvin/rin)*kms2kpcMyr**2 #(kpc/Myr)^2/kpc
    return acc #kpc/Myr^2

def acchalonfwnostars_MW_ac(rvin,Mvir,rvir,cin):
    rs   = rvir/cin
    Mdmg = (1.0 - GLfrac_Mbulge - GLfrac_Mdisk)*Mvir  #Ms    
    A    = -G4*(Mdmg)/(np.log(1.+cin)-cin/(1.+cin))
    rin  = (np.dot(rvin,rvin))**0.5
    acc  = -A*(rvin/rin*(1./(rin*(rin+rs)))-np.log(1.+rin/rs)/rin**2*rvin/rin)
    if (rin>=rvir): acc =-G4*Mdmg/(rin*rin)*(rvin/rin) #(kpc/Myr)^2/kpc
    return acc #kpc/Myr^2

def acchalokeplernostars_MW_ac(rvin,Mvir,rvir,cin):
    eps  =  10. #kpc
    rin  = (np.dot(rvin,rvin))**0.5
    acc  = -G4*Mvir/(rin**2+eps**2)**(3./2)*rvin #(kpc/Myr)^2/kpc
    return acc #kpc/Myr^2

def accdisc_MW_ac(rvin,Mviri):
    zh   = 0.3 #kpc #Ms BO15
    Rd   = 2.5 #kpc #Ms BO15
#     Md   = 6.6E+10 #Ms # from Bovy2015??
    Md   = GLfrac_Mdisk*Mviri #Ms BO15 = stellar disk+ gas disk
    return acc_miyamotonagai1(rvin,Md=Md,Rd=Rd,zh=zh,incl=0.,PA=90.)  #kpc/Myr^2

# def accdisc_MW_ac(rvin,Mviri):
#     zh   = 0.3 #kpc #Ms BO15
#     Rd   = 2.5 #kpc #Ms BO15
# #     Md   = 6.6E+10 #Ms # from Bovy2015??
#     Md   = GLfrac_Mdisk*Mviri #Ms BO15 = stellar disk+ gas disk
#     acc  = acc_miyamotonagai1(rvin,Md=Md,Rd=Rd,zh=zh,incl=0.,PA=90.) #acc_miyamotonagai1(rvin,Md=1e11,Rd=1e1,zh=1e0,incl=0.,PA=90.)
#     xi   = rvin[0]
#     yi   = rvin[1]
#     zi   = rvin[2]
#     ri   = (xi**2 + yi**2 + zi**2)**0.5
#     Ri   = (xi**2 + yi**2)**0.5
#     A    = Ri**2 + (Rd + (zi**2 + zh**2)**0.5)**2
#     B    = Ri + (zi**2 + zh**2)**0.5
#     B    = Rd + (zi**2 + zh**2)**0.5
#     acc  = -G4*Md/A**(3./2)*np.array([xi,yi,B*zi/(zi**2 + zh**2)**0.5])
#     return acc #kpc/Myr^2

def accbulgeplum_MW_ac(rvin,Mviri):
    Mpl   = GLfrac_Mbulge*Mviri
    rpl   = 0.3      # BO2015 (0-25%CB)
    rin  = (np.dot(rvin,rvin))**0.5
    acc   = -G4*Mpl/(rin**2+rpl**2)**1.5*rvin
    return acc #kpc/Myr^2


####################################################################
########## NFW FOR M31 (ANDROMEDA) #################################
def acchalonfw_M31(rvin):
    cin  = 8.74603379
    rvir = 266.3    #kpc
    rs   = rvir/cin
    Mvir = 1.1E+12 #Ms 0.8-1.1 x10^10Msun Tamm 2012, 
    Om   = 0.3
    A    = -G4*Mvir/(np.log(1.+cin)-cin/(1.+cin))
    rin  = (np.dot(rvin,rvin))**0.5
    acc  = -A*(rvin/rin*(1./(rin*(rin+rs)))-np.log(1.+rin/rs)/rin**2*rvin/rin)
    if (rin>=rvir): acc =-G4*Mvir/(rin*rin)*(rvin/rin) #(kpc/Myr)^2/kpc
    return acc #kpc/Myr^2


##############################################################################
########### POTENTIALS FOR  MALIN 1 ##########################################
##############################################################################
def acc_Malin1_nfw(rvin,Mvir,rvir,cin):
    rs   = rvir/cin
    Mdmg = Mvir  #Ms    
    A    = -G4*(Mdmg)/(np.log(1.+cin)-cin/(1.+cin))
    rin  = (np.dot(rvin,rvin))**0.5
    acc  = -A*(rvin/rin*(1./(rin*(rin+rs)))-np.log(1.+rin/rs)/rin**2*rvin/rin)
    if (rin>=rvir): acc =-G4*Mdmg/(rin*rin)*(rvin/rin) #(kpc/Myr)^2/kpc
    return acc #kpc/Myr^2

def acc_Malin1_disk(rvin,Mviri,incl=0.,PA=90.):
    zh   = 1. #kpc #  SAHA?
    Rd   = 10. #kpc # ?? SAHA??
    Md   = 1.0E+11 #Ms # from ??
    # Md   = GLfrac_Mdisk*Mviri #Ms BO15 = stellar disk+ gas disk
    PAang   = PA-90.
    inclang = incl
    
    xi   = rvin[0]
    yi   = rvin[1]
    zi   = rvin[2]
    if (inclang!=0.) or (PAang!=0.):
        xir,yir,zir    = Rotz_v2(xi,yi,zi,PAang)
        xirr,yirr,zirr = Rotx_v2(xir,yir,zir,inclang)
        xi,yi,zi       = xirr,yirr,zirr
        ri             = (xi**2 + yi**2 + zi**2)**0.5
        Ri             = (xi**2 + yi**2)**0.5
        A              = Ri**2 + (Rd + (zi**2 + zh**2)**0.5)**2
        B              = Rd + (zi**2 + zh**2)**0.5
        acc            = -G4*Md/A**(3./2)*np.array([xi,yi,B*zi/(zi**2 + zh**2)**0.5])
        accx_r,accy_r,accz_r    = Rotx_v2(acc[0],acc[1],acc[2],-inclang)
        accx_rr,accy_rr,accz_rr = Rotz_v2(accx_r,accy_r,accz_r,-PAang)
        acc = np.array([accx_rr,accy_rr,accz_rr])
    else: # avoids rotations in case PA,incl = 0 to make it faster
        ri             = (xi**2 + yi**2 + zi**2)**0.5
        Ri             = (xi**2 + yi**2)**0.5
        A              = Ri**2 + (Rd + (zi**2 + zh**2)**0.5)**2
        B              = Rd + (zi**2 + zh**2)**0.5
        acc            = -G4*Md/A**(3./2)*np.array([xi,yi,B*zi/(zi**2 + zh**2)**0.5])        
    
    return acc #kpc/Myr^2

def accbulgeplum_Malin1(rvin):
    Mpl   = 1.E10
    rpl   = 1.      
    rin  = (np.dot(rvin,rvin))**0.5
    acc   = -G4*Mpl/(rin**2+rpl**2)**1.5*rvin
    return acc #kpc/Myr^2


def acc_Malin1(rvin,Mviri,rviri,ci,i=0.,PA=90.):
    acc = accbulgeplum_Malin1(rvin)+acc_Malin1_disk(rvin,Mviri,i,PA) + acc_Malin1_nfw(rvin,Mviri,rviri,ci)
    return acc

# def acc_Malin1(rvin,Mviri,rviri,ci,i=0):
#     acc = acc_Malin1_disk(rvin,Mviri,i)
#     return acc

# def acc_Malin1(rvin,Mviri,rviri,ci,i=0):
#     acc = acc_Malin1_nfw(rvin,Mviri,rviri,ci)
#     return acc


########################################################################
########## POTENTIAL GRADIENTS FOR THE FORNAX CLUSTER ##################
def acc_Fornax_nfw_Mast(rvin):
    rin  = (np.dot(rvin,rvin))**0.5
    rs   = 120.
#     rho_in= 1.157707E+09
#     MnfwR = 4.*np.pi*rho_in*(rs_in**3.)*(np.log(1.0+r_in/rs_in)-(r_in/(rs_in+r_in)))
#     rvir = 266.3    #kpc
#     rs   = rvir/cin
    rvir = 978.0 # Mastropietro+2021
    Mvir = 1.0E+14 # Mastropietro+2021  
    cin  = 8.15 # Mastropietro+2021    
    Om   = 0.31
    A    = -G4*Mvir/(np.log(1.+cin)-cin/(1.+cin))
    acc  = -A*(rvin/rin*(1./(rin*(rin+rs)))-np.log(1.+rin/rs)/rin**2*rvin/rin)
    if (rin>=rvir): acc =-G4*Mvir/(rin*rin)*(rvin/rin) #(kpc/Myr)^2/kpc
    #acc  = -G4*Mvir/rin**3*rvin
    return acc #kpc/Myr^2

def acc_Fornax_nfw(rvin):
    rin  = (np.dot(rvin,rvin))**0.5
    rs   = 120.
#     rho_in= 1.157707E+09
#     MnfwR = 4.*np.pi*rho_in*(rs_in**3.)*(np.log(1.0+r_in/rs_in)-(r_in/(rs_in+r_in)))
#     rvir = 266.3    #kpc
#     rs   = rvir/cin
    rvir = 978.0 
    Mvir = 7.0E+13   
    cin  = 5.78330723     
    Om   = 0.31
    A    = -G4*Mvir/(np.log(1.+cin)-cin/(1.+cin))
    acc  = -A*(rvin/rin*(1./(rin*(rin+rs)))-np.log(1.+rin/rs)/rin**2*rvin/rin)
    if (rin>=rvir): acc =-G4*Mvir/(rin*rin)*(rvin/rin) #(kpc/Myr)^2/kpc
    #acc  = -G4*Mvir/rin**3*rvin
    return acc #kpc/Myr^2

def acc_Fornax_nfw_ac(rvin,Mvir,rvir,cin):
    rs   = rvir/cin
    Mdmg = Mvir  #Ms    
    A    = -G4*(Mdmg)/(np.log(1.+cin)-cin/(1.+cin))
    rin  = (np.dot(rvin,rvin))**0.5
    acc  = -A*(rvin/rin*(1./(rin*(rin+rs)))-np.log(1.+rin/rs)/rin**2*rvin/rin)
    if (rin>=rvir): acc =-G4*Mdmg/(rin*rin)*(rvin/rin) #(kpc/Myr)^2/kpc
    return acc #kpc/Myr^2


  
########################################################################################################
## LOADING PRE-COMPUTED REDSHIFT DEPENDENT PARAMETERS Mvir,Rvir, z , c from Cosmological parameters ####
########################################################################################################
fdatcos       = path_inputs+"Data_Cosmo/"
## Loading cosmological model from astropy Planck2015
z_red = 0
# print("den_DE("+str(z_red)+")=",Planck15.Ode(z_red)*Planck15.critical_density(z_red).to(u.solMass/u.kpc**3))
# print("den_M("+str(z_red)+")=",(Planck15.Om(z_red))*Planck15.critical_density(z_red).to(u.solMass/u.kpc**3))
# print("OM_DM("+str(z_red)+")=",Planck15.Om(z_red))
###################################################################################  
## Loading PARAMETERS Mvir,Rvir, z , c from Cosmological c-M Correa et al. COMMAH program

# fdatcos       = path_inputs+"/Data_Cosmo/"
# hfr           = h5py.File(fdatcos+'datacosmo13Gyr.h5', 'r')
# redshift_out  = np.array(hfr.get('redshift'))
# time_out      = np.array(hfr.get('time'))
# a_out         = Planck15.scale_factor(redshift_out)
# H_out         = (Planck15.H(redshift_out).to(1./u.myr)).value
# Hkmskpc_out   = Planck15.H(redshift_out).value
# hfr.close()

hfr           = h5py.File(fdatcos+'datacosmo13Gyr.h5', 'r')
redshift_out  = np.array(hfr.get('redshift'))
time_out      = np.array(hfr.get('time'))
Mvir_out      = np.array(hfr.get('Mz'))
Rvir_out      = np.array(hfr.get('Rvir'))
c_out         = np.array(hfr.get('c'))
redshift_out  = np.array(hfr.get('redshift'))
a_out         = Planck15.scale_factor(redshift_out)
H_out         = (Planck15.H(redshift_out).to(1./u.myr)).value
Hkmskpc_out   = Planck15.H(redshift_out).value
hfr.close()

deltavir_out = (18.*np.pi**2+82.*(Planck15.Om(redshift_out)-1)-39*(Planck15.Om(redshift_out)-1)**2)/(Planck15.Om(redshift_out))
denDMbar_out = (Planck15.Om(redshift_out))*Planck15.critical_density(redshift_out).to(u.solMass/u.kpc**3).value

########################################################################
## only available for dt=0.05,0.1,0.5,1 Myr: IN DELOREAN_COSMO_v1.ipynb you can compute tables for new time steps
hfr           = h5py.File(fdatcos+'data_DK_dt1.0Myr_v1.h5', 'r')
K_dt1Myr_out  = np.array(hfr.get('intK'))
D_dt1Myr_out  = np.array(hfr.get('intD'))
time_dt1Myr_out  = np.array(hfr.get('t_DK'))
hfr.close()

hfr                = h5py.File(fdatcos+'data_DK_dt0.5Myr_v1.h5', 'r')
K_dt0pt5Myr_out    = np.array(hfr.get('intK'))
D_dt0pt5Myr_out    = np.array(hfr.get('intD'))
time_dt0pt5Myr_out = np.array(hfr.get('t_DK'))
hfr.close()

hfr                = h5py.File(fdatcos+'data_DK_dt0.1Myr_v1.h5', 'r')
K_dt0pt1Myr_out    = np.array(hfr.get('intK'))
D_dt0pt1Myr_out    = np.array(hfr.get('intD'))
time_dt0pt1Myr_out = np.array(hfr.get('t_DK'))
hfr.close()

hfr                = h5py.File(fdatcos+'data_DK_dt0.05Myr_v1.h5', 'r')
K_dt0pt05Myr_out   = np.array(hfr.get('intK'))
D_dt0pt05Myr_out   = np.array(hfr.get('intD'))
time_dt0pt05Myr_out= np.array(hfr.get('t_DK'))
hfr.close()

f_K_dt0pt1Myr      = interp1d(time_dt0pt1Myr_out, K_dt0pt1Myr_out, kind='cubic')
f_D_dt0pt1Myr      = interp1d(time_dt0pt1Myr_out, D_dt0pt1Myr_out, kind='cubic')

f_K_dt0pt05Myr     = interp1d(time_dt0pt05Myr_out, K_dt0pt05Myr_out, kind='cubic')
f_D_dt0pt05Myr     = interp1d(time_dt0pt05Myr_out, D_dt0pt05Myr_out, kind='cubic')

f_K_dt1Myr         = interp1d(time_dt1Myr_out, K_dt1Myr_out, kind='cubic')
f_D_dt1Myr         = interp1d(time_dt1Myr_out, D_dt1Myr_out, kind='cubic')

f_K_dt0pt5Myr      = interp1d(time_dt0pt5Myr_out, K_dt0pt5Myr_out, kind='cubic')
f_D_dt0pt5Myr      = interp1d(time_dt0pt5Myr_out, D_dt0pt5Myr_out, kind='cubic')

hfr                = h5py.File(fdatcos+'data_Kp2_dt1.0Myr_v1.h5', 'r')
Kp2_dt1Myr_out     = np.array(hfr.get('intKp2'))
hfr.close()
hfr                = h5py.File(fdatcos+'data_Kp2_dt0.5Myr_v1.h5', 'r')
Kp2_dt0pt5Myr_out  = np.array(hfr.get('intKp2'))
hfr.close()
hfr                = h5py.File(fdatcos+'data_Kp2_dt0.1Myr_v1.h5', 'r')
Kp2_dt0pt1Myr_out  = np.array(hfr.get('intKp2'))
hfr.close()
hfr                = h5py.File(fdatcos+'data_Kp2_dt0.05Myr_v1.h5', 'r')
Kp2_dt0pt05Myr_out = np.array(hfr.get('intKp2'))
hfr.close()

f_Kp2_dt0pt1Myr      = interp1d(time_dt0pt1Myr_out, Kp2_dt0pt1Myr_out, kind='cubic')
f_Kp2_dt0pt05Myr     = interp1d(time_dt0pt05Myr_out, Kp2_dt0pt05Myr_out, kind='cubic')
f_Kp2_dt1Myr         = interp1d(time_dt1Myr_out, Kp2_dt1Myr_out, kind='cubic')
f_Kp2_dt0pt5Myr      = interp1d(time_dt0pt5Myr_out, Kp2_dt0pt5Myr_out, kind='cubic')

f_redshift_t = interp1d(time_out, redshift_out, kind='cubic')
f_a_t        = interp1d(time_out, a_out, kind='cubic')
f_H_t        = interp1d(time_out, H_out, kind='cubic')
f_Hkmskpc_t  = interp1d(time_out, Hkmskpc_out, kind='cubic')
f_denDMbar_t = interp1d(time_out, denDMbar_out, kind='cubic')

f_deltavir_t = interp1d(time_out, deltavir_out, kind='cubic')


f_Mvir_t       = interp1d(time_out, Mvir_out, kind='cubic')
f_c_t          = interp1d(time_out, c_out, kind='cubic')
f_Rvir_t       = interp1d(time_out, Rvir_out, kind='cubic')

f_redshift_t   = interp1d(time_out, redshift_out, kind='cubic')
f_a_t          = interp1d(time_out, a_out, kind='cubic')
f_H_t          = interp1d(time_out, H_out, kind='cubic')

hfr2a          = h5py.File(fdatcos+'datacosmo13Gyr10E9Ms.h5', 'r')
redshift_out2a = np.array(hfr2a.get('redshift'))
time_out2a     = np.array(hfr2a.get('time'))
Mvir_out2a     = np.array(hfr2a.get('Mz'))
Rvir_out2a     = np.array(hfr2a.get('Rvir'))
c_out2a        = np.array(hfr2a.get('c'))
hfr2a.close()
f_Mvir_t2a     = interp1d(time_out2a, Mvir_out2a, kind='cubic')
f_c_t2a        = interp1d(time_out2a, c_out2a, kind='cubic')
f_Rvir_t2a     = interp1d(time_out2a, Rvir_out2a, kind='cubic')

hfr2b          = h5py.File(fdatcos+'datacosmo13Gyr16e9Ms.h5', 'r')
redshift_out2b = np.array(hfr2b.get('redshift'))
time_out2b     = np.array(hfr2b.get('time'))
Mvir_out2b     = np.array(hfr2b.get('Mz'))
Rvir_out2b     = np.array(hfr2b.get('Rvir'))
c_out2b        = np.array(hfr2b.get('c'))
hfr2b.close()
f_Mvir_t2b     = interp1d(time_out2b, Mvir_out2b, kind='cubic')
f_c_t2b        = interp1d(time_out2b, c_out2b, kind='cubic')
f_Rvir_t2b     = interp1d(time_out2b, Rvir_out2b, kind='cubic')

# COSMO DATA FOR FORNAX CLUSTER
hfrF2          = h5py.File(fdatcos+'cosmohalo13Gyr_Mvir7E13Ms.h5', 'r')
redshift_outF2 = np.array(hfrF2.get('redshift'))
time_outF2     = np.array(hfrF2.get('time'))
Mvir_outF2     = np.array(hfrF2.get('Mz'))
Rvir_outF2     = np.array(hfrF2.get('Rvir'))
c_outF2        = np.array(hfrF2.get('c'))
hfrF2.close()
f_Mvir_tF2     = interp1d(time_outF2, Mvir_outF2, kind='cubic')
f_c_tF2        = interp1d(time_outF2, c_outF2, kind='cubic')
f_Rvir_tF2     = interp1d(time_outF2, Rvir_outF2, kind='cubic')


##############################################################################################
## PRE COMPUTING ORBITS FOR M31 AND MW TO BE USED LATER FOR THE LG SATELLITES ################
##############################################################################################
def leapM31(objclass,case,ti,dti,ri,vi,riPOT,viPOT,Mviri,Rviri,ci,*argv):
    global integ,integmethod,vartimestep,nsteporb ## Parameters to control orbit integration method and accuracy
    if case=="case1_cosmo" or case=="case2_cosmo":
        if integmethod=='SYM': # argv=a_i, a_i+1, f2_i+1/2, f1_i, f2_i+1
            ai   = argv[0]
            ai1  = argv[1]
            Da05= argv[2]
            Kai = argv[3]
            Da1 = argv[4]
            pi  = ai**2*vi
            ri  = ri + pi*Da05 # Drift: kpc
            acc = accase(objclass,case,ti,dti,ri,vi,riPOT,viPOT,0.,Mviri,Rviri,ci) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2
            pi  = pi + acc*Kai   # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr  # + or - ??  
            vi  = pi/ai1**2
            ri  = ri + pi*Da1    # Drift: kpc        
        else:
            ai,Hi = argv[0],argv[1]
            ri  = ri+ (0.5*dti*vi)  # Drift: kpc
            acc = accase(objclass,case,ti,dti,ri,vi,riPOT,viPOT,0.,Mviri,Rviri,ci)
            vi  = vi*(1-Hi*dti)/(1+Hi*dti) + dti*acc/(1.+Hi*dti)/ai**3    # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr  # + or - ??  
            ri  = ri + (0.5*dti*vi) # Drift: kpc
    else:
        ri  = ri + (0.5*dti*vi) # Drift: kpc
        acc = accase(objclass,case,ti,dti,ri,vi,riPOT,viPOT,0.,Mviri,Rviri,ci) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2
        vi  = vi + dti*acc      # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr
        ri  = ri + (0.5*dti*vi) # Drift: kpc
    return ri,vi

def leapMW(objclass,case,ti,dti,ri,vi,riPOT,*argv):
    global integ,integmethod,vartimestep,nsteporb ## Parameters to control orbit integration method and accuracy
    if case=="case1_cosmo" or case=="case2_cosmo":
        if integmethod=='SYM':# argv=a_i, a_i+1, f2_i+1/2, f1_i, f2_i+1
            ai   = argv[0]
            ai1  = argv[1]
            Da05= argv[2]
            Kai = argv[3]
            Da1 = argv[4]
#             Kaip2 = argv[5]
            pi  = ai**2*vi
            ri  = ri + pi*Da05 # Drift: kpc
            acc = acchalonfw_M31(ri-riPOT) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2
            pi  = pi + acc*Kai   # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr  # + or - ??  
            vi  = pi/ai1**2
            ri  = ri + pi*Da1    # Drift: kpc
        else:
            ai,Hi = argv[0],argv[1]
            ri  = ri+ (0.5*dti*vi)  # Drift: kpc
            acc = acchalonfw_M31(ri-riPOT)
            vi  = vi*(1-Hi*dti)/(1+Hi*dti) + dti*acc/(1.+Hi*dti)/ai**3    # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr  # + or - ??  
            ri  = ri + (0.5*dti*vi) # Drift: kpc            
    else:    
        ri  = ri + (0.5*dti*vi) # Drift: kpc
        acc = acchalonfw_M31(ri-riPOT) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2
        vi  = vi + dti*acc      # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr
        ri  = ri + (0.5*dti*vi) # Drift: kpc        
    return ri,vi

def orbMWM31(case,tin,tf,dti,posoMW,veloMW,posoM31,veloM31):
    global integ,integmethod,vartimestep,nsteporb ## Parameters to control orbit integration method and accuracy
    integ,vartimestep, nsteporb,integmethod = 'DKD',False,100.,'SYM' ## Parameters to control orbit integration method and accuracy
    
    if (tf<0.): dti=-dti
    ti = np.arange(tin,tf+dti,dti)

    if integ =="DKD": # it neads acceleration at t+0.5dt
        Mviri,Rviri,ci = MRcfun(case,ti+dti*0.5)
    elif integ =="KDK":
        Mviri,Rviri,ci = MRcfun(case,ti)    
    
    if case=="case1_cosmo" or case=="case2_cosmo":
        if integmethod=='SYM':
            par1 = f_a_t(ti) # = a_i
            par2 = f_a_t(ti+dti) # =a_i+1
            par3,par4,par5,par6 = GetCosParDKD2(dti,ti,ti,ti+dti*0.5,ti)  # =D_i+1/2(dt/2), K_i(dt), D_i+1(dt/2)  Problem: is not instant value, but from integral and val=val(0)+dt, no initial val(t=0)=0
#             par3,par4,par5 = GetCosParFa(dti,ti,ti,ti+dti/2)  # =D_i+1/2(dt/2), K_i(dt), D_i+1(dt/2)  Problem: is not instant value, but from integral and val=val(0)+dt, no initial val(t=0)=0
        else:
            par1,par2 = GetCosPar(ti+dti*0.5) # par1,par2 = a_i+1/2,H_i+1/2

    GLfrac_Mdisk  = GLMdisk/Mviri[0]
    GLfrac_Mbulge = GLMbulge/Mviri[0]
    
    veloMW = veloMW*kms2kpcMyr
    veloM31 = veloM31*kms2kpcMyr
    t_it = tin
    pos_listMW,vel_listMW = np.zeros((len(ti),len(posoMW))),np.zeros((len(ti),len(veloMW)))
    pos_listMW[0,:],vel_listMW[0,:] = posoMW,veloMW
    pos_listM31,vel_listM31 = np.zeros((len(ti),len(posoM31))),np.zeros((len(ti),len(veloM31)))
    pos_listM31[0,:],vel_listM31[0,:] = posoM31,veloM31
    pos_itM31,vel_itM31     = posoM31,veloM31
    start = timer.time()
    for i in range(0,len(ti)):
        t_it                    = t_it+dti #Myr 
        if case=="case1_cosmo" or case=="case2_cosmo":
            if integmethod=='SYM':  
                pos_itMW,vel_itMW       = leapMW(objclass,case,t_it,dti,posoMW,veloMW,pos_itM31,par1[i],par2[i],par3[i],par4[i],par5[i],par6[i])
                pos_itM31,vel_itM31     = leapM31(objclass,case,t_it,dti,posoM31,veloM31,pos_itMW,vel_itMW,Mviri[i],Rviri[i],ci[i],par1[i],par2[i],par3[i],par4[i],par5[i],par6[i])
            else:
                pos_itMW,vel_itMW       = leapMW(objclass,case,t_it,dti,posoMW,veloMW,pos_itM31,par1[i],par2[i])
                pos_itM31,vel_itM31     = leapM31(objclass,case,t_it,dti,posoM31,veloM31,pos_itMW,vel_itMW,Mviri[i],Rviri[i],ci[i],par1[i],par2[i])
        else:
            pos_itMW,vel_itMW       = leapMW(case,t_it,dti,posoMW,veloMW,pos_itM31)
            pos_itM31,vel_itM31     = leapM31(case,t_it,dti,posoM31,veloM31,pos_itMW,vel_itMW,Mviri[i],Rviri[i],ci[i])
        
        pos_listMW[i,:],vel_listMW[i,:] = pos_itMW,vel_itMW
        posoMW,veloMW           = pos_itMW,vel_itMW
        
        pos_listM31[i,:],vel_listM31[i,:] = pos_itM31,vel_itM31
        posoM31,veloM31         = pos_itM31,vel_itM31
    print(timer.time() - start),"s"    
    return pos_listMW,vel_listMW*kpcMyr2kms,pos_listM31,vel_listM31*kpcMyr2kms,ti


# case1: MWconstant, case2: MW(time), case3: MWconstant+M31, case4: MW(time)+M31, case5: MW(time)+M31(orb2)
def orbMWM31case(case):
    if case=="case1" or case=="case1a" or case=="case1b" or case=="case1nbody"or case=="caseCOtest" or \
        case=="case1_cosmo" or case=="case2_cosmo" or case=="case2" or case=="case2a" or case=="case2b" \
        or case=="case5"or case=="case6" or case=="case1_kepler" or case=="Fornax_Setup1" or case=="case5" \
        or case=="case6" or case=="case1_kepler" or case=="Fornax_Setup2":
        posoMWf,veloMWf,posoM31f=orb0_posMW_f,orb0_velMW_f,orb1_posM31_f 
    elif case=="case3_cosmo" or case=="case4_cosmo":
        posoMWf,veloMWf,posoM31f=orb1_cosmo_posMW_f,orb1_cosmo_velMW_f,orb1_cosmo_posM31_f
    elif case=="case3" or case=="case4":
        posoMWf,veloMWf,posoM31f=orb1_posMW_f,orb1_velMW_f,orb1_posM31_f
    elif case=="case3a":
        posoMWf,veloMWf,posoM31f=orb1a_posMW_f,orb1a_velMW_f,orb1a_posM31_f
    elif case=="case3b":
        posoMWf,veloMWf,posoM31f=orb1b_posMW_f,orb1b_velMW_f,orb1b_posM31_f        
    elif case=="case4a":
        posoMWf,veloMWf,posoM31f=orb1a_posMW_f,orb1a_velMW_f,orb1a_posM31_f        
    elif case=="case4b":
        posoMWf,veloMWf,posoM31f=orb1b_posMW_f,orb1b_velMW_f,orb1b_posM31_f
    elif case=="case5fut" or case=="case1fut" or case=="case1_burk" or case=="case2_burk" or case=="case1_nfw":
        posoMWf,veloMWf,posoM31f=orb_fut0_posMW_f,orb_fut0_velMW_f,orb_fut0_posMW_f        
    else:
        posoMWf,veloMWf,posoM31f=orb_fut0_posMW_f,orb_fut0_velMW_f,orb_fut0_posMW_f  
        # posoMWf,veloMWf,posoM31f=orb0_posMW_f,orb0_velMW_f,orb1_posM31_f
    return posoMWf,veloMWf,posoM31f

useM31MWorbs = True
# useM31MWorbs = False

# calculateMWM31orbits = True
calculateMWM31orbits = False
if calculateMWM31orbits and useM31MWorbs:
    print("Calculating MW-M31 orbits: True")
    MWx  = 0.  # kpc
    MWy  = 0.  # kpc
    MWz  = 0.  # kpc
    MWpos= np.array([MWx,MWy,MWz])
    MWvx = 0.  # km/s
    MWvy = 0.  # km/s
    MWvz = 0.  # km/s
    MWvel= np.array([MWvx,MWvy,MWvz])
    # # reading data as class
    fnameIC  = path_inputs+"/Data_IC/M31_IC.pkl"
    
    with open(fnameIC, 'rb') as input:
        M31a = pickle.load(input)
        M31b = pickle.load(input) 
    
    M31pos    = np.array([M31a.GC.x.value,M31a.GC.y.value,M31a.GC.z.value])
    M31vel    = np.array([M31a.GC.v_x.value,M31a.GC.v_y.value,M31a.GC.v_z.value])
    M31pos_b  = np.array([M31b.GC.x.value,M31b.GC.y.value,M31b.GC.z.value])
    M31vel_b  = np.array([M31b.GC.v_x.value,M31b.GC.v_y.value,M31b.GC.v_z.value])

    t_init,t_fin,dt=0., -12200., 0.1 # Myr, Myr, Myr
    case = "case1_cosmo"
    orb1_cosmo_posMW,orb1_cosmo_velMW,orb1_cosmo_posM31,orb1_cosmo_velM31,orb1_cosmo_timeMWM31 = orbMWM31(case,t_init,t_fin,dt,MWpos,MWvel,M31pos,M31vel)    
    
    t_init,t_fin,dt=0., -12200., 0.1 # Myr, Myr, Myr
    case = "case1"
    orb1_posMW,orb1_velMW,orb1_posM31,orb1_velM31,orb1_timeMWM31 = orbMWM31(case,t_init,t_fin,dt,MWpos,MWvel,M31pos,M31vel)
    orb0_posMW,orb0_velMW = np.zeros((len(orb1_posMW),3)),np.zeros((len(orb1_velMW),3))    
    case = "case1a"
    orb1a_posMW,orb1a_velMW,orb1a_posM31,orb1a_velM31,orb1a_timeMWM31 = orbMWM31(case,t_init,t_fin,dt,MWpos,MWvel,M31pos,M31vel)
    case = "case1b"
    orb1b_posMW,orb1b_velMW,orb1b_posM31,orb1b_velM31,orb1b_timeMWM31 = orbMWM31(case,t_init,t_fin,dt,MWpos,MWvel,M31pos,M31vel)
    case = "case1"
    orb2_posMW,orb2_velMW,orb2_posM31,orb2_velM31,orb2_timeMWM31 = orbMWM31(case,t_init,t_fin,dt,MWpos,MWvel,M31pos_b,M31vel_b)
    t_init, t_fin, dt = 0., 12200., 0.1 # Myr, Myr, Myr
    case = "case1"
    orb1_fut_posMW,orb1_fut_velMW,orb1_fut_posM31,orb1_fut_velM31,orb1_fut_timeMWM31 = orbMWM31(case,t_init,t_fin,dt,MWpos,MWvel,M31pos,M31vel)
    case = "case1"
    orb2_fut_posMW,orb2_fut_velMW,orb2_fut_posM31,orb2_fut_velM31,orb2_fut_timeMWM31 = orbMWM31(case,t_init,t_fin,dt,MWpos,MWvel,M31pos_b,M31vel_b)
        
    fnameMWM31= path_inputs+"Data_Orbits_MWM31/data_orbitsMW_M31_v6.h5"
    os.system("rm "+fnameMWM31)
    hf = h5py.File(fnameMWM31, 'w')
    hf.create_dataset('orb0_posMW',data=orb0_posMW)
    hf.create_dataset('orb0_velMW',data=orb0_velMW)
    
    hf.create_dataset('orb1_posMW',data=orb1_posMW)
    hf.create_dataset('orb1_velMW',data=orb1_velMW)
    hf.create_dataset('orb1_posM31',data=orb1_posM31)
    hf.create_dataset('orb1_velM31',data=orb1_velM31)
    
    hf.create_dataset('orb1a_posMW',data=orb1a_posMW)
    hf.create_dataset('orb1a_velMW',data=orb1a_velMW)
    hf.create_dataset('orb1a_posM31',data=orb1a_posM31)
    hf.create_dataset('orb1a_velM31',data=orb1a_velM31)
    
    hf.create_dataset('orb1b_posMW',data=orb1b_posMW)
    hf.create_dataset('orb1b_velMW',data=orb1b_velMW)
    hf.create_dataset('orb1b_posM31',data=orb1b_posM31)
    hf.create_dataset('orb1b_velM31',data=orb1b_velM31)
    
    hf.create_dataset('orb2_posMW',data=orb2_posMW)
    hf.create_dataset('orb2_velMW',data=orb2_velMW)
    hf.create_dataset('orb2_posM31',data=orb2_posM31)
    hf.create_dataset('orb2_velM31',data=orb2_velM31)
    
    hf.create_dataset('orb1_fut_posMW',data=orb1_fut_posMW)
    hf.create_dataset('orb1_fut_velMW',data=orb1_fut_velMW)
    hf.create_dataset('orb1_fut_posM31',data=orb1_fut_posM31)
    hf.create_dataset('orb1_fut_velM31',data=orb1_fut_velM31)
    hf.create_dataset('orb1_fut_timeMWM31',data=orb1_fut_timeMWM31)
    
    hf.create_dataset('orb2_fut_posMW',data=orb2_fut_posMW)
    hf.create_dataset('orb2_fut_velMW',data=orb2_fut_velMW)
    hf.create_dataset('orb2_fut_posM31',data=orb2_fut_posM31)
    hf.create_dataset('orb2_fut_velM31',data=orb2_fut_velM31)
    hf.create_dataset('orb2_fut_timeMWM31',data=orb2_fut_timeMWM31)
    
    hf.create_dataset('orb1_timeMWM31',data=orb1_timeMWM31)
    hf.create_dataset('orb2_timeMWM31',data=orb2_timeMWM31)

    hf.create_dataset('orb1_cosmo_posMW',data=orb1_cosmo_posMW)
    hf.create_dataset('orb1_cosmo_velMW',data=orb1_cosmo_velMW)
    hf.create_dataset('orb1_cosmo_posM31',data=orb1_cosmo_posM31)
    hf.create_dataset('orb1_cosmo_velM31',data=orb1_cosmo_velM31)    
    hf.create_dataset('orb1_cosmo_timeMWM31',data=orb1_cosmo_timeMWM31)
    hf.close()

    
if not calculateMWM31orbits:
    # fnameMWM31  = path_inputs+'Data_Orbits_MWM31/data_orbitsMW_M31_v6.h5'
    # fnameMWM31  = path_inputs+'Data_Orbits_MWM31/data_orbitsMW_M31_v7.h5'
    fnameMWM31    = path_inputs+'Data_Orbits_MWM31/data_orbitsMW_M31_v8.h5'    
    hfr           = h5py.File(fnameMWM31, 'r')
    orb1_timeMWM31= np.array(hfr.get('orb1_timeMWM31'))
    orb2_timeMWM31= np.array(hfr.get('orb2_timeMWM31'))

    orb0_posMW    = np.array(hfr.get('orb0_posMW'))
    orb0_velMW    = np.array(hfr.get('orb0_velMW'))
    
    orb1_posMW    = np.array(hfr.get('orb1_posMW'))
    orb1_velMW    = np.array(hfr.get('orb1_velMW'))
    orb1_posM31   = np.array(hfr.get('orb1_posM31'))
    orb1_velM31   = np.array(hfr.get('orb1_velM31'))

    orb1a_posMW    = np.array(hfr.get('orb1a_posMW'))
    orb1a_velMW    = np.array(hfr.get('orb1a_velMW'))
    orb1a_posM31   = np.array(hfr.get('orb1a_posM31'))
    orb1a_velM31   = np.array(hfr.get('orb1a_velM31'))

    orb1b_posMW    = np.array(hfr.get('orb1b_posMW'))
    orb1b_velMW    = np.array(hfr.get('orb1b_velMW'))
    orb1b_posM31   = np.array(hfr.get('orb1b_posM31'))
    orb1b_velM31   = np.array(hfr.get('orb1b_velM31'))

    orb2_posMW    = np.array(hfr.get('orb2_posMW'))
    orb2_velMW    = np.array(hfr.get('orb2_velMW'))
    orb2_posM31   = np.array(hfr.get('orb2_posM31'))
    orb2_velM31   = np.array(hfr.get('orb2_velM31'))

    orb1_fut_posMW    = np.array(hfr.get('orb1_fut_posMW'))
    orb1_fut_velMW    = np.array(hfr.get('orb1_fut_velMW'))
    orb1_fut_posM31   = np.array(hfr.get('orb1_fut_posM31'))
    orb1_fut_velM31   = np.array(hfr.get('orb1_fut_velM31'))
    orb1_fut_timeMWM31= np.array(hfr.get('orb1_fut_timeMWM31'))

    orb2_fut_posMW    = np.array(hfr.get('orb2_fut_posMW'))
    orb2_fut_velMW    = np.array(hfr.get('orb2_fut_velMW'))
    orb2_fut_posM31   = np.array(hfr.get('orb2_fut_posM31'))
    orb2_fut_velM31   = np.array(hfr.get('orb2_fut_velM31'))
    orb2_fut_timeMWM31= np.array(hfr.get('orb2_fut_timeMWM31'))

    orb1_cosmo_posMW    = np.array(hfr.get('orb1_cosmo_posMW'))
    orb1_cosmo_velMW    = np.array(hfr.get('orb1_cosmo_velMW'))
    orb1_cosmo_posM31   = np.array(hfr.get('orb1_cosmo_posM31'))
    orb1_cosmo_velM31   = np.array(hfr.get('orb1_cosmo_velM31'))

    orb1_cosmo_timeMWM31 = np.array(hfr.get('orb1_cosmo_timeMWM31'))

    hfr.close()

# plotMWM31orbits = True
plotMWM31orbits = False
if plotMWM31orbits and useM31MWorbs:   
    stp=10
    orb1  = (orb1_posMW[::stp,:]-orb1_posM31[::stp,:])
    orb2  = (orb2_posMW[::stp,:]-orb2_posM31[::stp,:])
    t     = orb1_timeMWM31[::stp]
    plotorbMWM31(orb1,orb2,t)

    orb1_fut  = (orb1_fut_posMW[::stp,:]-orb1_fut_posM31[::stp,:])
    orb2_fut  = (orb2_fut_posMW[::stp,:]-orb2_fut_posM31[::stp,:])
    t         = orb2_fut_timeMWM31[::stp]
    plotorbMWM31(orb1_fut,orb2_fut,t)

    fig, ax = plt.subplots(1, 1,figsize=(12,5))
    stp=10
    ax.plot(orb1_posMW[::stp,0],orb1_posMW[::stp,1],"-b")
    ax.plot(orb1_posM31[::stp,0],orb1_posM31[::stp,1],"--b")
    ax.plot(orb2_posMW[::stp,0],orb2_posMW[::stp,1],"-g")
    ax.plot(orb2_posM31[::stp,0],orb2_posM31[::stp,1],"--g")
    plt.show()
    # fig, ax = plt.subplots(1, 1,figsize=(12,5))
    stp=10
    # ax.plot(orb1_fut_posMW[::stp,0],orb1_fut_posMW[::stp,1],"-b")
    # ax.plot(orb1_fut_posM31[::stp,0],orb1_fut_posM31[::stp,1],"--b")
    # ax.plot(orb2_fut_posMW[::stp,0],orb2_fut_posMW[::stp,1],"-g")
    # ax.plot(orb2_fut_posM31[::stp,0],orb2_fut_posM31[::stp,1],"--g")
    # plt.show()
    # fig, ax = plt.subplots(1, 1,figsize=(12,5))
    # stp=10
    # ax.plot(orb1_posMW[::stp,0],orb1_posMW[::stp,1],"-b")
    # ax.plot(orb1_cosmo_posM31[::stp,0],orb1_cosmo_posM31[::stp,1],"--b")
    # plt.show()
    orb1  = (orb1_posMW[::stp,:]-orb1_posM31[::stp,:])
    orb1_cosmo  = (orb1_cosmo_posMW[::stp,:]-orb1_cosmo_posM31[::stp,:])*np.array([f_a_t(orb1_timeMWM31[::stp]),f_a_t(orb1_timeMWM31[::stp]),f_a_t(orb1_timeMWM31[::stp])]).T
    t     = orb1_timeMWM31[::stp]
    plotorbMWM31(orb1,orb1_cosmo,t)
    
    
if useM31MWorbs: 
        # Setup Orbits MW and M31 for Satellite orbits
        rPot_time    = np.arange(-15000.,15000+1000,1000)
        rPot         = np.zeros((len(rPot_time),3))
        orb_fut0_posMW_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb_fut0_posMW_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb_fut0_posMW_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb_fut0_posMW_f(x):
            return np.array([orb_fut0_posMW_x(x),orb_fut0_posMW_y(x),orb_fut0_posMW_z(x)]).T

        vPot_time    = rPot_time.copy() #np.arange(0.,13000+1000,1000)
        vPot         = np.zeros((len(rPot_time),3))
        orb_fut0_velMW_x = interp1d(vPot_time, vPot[:,0], kind='cubic')
        orb_fut0_velMW_y = interp1d(vPot_time, vPot[:,1], kind='cubic')
        orb_fut0_velMW_z = interp1d(vPot_time, vPot[:,2], kind='cubic')
        def orb_fut0_velMW_f(x):
            return np.array([orb_fut0_velMW_x(x),orb_fut0_velMW_y(x),orb_fut0_velMW_z(x)]).T

        rPot_time    = orb1_timeMWM31
        rPot         = orb0_posMW
        orb0_posMW_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb0_posMW_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb0_posMW_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb0_posMW_f(x):
            return np.array([orb0_posMW_x(x),orb0_posMW_y(x),orb0_posMW_z(x)]).T

        vPOT         = orb0_velMW
        orb0_velMW_x = interp1d(rPot_time, vPOT[:,0], kind='cubic')
        orb0_velMW_y = interp1d(rPot_time, vPOT[:,1], kind='cubic')
        orb0_velMW_z = interp1d(rPot_time, vPOT[:,2], kind='cubic')
        def orb0_velMW_f(x):
            return np.array([orb0_velMW_x(x),orb0_velMW_y(x),orb0_velMW_z(x)]).T

        rPot         = orb1_posMW
        orb1_posMW_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1_posMW_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1_posMW_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1_posMW_f(x):
            return np.array([orb1_posMW_x(x),orb1_posMW_y(x),orb1_posMW_z(x)]).T

        vPOT         = orb1_velMW
        orb1_velMW_x = interp1d(rPot_time, vPOT[:,0], kind='cubic')
        orb1_velMW_y = interp1d(rPot_time, vPOT[:,1], kind='cubic')
        orb1_velMW_z = interp1d(rPot_time, vPOT[:,2], kind='cubic')
        def orb1_velMW_f(x):
            return np.array([orb1_velMW_x(x),orb1_velMW_y(x),orb1_velMW_z(x)]).T

        rPot         = orb1a_posMW
        orb1a_posMW_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1a_posMW_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1a_posMW_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1a_posMW_f(x):
            return np.array([orb1a_posMW_x(x),orb1a_posMW_y(x),orb1a_posMW_z(x)]).T

        vPOT         = orb1a_velMW
        orb1a_velMW_x = interp1d(rPot_time, vPOT[:,0], kind='cubic')
        orb1a_velMW_y = interp1d(rPot_time, vPOT[:,1], kind='cubic')
        orb1a_velMW_z = interp1d(rPot_time, vPOT[:,2], kind='cubic')
        def orb1a_velMW_f(x):
            return np.array([orb1a_velMW_x(x),orb1a_velMW_y(x),orb1a_velMW_z(x)]).T

        rPot         = orb1b_posMW
        orb1b_posMW_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1b_posMW_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1b_posMW_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1b_posMW_f(x):
            return np.array([orb1b_posMW_x(x),orb1b_posMW_y(x),orb1b_posMW_z(x)]).T
        
        vPOT         = orb1b_velMW
        orb1b_velMW_x = interp1d(rPot_time, vPOT[:,0], kind='cubic')
        orb1b_velMW_y = interp1d(rPot_time, vPOT[:,1], kind='cubic')
        orb1b_velMW_z = interp1d(rPot_time, vPOT[:,2], kind='cubic')
        def orb1b_velMW_f(x):
            return np.array([orb1b_velMW_x(x),orb1b_velMW_y(x),orb1b_velMW_z(x)]).T

        rPot         = orb2_posMW
        rPot_time    = orb2_timeMWM31
        orb2_posMW_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb2_posMW_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb2_posMW_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb2_posMW_f(x):
            return np.array([orb2_posMW_x(x),orb2_posMW_y(x),orb2_posMW_z(x)]).T

        vPOT         = orb2_velMW
        rPot_time    = orb2_timeMWM31
        orb2_velMW_x = interp1d(rPot_time, vPOT[:,0], kind='cubic')
        orb2_velMW_y = interp1d(rPot_time, vPOT[:,1], kind='cubic')
        orb2_velMW_z = interp1d(rPot_time, vPOT[:,2], kind='cubic')
        def orb2_velMW_f(x):
            return np.array([orb2_velMW_x(x),orb2_velMW_y(x),orb2_velMW_z(x)]).T

        rPot             = orb1_fut_posMW
        rPot_time        = orb1_fut_timeMWM31
        orb1_fut_posMW_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1_fut_posMW_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1_fut_posMW_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1_fut_posMW_f(x):
            return np.array([orb1_fut_posMW_x(x),orb1_fut_posMW_y(x),orb1_fut_posMW_z(x)]).T

        vPOT         = orb1_fut_velMW
        orb1_fut_velMW_x = interp1d(rPot_time, vPOT[:,0], kind='cubic')
        orb1_fut_velMW_y = interp1d(rPot_time, vPOT[:,1], kind='cubic')
        orb1_fut_velMW_z = interp1d(rPot_time, vPOT[:,2], kind='cubic')
        def orb1_fut_velMW_f(x):
            return np.array([orb1_fut_velMW_x(x),orb1_fut_velMW_y(x),orb1_fut_velMW_z(x)]).T

        rPot          = orb1_posM31
        rPot_time     = orb1_timeMWM31
        orb1_posM31_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1_posM31_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1_posM31_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1_posM31_f(x):
            return np.array([orb1_posM31_x(x),orb1_posM31_y(x),orb1_posM31_z(x)]).T

        rPot           = orb1a_posM31
        rPot_time      = orb1_timeMWM31
        orb1a_posM31_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1a_posM31_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1a_posM31_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1a_posM31_f(x):
            return np.array([orb1a_posM31_x(x),orb1a_posM31_y(x),orb1a_posM31_z(x)]).T

        rPot           = orb1b_posM31
        rPot_time      = orb1_timeMWM31
        orb1b_posM31_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1b_posM31_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1b_posM31_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1b_posM31_f(x):
            return np.array([orb1b_posM31_x(x),orb1b_posM31_y(x),orb1b_posM31_z(x)]).T

        rPot          = orb2_posM31
        rPot_time     = orb2_timeMWM31
        orb2_posM31_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb2_posM31_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb2_posM31_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb2_posM31_f(x):
            return np.array([orb2_posM31_x(x),orb2_posM31_y(x),orb2_posM31_z(x)]).T

        rPot          = orb1_fut_posM31
        rPot_time     = orb1_fut_timeMWM31
        orb1_fut_posM31_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1_fut_posM31_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1_fut_posM31_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1_fut_posM31_f(x):
            return np.array([orb1_fut_posM31_x(x),orb1_fut_posM31_y(x),orb1_fut_posM31_z(x)]).T


        rPot           = orb1_cosmo_posMW
        rPot_time      = orb1_cosmo_timeMWM31
        orb1_cosmo_posMW_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1_cosmo_posMW_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1_cosmo_posMW_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1_cosmo_posMW_f(x):
            return np.array([orb1_cosmo_posMW_x(x),orb1_cosmo_posMW_y(x),orb1_cosmo_posMW_z(x)]).T

        vPOT         = orb1_cosmo_velMW
        orb1_cosmo_velMW_x = interp1d(rPot_time, vPOT[:,0], kind='cubic')
        orb1_cosmo_velMW_y = interp1d(rPot_time, vPOT[:,1], kind='cubic')
        orb1_cosmo_velMW_z = interp1d(rPot_time, vPOT[:,2], kind='cubic')
        def orb1_cosmo_velMW_f(x):
            return np.array([orb1_cosmo_velMW_x(x),orb1_cosmo_velMW_y(x),orb1_cosmo_velMW_z(x)]).T

        rPot          = orb1_cosmo_posM31
        rPot_time     = orb1_cosmo_timeMWM31
        orb1_cosmo_posM31_x = interp1d(rPot_time, rPot[:,0], kind='cubic')
        orb1_cosmo_posM31_y = interp1d(rPot_time, rPot[:,1], kind='cubic')
        orb1_cosmo_posM31_z = interp1d(rPot_time, rPot[:,2], kind='cubic')
        def orb1_cosmo_posM31_f(x):
            return np.array([orb1_cosmo_posM31_x(x),orb1_cosmo_posM31_y(x),orb1_cosmo_posM31_z(x)]).T