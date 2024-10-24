#########################################################################################################################################
## ######################################################################################################################################
## CODE: ORBIT CALCULATOR DELOREAN (v2) (update 18.01.2024)
## AUTHOR: MATIAS A. BLANA D.
## CHILE, SANTIAGO JANUARY (18.01.2024) 2024
## VERSION : MW DWARF SATELLITES INCLUDING MW MULTIPLE POTENTIALS and M31, AND Fornax Cluster potentials, WITH ORBITS WITH COSMIC EXPANSION OPTIONS
## SCRIPT  : MAIN SCRIPT THAT LOADS INTEGRATIONS ROUTINES, AND WHERET TO DEFINE OBSERVABLES 
## #######################################################################################################################################
##########################################################################################################################################
from DELOREAN_pylib import *
from DELOREAN_UC_v1 import *    ## IMPORT CONSTANTS AND UNITS 
from DELOREAN_IC_v1 import *    ## IMPORT INITIAL CONDITIONS GENERATOR FOR DELOREAN INTEGRATOR
from DELOREAN_FUN_v1 import * 
from DELOREAN_SETUP_v1 import * 

#################################################################
# MAIN ORBIT INTEGRATIONS FUNCTIONS ##################################
#################################################################
## Cases with cosmo INCLUDE Cosmic expansion with a symplectic integrator
## in Comoving Coordinates from Quinn et al. 1997 http://arxiv.org/abs/astro-ph/9710043 (used in Gadget2 Springer 2005)
## which re-writes the bottom equations with the specific canonical linear momentum and potentials derived from Lagrangian
# new technique: symplectic
# ai1,ai2,Hi = argv[0],argv[1],argv[2]
# pi  = ai1**2*vi
# ri  = ri+ (0.5*dti*pi)*(1/ai2**2-1./ai1**2)  # Drift: kpc
# acc = accase(case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci,ai1,Hi)
# vi  = pi - dti*acc/(1.+Hi*dti)    # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr  # + or - ??  
# ri  = ri + (0.5*dti*vi) # Drift: kpc
########################################################################################################################
# Version used for cosmological eq. of motion versions with KDK, used in older Gadget 1 version (shown also in Quin et al. 1997) 
# r(t)  = a(t) x(t)
# dr/dt = x(t)da/dt + a(t)dx/dt
# v(t)  = x(t) H(z) + a(t) dx/dt
# a     = 1/(1+z)
# H(z)  = (da/dt)/a
#########################################################################################################################
def leapAdapM_v2(objclass,param_it,*argv):
    case,ti,dti,ri,vi,riPOT1,riPOT1f,viPOT1,viPOT1f,riPOT2,riPOT2f,Mviri,Rviri,ci = param_it
    if integ=="DKD": ## using Drift-Kick-Drift implementation
        if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
            dtau0,acc,rout = dtcrit(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci,*argv)
        else:
            dtau0,acc,rout = dtcrit(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci)
        Ndt   = 1
        dtiab = abs(dti)
        if not vartimestep: dtau0=1E100            
        if (dtau0>=dtiab):
            if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
                if integmethod=='SYM': # Symplectic method: argv=a_i, a_i+1, f2_i+1/2, f1_i, f2_i+1
                    ai  = argv[0]
                    ai1 = argv[1]
                    Da05= argv[2]
                    Kai = argv[3]
                    Da1 = argv[4]
                    Kp2 = argv[5]
                    ri  = rout            # Drift: kpc
                    pi  = ai**2*vi
                    pi  = pi + acc*Kai + Kp2*ri  # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr  # + or - ??  
                    vi  = pi/ai1**2
                    ri  = ri + pi*Da1    # Drift: kpc
                else: # Technique old: not symplectic # argv=a_i+1/2,H_i+1/2
                    ai  = argv[0] 
                    Hi  = argv[1] 
                    ri  = rout              # Drift: kpc
                    vi  = vi*(1-Hi*dti)/(1+Hi*dti) + dti*acc/(1.+Hi*dti)/ai**3    # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr  # + or - ?
                    ri  = ri + (0.5*dti*vi) # Drift: kpc
            else: # all  cases not cosmological
                ri  = rout              # Drift: kpc
                vi  = vi + dti*acc      # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr    
                ri  = ri + (0.5*dti*vi) # Drift: kpc
        elif(dtau0<dtiab):
            dtii   = dti
            dtiiab = abs(dtii)
            while (dtau0<dtiiab):
                Ndt    = Ndt +1
                dtii   = 0.5*dtii
                dtiiab = abs(dtii)
                tii    = ti + dtii
                riPOT1 = riPOT1f(tii)
                viPOT1 = viPOT1f(tii)*kms2kpcMyr
                riPOT2 = riPOT2f(tii)
                Mviri,Rviri,ci = objclass.setup.MRcfun_v2(objclass,np.array([tii]))
                if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
                    dtau0,acc,rout = dtcrit(objclass,case,tii,dtii,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci,*argv)
                else:
                    dtau0,acc,rout = dtcrit(objclass,case,tii,dtii,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci)
            dtii   = dti/Ndt
#             dtii   = dtau0*np.sign(dti)
#             Ndt    = nsteporb
            tii    = ti
            tiiv   = np.arange(tii,tii+dtii*Ndt,dtii) + dtii*0.5
            riPOT1 = riPOT1f(tiiv)
            viPOT1 = viPOT1f(tiiv)*kms2kpcMyr
            riPOT2 = riPOT2f(tiiv)
            Mviri,Rviri,ci = objclass.setup.MRcfun_v2(objclass,tiiv)
            if case =="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
                aii,Hii=GetCosPar(tiiv)
            for i in range(0,Ndt):
                tii   =  tii + dtii             
                if case =="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
                    ri  = ri+ (0.5*dtii*vi)   # Drift: kpc
                    acc = objclass.setup.accase(objclass,case,tii,dtii,ri,vi,riPOT1[i],viPOT1[i],riPOT2[i],Mviri[i],Rviri[i],ci[i])
                    vi  = vi*(1-Hii[i]*dtii)/(1+Hii[i]*dtii) + dtii*acc/(1.+Hii[i]*dtii)/aii[i]**3    # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr  # + or - ??  
                    ri  = ri + (0.5*dtii*vi) # Drift: kpc                    
                else:
                    ri  = ri + (0.5*dtii*vi) # Drift: kpc            
                    acc = objclass.setup.accase(objclass,case,tii,dtii,ri,vi,riPOT1[i],viPOT1[i],riPOT2[i],Mviri[i],Rviri[i],ci[i]) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2
                    vi  = vi + dtii*acc      # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr
                    ri  = ri + (0.5*dtii*vi) # Drift: kpc
        return ri,vi
    elif integ=="KDK": ## using Kick-Drift-Kick implementation
        if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
            ai,Hi = argv[0],argv[0]         
            dtau0,acc,vout = dtcrit(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci,ai,Hi)
        else:
            dtau0,acc,vout = dtcrit(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci)
        
        if not vartimestep: dtau0=1E100
        Ndt   = 1
        dtiab = abs(dti)
        if (dtau0>=dtiab):
            if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
            # Technique old: not symplectic
                vi  = vout               # Kick: kpc/Myr
                ri  = ri + dti*vi        # Drift: kpc
                riPOT1 = riPOT1f(ti+dti)
                viPOT1 = viPOT1f(ti+dti)*kms2kpcMyr
                riPOT2 = riPOT2f(ti+dti)
                Mviri,Rviri,ci,ai,Hi = objclass.setup.MRcfun_v2(objclass,np.array([ti+dti]))
                acc = objclass.setup.accase(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci)
                vi  = vi*(1-Hi*0.5*dti)/(1+Hi*0.5*dti) + 0.5*dti*acc/(1.+Hi*0.5*dti)/ai**3 # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr 
            else:
                vi  = vout              # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr    
                ri  = ri + dti*vi       # Drift: kpc
                riPOT1 = riPOT1f(ti+dti)
                viPOT1 = viPOT1f(ti+dti)*kms2kpcMyr
                riPOT2 = riPOT2f(ti+dti)
                Mviri,Rviri,ci = objclass.setup.MRcfun_v2(objclass,np.array([ti+dti]))
                acc = objclass.setup.accase(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci)
                vi  = vi + 0.5*dti*acc # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr                 
        elif(dtau0<dtiab):
            dtii   = dti
            dtiiab = abs(dtii)
            while (dtau0<dtiiab):            
                Ndt    = Ndt +1
                dtii   = 0.5*dtii
                dtiiab = abs(dtii)
                tii    = ti + dtii
                Mviri,Rviri,ci = objclass.setup.MRcfun_v2(objclass,np.array([tii]))
                riPOT1 = riPOT1f(tii)
                viPOT1 = viPOT1f(tii)*kms2kpcMyr
                riPOT2 = riPOT2f(tii)            
                dtau0,acc,vout = dtcrit(objclass,case,tii,dtii,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci)
            dtii = dti/Ndt
            tii  = ti
            tiiv   = np.arange(tii,tii+dtii*Ndt+dtii,dtii)
            riPOT1 = riPOT1f(tiiv)
            viPOT1 = viPOT1f(tiiv)*kms2kpcMyr
            riPOT2 = riPOT2f(tiiv)
            Mviri,Rviri,ci = objclass.setup.MRcfun_v2(objclass,tiiv)
            if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
                aii,Hii = GetCosPar(tiiv)
            for i in range(0,Ndt):
                tii    = ti + dtii
                if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
                    acc   = objclass.setup.accase(objclass,case,tii,dtii,ri,vi,riPOT1[i],viPOT1[i],riPOT2[i],Mviri[i],Rviri[i],ci[i])
                    vi    = vi*(1-Hii[i]*0.5*dtii)/(1+Hii[i]*0.5*dtii) + 0.5*dtii*acc/(1.+Hii[i]*0.5*dtii)/aii[i]**3        # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr
                    ri    = ri+ dtii*vi                                                                               # Drift: kpc
                    acc   = objclass.setup.accase(objclass,case,tii,dtii,ri,vi,riPOT1[i+1],viPOT1[i+1],riPOT2[i+1],Mviri[i+1],Rviri[i+1],ci[i+1])                    
                    vi    = vi*(1-Hii[i+1]*0.5*dtii)/(1+Hii[i+1]*0.5*dtii) + 0.5*dtii*acc/(1.+Hii[i+1]*0.5*dtii)/aii[i+1]**3 # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr
                else:
                    acc   = objclass.setup.accase(objclass,case,tii,dtii,ri,vi,riPOT1[i],viPOT1[i],riPOT2[i],Mviri[i],Rviri[i],ci[i]) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2                    
                    vi    = vi + 0.5*dtii*acc      # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr
                    ri    = ri + dtii*vi       # Drift: kpc            
                    acc   = objclass.setup.accase(objclass,case,tii,dtii,ri,vi,riPOT1[i+1],viPOT1[i+1],riPOT2[i+1],Mviri[i+1],Rviri[i+1],ci[i+1])
                    vi    = vi + 0.5*dtii*acc      # Kick: Myr*[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2 = kpc/Myr                
        return ri,vi


#########################################################################
########### FUNCTION GET CRITERIA TO INCRESE dt RESOLUTION ##############
def dtcrit(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci,*argv):
#     global integ
    if integ=='DKD': ## using Drift-Kick-Drift implementation
        if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
            if integmethod=='SYM': # using symplectic method: argv=a_i, a_i+1, f2_i+1/2, f1_i, f2_i+1
                ai  = argv[0]
                Da05= argv[2]
                pi  = ai**2*vi
                ri  = ri + pi*Da05 # Drift: kpc
                acc = objclass.setup.accase(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2
                dtau0 = 1E100 # no variable time step yet
            else: # not-using symplectic method: 
                ai  = argv[0]
                ri  = ri + (0.5*dti*vi) # Drift: kpc
                acc = objclass.setup.accase(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2
                di  = (np.dot(ri-riPOT1,ri-riPOT1))**0.5 + 1.E-10
                ari = (np.dot(acc,acc))**0.5/ai**2 + 1.E-10
                dtau0 = 1E100 # no variable time step yet
        else: # all  cases not cosmological
            ri  = ri + (0.5*dti*vi) # Drift: kpc
            acc = objclass.setup.accase(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2
            di  = (np.dot(ri-riPOT1,ri-riPOT1))**0.5 + 1.E-10
            ari = (np.dot(acc,acc))**0.5 + 1.E-10
            vci2= (ari*di)
            vi2 = np.dot(vi-viPOT1,vi-viPOT1)
            vmi = (vci2 + vi2)**0.5
            dtau0 = 2.*np.pi*di/vmi/nsteporb # criteria for very radial orbits
        return dtau0,acc,ri
    elif integ=='KDK': ## using Kick-Drift-Kick implementation
        if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
            ai,Hi = argv[0],argv[1]
            acc = objclass.setup.accase(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci) #kpc/Myr^2
            vi  = vi*(1-Hi*0.5*dti)/(1+Hi*0.5*dti) + 0.5*dti*acc/(1.+Hi*0.5*dti)/ai**3 # Kick: kpc/Myr
            di  = (np.dot(ri-riPOT1,ri-riPOT1))**0.5 + 1.E-10
            ari = (np.dot(acc,acc))**0.5/ai**2 + 1.E-10
            dtau0 = 2.*np.pi*np.sqrt(abs(di/ari))/nsteporb  # criteria for circular orbits         
        else:
            acc = objclass.setup.accase(objclass,case,ti,dti,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci) #[(km/s)^2/kpc]*((kpc/Myr)/(km/s))^2=kpc/Myr^2
            vi  = vi + 0.5*dti*acc # Kick: kpc/Myr
            di  = (np.dot(ri-riPOT1,ri-riPOT1))**0.5 + 1.E-10
            ari = (np.dot(acc,acc))**0.5 + 1.E-10
            dtau0 = 2.*np.pi*np.sqrt(abs(di/ari))/nsteporb  # criteria for circular orbits
        return dtau0,acc,vi  

def GetCosPar(t_in):
    Hi        = f_H_t(t_in)
    ai        = f_a_t(t_in)
    return ai,Hi

def GetCosPar2(t_in):
    f1        = f_intf1a(t_in)
    f2        = f_intf2a(t_in)
    f105      = f_intf1a05(t_in)
    f205      = f_intf2a05(t_in)
    return ai,Hi

def GetCosParFa(dt,t_in1,t_in2,t_in3): # =D_i+1/2, K_i, D_i+1
    if abs(dt)==0.1:
        Ddt05  = f_D_dt0pt05Myr(t_in1)
        Kdt1   = f_K_dt0pt1Myr(t_in2)
        Ddt1   = f_D_dt0pt05Myr(t_in3)
    if abs(dt)==1.0:
        Ddt05  = f_D_dt0pt5Myr(t_in1)
        Kdt1   = f_K_dt1Myr(t_in2)
        Ddt1   = f_D_dt0pt5Myr(t_in3)        
    return Ddt05,Kdt1,Ddt1

def GetCosParDKD2(dt,t_in1,t_in2,t_in3,t_in4): # =D_i+1/2, K_i, D_i+1,Kp2_i
    if abs(dt)==0.1:
        Ddt05  = f_D_dt0pt05Myr(t_in1)
        Kdt1   = f_K_dt0pt1Myr(t_in2)
        Ddt1   = f_D_dt0pt05Myr(t_in3)
        Kp2dt1 = f_Kp2_dt0pt1Myr(t_in4)         
    if abs(dt)==1.0:
        Ddt05  = f_D_dt0pt5Myr(t_in1)
        Kdt1   = f_K_dt1Myr(t_in2)
        Ddt1   = f_D_dt0pt5Myr(t_in3)
        Kp2dt1 = f_Kp2_dt1Myr(t_in4)        
    return Ddt05,Kdt1,Ddt1,Kp2dt1
    
###############################################################################
#### FUNCTION THAT INTEGRATE ORBITAL TIME STEPS BY CALLING LEAPFROG ###########
###############################################################################
def IntegrateOrb_v2(objclass):
    case,tv,dti = objclass.setup.case,objclass.paramtime[3],objclass.paramtime[2]
    objcoords   = objclass.objcoords
    posMWif,velMWif,posM31if = objclass.paramfixpot
    posMWi      = posMWif(tv)
    posM31i     = posM31if(tv)
    velMWi      = velMWif(tv)*kms2kpcMyr
    
    if integ =="DKD": # it neads acceleration at t+0.5dt
        Mviri,Rviri,ci = objclass.setup.MRcfun_v2(objclass,tv+dti*0.5)
    elif integ =="KDK":
        Mviri,Rviri,ci = objclass.setup.MRcfun_v2(objclass,tv)
   
    if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
        if integ =="DKD": # it neads acceleration at t+0.5dt
            if integmethod=='SYM':
                par1 = f_a_t(tv) # = a_i
                par2 = f_a_t(tv+dti) # =a_i+1
                par3,par4,par5,par6 = GetCosParDKD2(dti,tv,tv,tv+dti*0.5,tv)  # =D_i+1/2(dt/2), K_i(dt), D_i+1(dt/2)  Problem: is not instant value, but from integral and val=val(0)+dt, no initial val(t=0)=0
            else:
                par1,par2 = GetCosPar(tv+dti*0.5) # par1,par2 = a_i+1/2,H_i+1/2
        elif integ =="KDK":
            return "KDK not available yet for cosmological orbit"
            par1,par2 = GetCosPar(tv+dti*0.5)

    global GLfrac_Mdisk 
    global GLfrac_Mbulge
    GLfrac_Mdisk  = GLMdisk/Mviri[0]
    GLfrac_Mbulge = GLMbulge/Mviri[0]

    posi,veli = objcoords
    veli = veli*kms2kpcMyr
    pos_list,vel_list = np.zeros((len(tv),len(posi))),np.zeros((len(tv),len(veli)))
    pos_list[0,:],vel_list[0,:] = posi+veli*(-0.5)*dti,veli
    pos_it,vel_it = posi + veli*(-0.5)*dti, veli

    for i in range(1,len(tv)): # WHERE ORBIT INTEGRATION IS PERFORMED
        if case=="case1_cosmo" or case=="case2_cosmo" or case=="case3_cosmo" or case=="case4_cosmo":
            if integmethod=='SYM':
                param_it            = case,tv[i],dti,pos_it,vel_it,posMWi[i],posMWif,velMWi[i],velMWif,posM31i[i],posM31if,Mviri[i],Rviri[i],ci[i]
                pos_it,vel_it       = leapAdapM_v2(objclass,param_it,par1[i],par2[i],par3[i],par4[i],par5[i],par6[i])
            else:
                param_it            = case,tv[i],dti,pos_it,vel_it,posMWi[i],posMWif,velMWi[i],velMWif,posM31i[i],posM31if,Mviri[i],Rviri[i],ci[i]
                pos_it,vel_it       = leapAdapM_v2(objclass,param_it,par1[i],par2[i])
        else:
            param_it                = case,tv[i],dti,pos_it,vel_it,posMWi[i],posMWif,velMWi[i],velMWif,posM31i[i],posM31if,Mviri[i],Rviri[i],ci[i]
            pos_it,vel_it           = leapAdapM_v2(objclass,param_it)
        pos_list[i,:],vel_list[i,:] = pos_it,vel_it    
        if case=="case1nbody":
            if (i % 10)==0: print("t=",tv[i],"Myr")
    
    orb    = np.hstack([pos_list, vel_list*kpcMyr2kms])
    try: 
        if objclass.orbprop:
            orbprop = OrbProp_v1(objclass,orb)
            return orb,orbprop
    except:
        return [orb]
    #         return [orb,orbprop]
    # except:
    #     return [orb]

    
######################################################################################################
#### FUNCTIONS THAT READ INPUTS TO BE USED TO BY INTEGRATION FUNCTIONS AND GENERATE OUTPUTS ##########
def FunInputOutput_v2(paramvar,objclass): 
    objclass.objcoords       = CoordIN2CoordORB_v2(objclass,paramvar) #objcoords
    orb                      = IntegrateOrb_v2(objclass)
    return orb

def FunInputOutput_chunk_v2(paramvars,objclass): # For when running in parallel with mp.Pool
    chunk_size  = len(paramvars)
    fun_chunk1  = []
    for i in range(chunk_size):
        orb = FunInputOutput_v2(paramvars[i],objclass)
        fun_chunk1.append(orb)

    return fun_chunk1
    
# ######################################################################################################
# #### FUNCTIONS THAT READ INPUTS TO BE USED TO BY INTEGRATION FUNCTIONS AND GENERATE OUTPUTS ##########
# def FunInputOutput_v2(paramvar,objclass): 
#     objclass.objcoords       = CoordIN2CoordORB_v2(objclass,paramvar) #objcoords
#     orb                      = IntegrateOrb_v2(objclass)
#     if objclass.orbprop: 
#         orbprop = OrbProp_v1(objclass,orb)    
#         return orb,orbprop
#     else:
#         return orb

# def FunInputOutput_chunk_v2(paramvars,objclass): # For when running in parallel with mp.Pool
#     chunk_size  = len(paramvars)
#     fun_chunk1  = []
#     if objclass.orbprop : fun_chunk2  = []
#     for i in range(chunk_size):
#         if objclass.orbprop:
#             orb,orbprop = FunInputOutput_v2(paramvars[i],objclass)
#             fun_chunk1.append(orb)
#             fun_chunk2.append(orbprop)
    
#     if objclass.orbprop: 
#         return fun_chunk1,fun_chunk2
#     else:
#         return fun_chunk1

######################################################################################################
### FUNCTIONS TO MEASURE ORBITAL PROPERTIES: Rperi, Rapo, etc.
def OrbProp_v1(objclass,orb):
    # pos      = orb[:,] # kpc, km/s
    # vel      = orb[:,3]
    time     = objclass.paramtime[3]
    x,y,z    = orb[:,0],orb[:,1],orb[:,2]
    vx,vy,vz = orb[:,3],orb[:,4],orb[:,5]
    R        = (x**2 + y**2)**0.5
    r        = (R**2 +z**2)**0.5
    phi      = np.arctan2(y,x)*180./np.pi
    phi[phi<0]= phi[phi<0] + 360.
    phi0     = 0
    phiswitch= True
    # determine minimum pericentre (Rp_min) and maximum Apocentre (Ra_max)
    Rpmin    = np.amin(r)
    TRpmin   = time[np.argmin(r)]
    Ramax    = np.amax(r)
    TRamax   = time[np.argmax(r)]
    amajor   = (Rpmin + Ramax)*0.5
    bminor   = (Rpmin*Ramax)**0.5
    el       = 1.-bminor/amajor

    Ra       = [] # array with apocentres
    TRa      = [] # array with time of apocentres
    Rp       = [] # array with pericentres = dmin
    TRp      = [] # array with time of pericentres
    Tphi     = [] # array with 2pi angles

    getorbperiod=False
    getorbperiod=False
    if getorbperiod:
        for i in range(0,len(time)-2):
            dr1 = r[i+1]-r[i]
            dr2 = r[i+2]-r[i+1]
            if (dr1>=0)and(dr2<=0): ## from increasing r to decreasing (apocentres), the (=) also works for circles (makes switch=False)
                Ra.append(r[i+1])
                TRa.append(time[i+1])
                if phiswitch: 
                    phi0 = phi[i+1]
                    timephi0 = time[i+1]
                    # Tphi.append(time[i+1])
                    # print('phi switch1')
                phiswitch = False
            elif (dr1<=0)and(dr2>=0): ## from decreasing r to increasing (Pericentres), the (=) also works for circles (makes switch=False)
                Rp.append(r[i+1])
                TRp.append(time[i+1])
                if phiswitch: 
                    phi0 = phi[i+1]
                    timephi0 = time[i+1]
                    # Tphi.append(time[i+1])
                    # print('phi switch2')
                phiswitch = False

            if not phiswitch:
                if time[i+1]!=timephi0:
                    rphi  = (phi - phi0) % 360
                    drphi1 = rphi[i+1] - rphi[i]
                    drphi2 = rphi[i+2] - rphi[i+1]
                    if drphi2<0 and drphi1>=0: # anti-clockwise orbits 
                        # print('anti-clockwise')
                        Tphi.append(time[i+1])
                    # elif drphi2>0 and drphi1<0: # clockwise orbits 
                    #     Tphi.append(time[i+1])
                    #     print('clockwise')
    orbprop = [[Rpmin,TRpmin,Ramax,TRamax,amajor,bminor,el], [Ra,TRa,Rp,TRp,Tphi] ]
    return orbprop


#################################################################################################################################
#################################################################################################################################
#### MAIN PROGRAM: CALLS ALL FUNCTIONS AND RUNS THE ORBIT INTEGRATION AND DOES THE PARAMETER EXPLORATION DEFINED BY THE USER ####
#################################################################################################################################
def ExploreParam_v2(listobjs): 
    ## VERSION FOR EMCEE #####################
    ## THINGS TO DO IN THE FUTURE: (if gravity between satellites is allowed: THEN VARY PARAMS BETWEEN DIFFERENT OBJECTS MEANS NPARS_SAT1 X NPARS_SAT2 X ... X NPARS_SATN  
    global integ,integmethod,vartimestep,nsteporb ## Parameters to control orbit integration method and accuracy
    integ,vartimestep, nsteporb,integmethod = 'DKD',True,100.,'SYM' ## Parameters to control orbit integration method and accuracy
    print("integration parameters: integ,vartimestep,nsteporb,integmethod=",integ,vartimestep,nsteporb,integmethod)

    Nobj = len(listobjs)
    listobjs_out = []
        
    for i in range(0,Nobj):
        objclass  = listobjs[i]
        orbs      = []
        chi2v     = []
        try: # CHECK IF emcee ATTRIBUTE EXISTS, IN CASE IS USING A SAVED OLD OBJECT FORMAT
            objclass.emcee.use 
        except:
            objclass.emcee = Emcee()
            objclass.emcee.use=False
        
        if (not objclass.emcee.use): # THIS BLOCK RUNS DELOREAN FOR A PRE-DEFINED SET OF PARAMETERS (eg. grid of parameters, etc)
            try:
                ncpu      = objclass.ncpu
            except:
                print('ncpu not found, assuming ncpu=1')
                ncpu = 1
            # case      = objclass.case
            npars     = len(objclass.paramvars)         
            print(f"Selected ncpu: {ncpu} of {os.cpu_count()}")
            if ncpu>npars: 
                ncpu = npars       
                print(f"ncpu cannot be more than Npars: reset to {ncpu} cpu")
            print(f"npars={npars}")
            ti,tf,dt = objclass.paramtime
            if (tf<ti) and (dt>0): dt=-dt
            t    = np.arange(ti,tf+dt,dt)
            objclass.paramtime[2] = dt
            objclass.paramtime.append(t)
            paramvars = objclass.paramvars
            objclass.paramfixpot  = orbMWM31case(objclass.setup.case)

            if ncpu==1: # OPTION TO RUN ORBIT CALCULATIONS IN 1 CPU, WITHOUT CALLING mp.Pool
                start_timer  = timer.time()
                for j in range(0,npars):
                    if (j % max(int(0.1*npars),1))or(j==0)==0: print(f"orbit n={j+1}/{npars} ({int(100*(j+1)/npars)}%)")
                    paramvar = objclass.paramvars[j]
                    # if objclass.orbprop:
                    #     orb      = FunInputOutput_v2(paramvar,objclass)
                    # else:
                    orb      = FunInputOutput_v2(paramvar,objclass)
                    if j==0:
                        print(f" One orbit took {timer.time() - start_timer} s")
                        print(f" All will take {npars*(timer.time() - start_timer)/60.} min")
                    if (j % max(int(0.1*npars),1))==0:print(f" Running time {timer.time() - start_timer} s")
                    orbs.append(orb)
                objclass.orbs = orbs                

            elif ncpu>1: # OPTION TO RUN ORBIT CALCULATIONS IN PARALLEL USING mp.Pool from multiprocessing library
                chunk_size = len(paramvars) // ncpu
                chunks     = [paramvars[i:i+chunk_size] for i in range(0, len(paramvars), chunk_size)]
                # print(f"chunks {chunks}")
                print(f"Starting time estimation of calculations for {npars} orbits")
                nspeed       = chunks[0][0:max(min(int(float(len(chunks[0]))/10.),10),1)]
                print(f"len(nspeed) {len(nspeed)}")
                start_timer  = timer.time()
                resultspeed  = FunInputOutput_chunk_v2(nspeed,objclass)             
                elapsed_time = (timer.time()- start_timer)/len(nspeed)
                print(f"One orbit calculation took {elapsed_time} seconds to run in 1 cpu.")
                print(f"1E6 orbits would take a total of {format(elapsed_time*1E6/60.,'.1f')} min or {format(elapsed_time*1E6/60./60.,'.2f')} hrs to run in 1 cpu.")
                print(f"{npars} orbits could take a total of {format(elapsed_time*npars/60.,'.2f')} min or {format(elapsed_time*npars/60./60.,'.2f')} hrs to run  in 1 cpu.")
                print(f"{npars} orbits could take a total of {format(elapsed_time*npars/60./ncpu,'.2f')} min or {format(elapsed_time*npars/60./60./ncpu,'.2f')} hrs to run  in {ncpu} cpu.")        

                print(f"Starting parallel calculations for {npars} orbits") 
                pool       = mp.Pool(processes=ncpu)
                results    = [pool.apply_async(FunInputOutput_chunk_v2, args=(chunk,objclass)) for chunk in chunks]
                pool.close()
                pool.join()   

                orbs_chunk    = [result.get() for result in results]   #orbs_chunk, variable2 = zip(*[result.get() for result in results]) #orbs,chi2v = np.concatenate(orbs_chunk), np.concatenate(chi2_chunk)
                # orbs_chunk    = result.get() for result in results   #orbs_chunk, variable2 = zip(*[result.get() for result in results]) #orbs,chi2v = np.concatenate(orbs_chunk), np.concatenate(chi2_chunk)                
                print("len(orbs_chunk)=",len(orbs_chunk))
                print("len(orbs_chunk[0])=",len(orbs_chunk[0])) 
                print("len(orbs_chunk[0][0])=",len(orbs_chunk[0][0])) 
                # orbs          = np.concatenate(orbs_chunk)
                orbs_chunk_new  = []

                for j1 in range(0,len(orbs_chunk)):
                    for j2 in range(0,len(orbs_chunk[j1])):
                        orb      = orbs_chunk[j1][j2][0] #Testorbit1.orbs[0][0][0]
                        try:
                            if objclass.orbprop: 
                                orbprop  = orbs_chunk[j1][j2][1]
                                orbs_chunk_new.append([orb,orbprop])
                        except:
                            orbs_chunk_new.append([orb])        
                        
                        
                objclass.orbs = orbs_chunk_new   
                
        elif (objclass.emcee.use): ## THIS BLOCK RUNS EMCEE TO FIT SOME PRE-DEFINED OBSERVABLE (if objclass.emcee.use=True) ###########
            start_timer  = timer.time()
            try:
                ncpu      = objclass.emcee.ncpu
            except:
                try:
                    print('emcee.ncpu not found, using objclass.ncpu')
                    ncpu      = objclass.ncpu
                except:
                    print('ncpu not found, assuming emcee ncpu=1')
                    ncpu = 1
            npars     = len(objclass.paramvars)         
            print(f"Selected emcee ncpu: {ncpu} of {os.cpu_count()} available in the server")
            print(f"npars={npars}")
            ti,tf,dt = objclass.paramtime
            if (tf<ti) and (dt>0): dt=-dt
            t         = np.arange(ti,tf+dt,dt)
            objclass.paramtime[2] = dt            
            objclass.paramtime.append(t)
            # objclass.paramfixpot  = [orbMWM31case(objclass.case)]
            objclass.paramfixpot  = orbMWM31case(objclass.setup.case)
            paramvar              = objclass.paramvars

            # testemceefuns       = False
            testemceefuns     = True                
            if testemceefuns:
                lpost = log_posterior_v1(paramvar, objclass)
                print("testing emcee: lpost=",lpost)
                
            runemcee   = True
            # runemcee = False
            if runemcee:
                import emcee
                import corner
                ncpu          = objclass.emcee.ncpu  # Set this to the number of CPUs you want to use
                pool          = mp.Pool(processes=ncpu)
                # Set up the emcee sampler
                nwalkers      = objclass.emcee.nwalkers
                initial_guess = generate_initial_guess(nwalkers, npars, objclass.paramvarlim)
                sampler       = emcee.EnsembleSampler(nwalkers, npars, log_posterior_v1, pool=pool, args=([objclass]))             
                # Run the sampler
                n_steps       = objclass.emcee.nsteps
                sampler.run_mcmc(initial_guess, n_steps, progress=True);       
                # Get the samples
                # samples = sampler.get_chain(flat=True)
                # Get the chain of samples
                chain = sampler.chain
                plotemcee = True
                # plotemcee = False
                if plotemcee:
                    # Plot trace plots for each parameter
                    plt.figure(figsize=(8, 8))

                    for i in range(npars):
                        plt.subplot(npars, 1, i + 1)
                        for walker in range(nwalkers):
                            plt.plot(chain[walker, :, i], alpha=0.7)
                        plt.ylabel(f'Param {i + 1}')
                        plt.xlabel('Step')

                    plt.tight_layout()
                    plt.show()

                    # Get the samples
                    samples = sampler.chain[:, :, :].reshape((-1, npars))
                    # # Plot the results
                    samples = sampler.chain[:, :, :].reshape((-1, npars))

                    try:
                        labels = [f"{objclass.paramvarsunit[i]}" for i in range(npars)]                    
                    except:
                        print('labels not found')
                        labels = [f"Param {i+1}" for i in range(npars)]

                    # Plot a corner plot
                    figure = corner.corner(samples, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
                    plt.show()
                    objclass.emcee.samples = samples
                    
                best_fit_mostprob = samples[np.argmax(sampler.flatlnprobability)]
                print("Best-Fitting Parameters:", best_fit_mostprob)
                # Calculate the median of each parameter from the samples
                # parameter_medians = np.median(sampler.chain, axis=(0, 1))
                # print("Best-Fitting Parameter Medians:", parameter_medians)
                parameter_best50 = np.percentile(sampler.chain,[50], axis=(0, 1))[0]
                print("Best-Fitting Parameter 50% Percentile:",parameter_best50)
                
                calcbestorb = True
                # calcbestorb = False
                if calcbestorb:
                    bestparam_list = []
                    for i in range(0,2):
                        if i==0: bestparam = best_fit_mostprob.tolist()
                        if i==1: bestparam = parameter_best50.tolist()
                        paramvar = bestparam # it should be just 1 list with N parameters that emcee will vary           
                        orb      = FunInputOutput_v2(paramvar,objclass)
                        lpost    = log_posterior_v1(paramvar, objclass)
                        bestparam_list.append([bestparam,lpost])                        
                        orbs.append(orb)
                    objclass.bestparam =  bestparam_list
                    objclass.orbs = orbs
                    
        listobjs_out.append(objclass)
    print(f"All calculations finished after {format((timer.time()- start_timer),'.2f')} s, {format((timer.time()- start_timer)/60.,'.2f')} min, {format((timer.time()- start_timer)/60./60.,'.2f')} hrs")
  
    return listobjs_out




#############################################################################
### FOR EMCEE ###############################################################
#############################################################################
# Define the log-posterior function
# from DELOREAN_ORBINT_v1 import FunInputOutput_v2,FunInputOutput_chunk_v2

def log_posterior_v1(paramvar,objclass):
    log_prior = log_prior_v1(paramvar, objclass.paramvarlim,objclass)
    if not np.isfinite(log_prior):
        return -np.inf
    return log_prior + log_likelihood_v1(paramvar,objclass)

# Define the log-prior function with limits for arbitrary number of parameters
## Function designed for Malin-1 problem
def log_prior_v1(params,limits,objclass):
    try:
        R,vlos = objclass.paramfix2
        Mvir,Rvir,ci = MRcfun_v2(objclass,objclass,0.)
        z,vt   = params[0],params[1]
        vtot2  = vt**2 + vlos**2
        r      = (R**2 + z**2)**0.5
        ek     = 0.5*vtot2
        eu     = - G4*Mvir/r *kpcMyr2kms**2
        etot   = ek + eu
        # erat   = ek/etot
        erat   = ek/eu  # virial ratio expects K/U = -0.5, so range for bound systems should be : 0(super bound) ... -0.5 (virial) ...  -1(limit bound-unbound)      
        params.append(erat)
    except:
        pass
    for p, lim in zip(params, limits):
        if not lim[0] < p < lim[1]:
            return -np.inf
    return 0.0

def log_likelihood_v1(paramvar, objclass):
    orb              = FunInputOutput_v2(paramvar,objclass)[0] # orbs, orbsprop
    moddata          = [orb,objclass.paramtime[3]] 
    redchi2          = Chi2_v2(moddata,objclass.obsdata,objclass)
    logprob          = -0.5*(redchi2)
    return logprob

def Chi2_v2(moddata,obsdata,objclass):
    # Observations #################################
    xobs,yobs,Robserr = obsdata # Position observation and error (pixel size for example (in kpc), )
    nobsdata          = len(xobs) 
    Projdist_min_all  = np.empty((nobsdata), dtype=float)
    # Model observables ############################
    orb,time = moddata
    x       = orb[:,0] 
    y       = orb[:,1]
    for i in range(0,nobsdata): # iterations over each observation data point
        Projdist         = ((x-xobs[i])**2 + (y-yobs[i])**2)**0.5 # Distances to each observation data point
        Projdist_min     = np.min(Projdist)
        # Projdist_min_arg = np.argmin(Projdist)
        # Projdist_f       = interp1d(time,Projdist, kind='cubic') ## looks smaller distance separations for higher timesteps
        # if Projdist_min_arg>0 or Projdist_min_arg<len(Projdist)-1: # avoids orbits interpolation where orbits go away from stream.
        #     New_steps       = 100. #subdivision new timesteps 
        #     New_dt          = (time[1]-time[0])/New_steps
        #     time_inter      = np.arange(time[Projdist_min_arg-1],time[Projdist_min_arg+1],New_dt)
        #     Projdist_interp = Projdist_f(time_inter)
        #     Projdist_min    = np.min(Projdist_interp)
        Projdist_min_all[i] = Projdist_min
            
    redchi2tot = np.sum((Projdist_min_all/Robserr)**2)/nobsdata
    return redchi2tot

def generate_initial_guess(nwalkers, npars, param_limits):
    """
    Generate initial guesses for emcee walkers.

    Parameters:
    - nwalkers (int): Number of walkers.
    - npars (int): Number of parameters.
    - param_limits (list of tuples): Parameter limits. Each tuple represents the
                                     lower and upper limit for a parameter.

    Returns:
    - initial_guess (ndarray): Array of shape (nwalkers, npars) containing
                               initial guesses for all walkers.
    """
    initial_guess = np.zeros((nwalkers, npars))

    for i in range(npars):
        lower_limit, upper_limit = param_limits[i]
        initial_guess[:, i] = np.random.uniform(lower_limit, upper_limit, nwalkers) # you can change for a gaussian distribution, for example.

    return initial_guess
#############################################################################
#############################################################################