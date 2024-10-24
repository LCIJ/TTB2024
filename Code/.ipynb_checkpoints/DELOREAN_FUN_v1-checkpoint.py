#########################################################################################################################################
## ######################################################################################################################################
## CODE: ORBIT CALCULATOR DELOREAN (v2)
## AUTHOR: MATIAS A. BLANA D.
## CHILE, SANTIAGO JULY 2023
## VERSION : MW DWARF SATELLITES INCLUDING MW MULTIPLE POTENTIALS and M31, AND Fornax Cluster potentials, WITH ORBITS WITH COSMIC EXPANSION OPTIONS
## SCRIPT  : TOOL FUNCTIONS: PLOTTING, TRIGONOMETRIC, ETC
## #######################################################################################################################################
##########################################################################################################################################
from DELOREAN_pylib import *
from DELOREAN_UC_v1 import *    ## IMPORT CONSTANTS AND UNITS 
from DELOREAN_IC_v1 import *    ## IMPORT INITIAL CONDITIONS GENERATOR FOR DELOREAN INTEGRATOR
# from DELOREAN_SETUP_v1 import * 

##########################################################
## FUNCTIONS: PLOTTING ###################################
##########################################################

def multiline(xs, ys, c, ax=None, **kwargs):
    ax = plt.gca() if ax is None else ax
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)
    lc.set_array(np.asarray(c))
    ax.add_collection(lc)
    ax.autoscale()
    return lc

def plotorbMWM31(dv,dvb,time):
    fig, ax = plt.subplots(1, 1,figsize=(12,5))
    d = (dv[:,0]**2+dv[:,1]**2+dv[:,2]**2)**0.5
    db = (dvb[:,0]**2+dvb[:,1]**2+dvb[:,2]**2)**0.5
    ax.plot(time,d,"-b")
    ax.plot(time,db,"-g")
    plt.show
    return

def plotorbplaneXY(obj,L,**kwargs):
    orbits = obj.orbs[:]
    leg    = f"T={obj.paramtime[3][-1]}[Myr]"
    no = len(orbits)
#     st,va,no = orbits.shape 
#     print('st,va,no =',st,va,no)
    
    evenly_spaced_interval = np.linspace(0.5, 1, no)
    colors = [cm.Purples(x) for x in evenly_spaced_interval]
    
    fig, axs = plt.subplots(1, 1,figsize=(6,6))
    ymin,ymax=-L,L
    ls = 20    
    ax = axs

    for i in range(0,no):
        orbit = orbits[i][0] # 0: x,y,z, ..., etc. 1: orbsprop
        # if no==1: orbit = orbits[i]
        # if no>1: orbit = orbits[i]
        if i==0: ax.plot(orbit[:,0],orbit[:,1],color = colors[i],alpha=0.7,label=leg)
        if i>0: ax.plot(orbit[:,0],orbit[:,1],color = colors[i],alpha=0.7)           
        ax.plot(orbit[0,0],orbit[0,1],"*r",markersize=20)
        ax.tick_params(axis='y', which='major', direction='in', labelsize=16)
        ax.tick_params(axis='y', which='minor', direction='in', labelsize=16)
        ax.tick_params(axis='x', which='major', direction='in', labelsize=16)
        ax.tick_params(axis='x', which='minor', direction='in', labelsize=16)
    ax.set_xlabel(r"$X$ [kpc]",size=ls)
    ax.set_ylabel(r"$Y$ [kpc]",size=ls)
#     ax.set_xlabel(r"$X_{\rm GC}$ [kpc]",size=ls)
#     ax.set_ylabel(r"$Y_{\rm GC}$ [kpc]",size=ls)    
    ax.set_xlim((ymin,ymax))
    ax.set_ylim((ymin,ymax))
    ax.legend(loc='upper right',ncol=1,fontsize=18)
    try:
        fig.savefig(kwargs['path']+obj.name+".pdf", bbox_inches='tight')
        cwd = os.getcwd()
        print("fig saved in ",cwd)
    except:
        print("fig NOT saved")
        return
    plt.show
    return

def plotorbplaneXYVlos(obj,**kwargs):
    orbits = obj.orbs[:]
    leg    = f"T={obj.paramtime[3][-1]}[Myr]"
    no = len(orbits)
    print(f"Found {no} orbits")
    varmax,varmin = 0,0
    xmin,xmax,ymin,ymax = 0,0,0,0
    for i in range(0,no):
        orbit  = orbits[i][0]
        var    = -1.*orbit[:,5] # Vlos = -1 Vz
        varmax = np.amax((varmax,np.amax(var)))
        varmin = np.amin((varmin,np.amin(var)))
        xmin,xmax=np.amin((xmin,np.amin(orbit[:,0]))),np.amax((xmax,np.amax(orbit[:,0])))
        ymin,ymax=np.amin((ymin,np.amin(orbit[:,1]))),np.amax((ymax,np.amax(orbit[:,1])))
    frac = 0.1
    ymin,ymax,xmin,xmax = ymin-(ymax-ymin)*frac,ymax+(ymax-ymin)*frac,xmin-(xmax-xmin)*frac,xmax+(xmax-xmin)*frac
    try:
        L = kwargs.get('L', None) #kwargs['L']
        if len(L)>1:
            xmin,xmax,ymin,ymax =  L
        else:
            xmin,xmax,ymin,ymax = -L,L,-L,L
    except Exception:
        pass

    # print('vmax,vmin=',varmax,varmin)
    fig, axs = plt.subplots(1, 1,figsize=(6,6))
    ls = 20    
    ax = axs
    try:
        from matplotlib import image
        print('Image')
        xminI,xmaxI,yminI,ymaxI = kwargs['sizeimage']
        print("kwargs['sizeimage']=",kwargs['sizeimage'],xminI,xmaxI,yminI,ymaxI)
        image = image.imread('Data/Data_Inputs/Data_Images/Malin1v2.jpeg')
        print('Plotting Image')
        ax.imshow(image,extent = ((xminI,xmaxI,yminI,ymaxI)))
    except Exception:
        pass
#     st,va,no = orbits.shape 
#     print('st,va,no =',st,va,no)
    
    
    ## plot data if any
    try:
        xobs,yobs,Robserr= obj.obsdata
        labobs = 'Stream'
        ax.scatter(xobs,yobs,s=45,color='magenta',marker="x",label=labobs,zorder=1E5)
    except Exception:
        pass
    
    evenly_spaced_interval = np.linspace(0.5, 1, no)
    colors = [cm.Purples(x) for x in evenly_spaced_interval]
    
    colormap =  plt.cm.get_cmap('jet')
    al = 0.3
    for i in range(0,no):
        orbit = orbits[i][0]
        var    = -1.*orbit[:,5] # Vlos = -1 Vz        
        # if no==1: orbit = orbits[i]
        # if no>1: orbit = orbits[i]
        if i==0: 
            im2    = ax.scatter(orbit[:,0],orbit[:,1],s=10,c=var,cmap =colormap,alpha=al,vmax=varmax,vmin=varmin,marker="o",label=leg)
        if i>0: 
            im2    = ax.scatter(orbit[:,0],orbit[:,1],s=10,c=var,cmap =colormap,alpha=al,vmax=varmax,vmin=varmin,marker="o")        
    ax.plot(orbit[0,0],orbit[0,1],"*r",markersize=15,alpha=0.45)
    ax.tick_params(axis='y', which='major', direction='in', labelsize=16)
    ax.tick_params(axis='y', which='minor', direction='in', labelsize=16)
    ax.tick_params(axis='x', which='major', direction='in', labelsize=16)
    ax.tick_params(axis='x', which='minor', direction='in', labelsize=16)
    ax.set_xlabel(r"$X$ [kpc]",size=ls)
    ax.set_ylabel(r"$Y$ [kpc]",size=ls)
#     ax.set_xlabel(r"$X_{\rm GC}$ [kpc]",size=ls)
#     ax.set_ylabel(r"$Y_{\rm GC}$ [kpc]",size=ls)    
    ax.set_xlim((xmin,xmax))
    ax.set_ylim((ymin,ymax))
    ax.legend(loc='upper right',ncol=1,fontsize=18)
    
    cbaraxe = fig.add_axes([0.902, 0.125, 0.04, 0.757])  # [left, bottom, width, height]
    sm1 = ScalarMappable(cmap=colormap, norm=Normalize(vmin=varmin, vmax=varmax))
    axcb = fig.colorbar(sm1, cax = cbaraxe)
    axcb.ax.tick_params(labelsize=22, direction='in',which='major', length=8, width=1.5) 
    axcb.ax.tick_params(labelsize=22, direction='in',which='minor', length=4, width=1.2) 
    labelbar = r'$v_{\rm los}\,[{\rm km/s}]$'
    axcb.ax.set_ylabel(labelbar,size=28,labelpad=-6,rotation=90)
    
    try:
        fig.savefig(kwargs['path']+obj.name+".pdf", bbox_inches='tight')
        cwd = os.getcwd()
        print("fig saved in ",cwd)
    except:
        print("fig NOT saved")
        return
    plt.show
    return


def plotsnapshotplaneXY(orbits,snap,L,leg,**kwargs):
    no = len(orbits)
#     st,va,no = orbits.shape 
#     print('st,va,no =',st,va,no)
    evenly_spaced_interval = np.linspace(0.5, 1, no)
    colors = [cm.Purples(x) for x in evenly_spaced_interval]
    
    fig, axs = plt.subplots(1, 1,figsize=(6,6))
    ymin,ymax=-L,L
    ls = 20    
    ax = axs
    markersize = 5
    for i in range(0,no):
        if no==1: orbit = orbits[0]
        if no>1: orbit = orbits[i][0]
        if i==0: ax.plot(orbit[snap,0],orbit[snap,1],"*b",markersize=markersize,alpha=0.7,label=leg)
        if i>0: ax.plot(orbit[snap,0],orbit[snap,1],"*b",markersize=markersize,alpha=0.7)           
        # ax.plot(orbit[0,0],orbit[0,1],"*r",markersize=20)
        ax.tick_params(axis='y', which='major', direction='in', labelsize=16)
        ax.tick_params(axis='y', which='minor', direction='in', labelsize=16)
        ax.tick_params(axis='x', which='major', direction='in', labelsize=16)
        ax.tick_params(axis='x', which='minor', direction='in', labelsize=16)
    ax.set_xlabel(r"$X$ [kpc]",size=ls)
    ax.set_ylabel(r"$Y$ [kpc]",size=ls)
#     ax.set_xlabel(r"$X_{\rm GC}$ [kpc]",size=ls)
#     ax.set_ylabel(r"$Y_{\rm GC}$ [kpc]",size=ls)    
    ax.set_xlim((ymin,ymax))
    ax.set_ylim((ymin,ymax))
    ax.legend(loc='upper right',ncol=1,fontsize=18)
    try:
        fig.savefig(kwargs.get('path')+"fig_orbitXY_"+kwargs.get('name')+".pdf", bbox_inches='tight')
        cwd = os.getcwd()
        print("fig saved in ",cwd)
    except:
        print("fig NOT saved")
        return
    plt.show
    return


def plotorbplanesXYZ(orbits,L,leg,**kwargs):
    no = len(orbits)
#     st,va,no = orbits.shape 
#     print('st,va,no =',st,va,no)
    
    evenly_spaced_interval = np.linspace(0.5, 1, no)
    colors = [cm.Purples(x) for x in evenly_spaced_interval]
    
    # fig, axs = plt.subplots(1, 1,figsize=(6,6))
    fig, axs = plt.subplots(1, 2,figsize=(13,6), subplot_kw={'aspect': 'equal'})
    
    ymin,ymax=-L,L
    ls = 20    
    markersize=5
    for i in range(0,no):
        orbit = orbits[i][0]
        ind = 0
        if i==0: axs[ind].plot(orbit[:,0],orbit[:,1],color = colors[i],alpha=0.7,label=leg)
        if i>0: axs[ind].plot(orbit[:,0],orbit[:,1],color = colors[i],alpha=0.7)           
        ind = 0
        axs[ind].plot(orbit[0,0],orbit[0,1],"*r",markersize=markersize)
        axs[ind].tick_params(axis='y', which='major', direction='in', labelsize=16)
        axs[ind].tick_params(axis='y', which='minor', direction='in', labelsize=16)
        axs[ind].tick_params(axis='x', which='major', direction='in', labelsize=16)
        axs[ind].tick_params(axis='x', which='minor', direction='in', labelsize=16)
        ind = 1        
        if i==0: axs[ind].plot(orbit[:,0],orbit[:,2],color = colors[i],alpha=0.7,label=leg)
        if i>0: axs[ind].plot(orbit[:,0],orbit[:,2],color = colors[i],alpha=0.7)           
        axs[ind].plot(orbit[0,0],orbit[0,2],"*r",markersize=markersize)
        axs[ind].tick_params(axis='y', which='major', direction='in', labelsize=16)
        axs[ind].tick_params(axis='y', which='minor', direction='in', labelsize=16)
        axs[ind].tick_params(axis='x', which='major', direction='in', labelsize=16)
        axs[ind].tick_params(axis='x', which='minor', direction='in', labelsize=16)
        
#     ax.set_xlabel(r"$X_{\rm GC}$ [kpc]",size=ls)
#     ax.set_ylabel(r"$Y_{\rm GC}$ [kpc]",size=ls)    
    # plt.gca().set_aspect('equal')
    axs[0].set_xlabel(r"$X$ [kpc]",size=ls)
    axs[0].set_ylabel(r"$Y$ [kpc]",size=ls)
    axs[0].set_xlim((ymin,ymax))
    axs[0].set_ylim((ymin,ymax))

    axs[1].set_xlabel(r"$X$ [kpc]",size=ls)
    axs[1].set_ylabel(r"$Z$ [kpc]",size=ls,labelpad=-10)
    axs[1].set_xlim((ymin,ymax))
    axs[1].set_ylim((ymin,ymax))
    
    axs[1].legend(loc='upper right',ncol=1,fontsize=18)
    try:
        pathsave = kwargs['path']+"/fig_orbitXYZ_"+kwargs['name']+".pdf"
        fig.savefig(pathsave, bbox_inches='tight')
        print("fig saved in ",pathsave)
    except:
        print("fig NOT saved")
        return
    plt.show
    return

def Planeorb(xi,yi,zi):
    anXY= 34.0526815627
    anZXR= 43.1612255488
    anLZLY= 139.029326425
    zer      = 0.*(xi.copy())
    xR,yR,zR = zer,zer,zer
    xRR,yRR,zRR =  zer,zer,zer
    xRRR,yRRR,zRRR =  zer,zer,zer   
    for i in range(0,len(xi)):
        xR[i],yR[i],zR[i]    = Rotz(xi[i],yi[i],zi[i],-anXY)
        xRR[i],yRR[i],zRR[i]    = Roty(xR[i],yR[i],zR[i],-anZXR)
        xRRR[i],yRRR[i],zRRR[i] = Rotx(xRR[i],yRR[i],zRR[i],-anLZLY)
    return xRRR,yRRR,zRRR

def Planeorb2(xi,yi,zi):
    anXY= 34.0526815627
    anZXR= 43.1612255488
    anLZLY= 139.029326425
    zer      = 0.*(xi.copy())
    xR,yR,zR = zer,zer,zer
    xRR,yRR,zRR =  zer,zer,zer
    xRRR,yRRR,zRRR =  zer,zer,zer
    xR,yR,zR    = Rotz(xi,yi,zi,-anXY)
    xRR,yRR,zRR    = Roty(xR,yR,zR,-anZXR)
    xRRR,yRRR,zRRR = Rotx(xRR,yRR,zRR,-anLZLY)
    return xRRR,yRRR,zRRR

def plotorb(dv,time,scale,ylabel):
    fig, ax = plt.subplots(1, 1,figsize=(12,5))
    d = (dv[:,0]**2+dv[:,1]**2+dv[:,2]**2)**0.5*scale
    ax.plot(time,d,"-b")
#     ax.plot(time,dv[:,10],"--k")
    ls = 20
    ax.set_xlabel("Lookback time[Myr]",size=ls)
    ax.set_ylabel(ylabel,size=ls)
#     ax.set_ylim((0,2.5E2))
#     ax.set_ylim((0.01,4E3))
#     ax.set_yscale('log')
    plt.show
    return

def plotorbvel(dv,time,scale,ylabel):
    fig, ax = plt.subplots(1, 1,figsize=(12,5))
    d = (dv[:,3]**2+dv[:,4]**2+dv[:,5]**2)**0.5
    ax.plot(time,d,"-b")
#     ax.plot(time,dv[:,10],"--k")
    ls = 20
    ax.set_xlabel("Lookback time[Myr]",size=ls)
    ax.set_ylabel(ylabel,size=ls)
#     ax.set_ylim((0,2.5E2))
#     ax.set_ylim((0.01,4E3))
#     ax.set_yscale('log')
    plt.show
    return


def plotatime(time):
    fig, ax = plt.subplots(1, 1,figsize=(8,4))
    ax.plot(time,f_a_t(time),"-b")
    ax.set_yscale('log')
    ax.set_ylabel("a(t)",size=20)
    ax.set_xlabel("lookback time[Myr]",size=20)    
    plt.show
    return

def plotztime(time):
    fig, ax = plt.subplots(1, 1,figsize=(8,4))
    ax.plot(time,f_redshift_t(time),"-b")
#     ax.set_yscale('log')
    ax.set_ylabel("z(t)",size=20)
    ax.set_xlabel("lookback time[Myr]",size=20)    
    plt.show
    return

def plotHtime(time):
    fig, ax = plt.subplots(1, 1,figsize=(8,4))
    ax.plot(time,f_Hkmskpc_t(time),"-b")
#     ax.set_yscale('log')
    ax.set_ylabel("H(t)",size=20)
    ax.set_xlabel("lookback time[Myr]",size=20)    
    plt.show
    return

###############################################################
## FUNCTIONS: TOOLS  ##########################################
###############################################################

def SinkGyr(xi,yi,zi,outputfile):
    fid = open(outputfile, 'w' )
    tosave = np.c_[xi,yi,zi]
    np.savetxt(fid,tosave, fmt ='%0.6e %0.6e %0.6e')
    fid.close()
    return

def f_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

def savelistclassobj(listobj,filename):
    os.system("rm "+filename)
    with open(filename, 'wb') as output:
        for i in range(0,len(listobj)):
            pickle.dump(listobj[i], output, pickle.HIGHEST_PROTOCOL)
    return

def findsnap(timesnap, objclass):
    time = objclass.paramtime[3]
    snapshot = np.argmin(abs(time-timesnap))
    return snapshot

def printdelorean():
    print(" ")
    print("DELOREAN LOADED WITH 1.21GW: GODSPEED USER")    
    print("             _           _           ")
    print("    ________|_\         /_|________  ")
    print("    \_______-''=.     .=''-_______/  ")
    print("               \\\    //             ")
    print("           _____\\\__//_____         ")
    print("         .'                '.        ")
    print("       _'______________.- -._'_      ")
    print("      /  _____________________ \     ")
    print("      |/[_][_]____DMC___[_][_]\|     ")
    print("      \___<_>____________<_>___/     ")
    print("      |\___================___/|     ")
    print("      |__|'                '|__|     ")
    print("                                     ")    
    return
    

def Rotz(xi,yi,zi,anglezi):
    anglezi = np.pi/180.*(anglezi)
    Rz      = np.array([[np.cos(anglezi),-np.sin(anglezi),0],[np.sin(anglezi),np.cos(anglezi),0],[0,0,1]])
    xyzi    = np.array([xi,yi,zi]).T
    xyzRi   = np.dot(Rz,xyzi.T)
    xri, yri, zri = xyzRi[0,:], xyzRi[1,:], xyzRi[2,:]
    return xri,yri,zri

def Rotx(xi,yi,zi,anglexi):
    anglexi = np.pi/180.*(anglexi)
    Rx      = np.array([[1,0,0],[0,np.cos(anglexi),-np.sin(anglexi)],[0,np.sin(anglexi),np.cos(anglexi)]])
    xyzi    = np.array([xi,yi,zi]).T
    xyzRi   = np.dot(Rx,xyzi.T)
    xri, yri, zri = xyzRi[0,:], xyzRi[1,:], xyzRi[2,:]
    return xri,yri,zri

def Roty(xi,yi,zi,angleyi):
    angleyi = np.pi/180.*(angleyi)
    Ry      = np.array([[np.cos(angleyi),0,np.sin(angleyi)],[0,1,0],[-np.sin(angleyi),0,np.cos(angleyi)]])
    xyzi    = np.array([xi,yi,zi]).T
    xyzRi   = np.dot(Ry,xyzi.T)
    xri, yri, zri = xyzRi[0,:], xyzRi[1,:], xyzRi[2,:]
    return xri,yri,zri

def Rotz_v2(xi,yi,zi,anglezi):
    anglezi = np.pi/180.*(anglezi)
    Rz      = np.array([[np.cos(anglezi),-np.sin(anglezi),0],[np.sin(anglezi),np.cos(anglezi),0],[0,0,1]])
    xyzi    = np.array([xi,yi,zi]).T
    xyzRi   = np.dot(Rz,xyzi.T)
    # print('xyzRi.shape=',xyzRi.shape)
    xri, yri, zri = xyzRi[0], xyzRi[1], xyzRi[2]
    return xri,yri,zri

def Rotx_v2(xi,yi,zi,anglexi):
    anglexi = np.pi/180.*(anglexi)
    Rx      = np.array([[1,0,0],[0,np.cos(anglexi),-np.sin(anglexi)],[0,np.sin(anglexi),np.cos(anglexi)]])
    xyzi    = np.array([xi,yi,zi]).T
    xyzRi   = np.dot(Rx,xyzi.T)
    xri, yri, zri = xyzRi[0], xyzRi[1], xyzRi[2]
    return xri,yri,zri

def Roty_v2(xi,yi,zi,angleyi):
    angleyi = np.pi/180.*(angleyi)
    Ry      = np.array([[np.cos(angleyi),0,np.sin(angleyi)],[0,1,0],[-np.sin(angleyi),0,np.cos(angleyi)]])
    xyzi    = np.array([xi,yi,zi]).T
    xyzRi   = np.dot(Ry,xyzi.T)
    xri, yri, zri = xyzRi[0], xyzRi[1], xyzRi[2]
    return xri,yri,zri

def Rotysc(xi,yi,zi,angleyi):
    angleyi = np.pi/180.*(angleyi)
    Ry      = np.array([[np.cos(angleyi),0,np.sin(angleyi)],[0,1,0],[-np.sin(angleyi),0,np.cos(angleyi)]])
    xyzi    = np.array([xi,yi,zi]).T
    xyzRi   = np.dot(Ry,xyzi.T)
    xri, yri, zri = xyzRi[0], xyzRi[1], xyzRi[2]
    return xri,yri,zri

# def Rotz(xi,yi,zi,anglezi):
#     anglezi = np.pi/180.*(anglezi)
#     Rz      = np.array([[np.cos(anglezi),-np.sin(anglezi),0],[np.sin(anglezi),np.cos(anglezi),0],[0,0,1]])
#     xyzi    = np.array([xi,yi,zi]).T
#     xyzRi   = np.dot(Rz,xyzi.T)
#     xri, yri, zri = xyzRi[0,:], xyzRi[1,:], xyzRi[2,:]
#     return xri,yri,zri

# def Rotx(xi,yi,zi,anglexi):
#     anglexi = np.pi/180.*(anglexi)
#     Rx      = np.array([[1,0,0],[0,np.cos(anglexi),-np.sin(anglexi)],[0,np.sin(anglexi),np.cos(anglexi)]])
#     xyzi    = np.array([xi,yi,zi]).T
#     xyzRi   = np.dot(Rx,xyzi.T)
#     xri, yri, zri = xyzRi[0,:], xyzRi[1,:], xyzRi[2,:]
#     return xri,yri,zri

# def Roty(xi,yi,zi,angleyi):
#     angleyi = np.pi/180.*(angleyi)
#     Ry      = np.array([[np.cos(angleyi),0,np.sin(angleyi)],[0,1,0],[-np.sin(angleyi),0,np.cos(angleyi)]])
#     xyzi    = np.array([xi,yi,zi]).T
#     xyzRi   = np.dot(Ry,xyzi.T)
#     xri, yri, zri = xyzRi[0,:], xyzRi[1,:], xyzRi[2,:]
#     return xri,yri,zri

################################################
## FUNCTIONS: DISTANCES ########################
def PdistDeg2AU(pddist_obj_deg,rdist_obj_kpc) :
    kpc2AU  = 2.063e+8 #(kpc/AU)
    deg2rad = np.pi/180.
    rdist_obj_AU = rdist_obj_kpc*kpc2AU
    pddist_obj_AU = 2.*rdist_obj_AU*np.tan(pddist_obj_deg*deg2rad/2)
    print("Pdist=",pddist_obj_AU,"[AU]")
    return

def PdistDeg2kpc(pddist_obj_deg,rdist_obj_kpc) :
#     kpc2AU  = 2.063e+8 #(kpc/AU)
    deg2rad = np.pi/180.
#     rdist_obj_AU = rdist_obj_kpc*kpc2AU
    pddist_obj_kpc = 2.*rdist_obj_kpc*np.tan(pddist_obj_deg*deg2rad/2)
    print("Pdist=",pddist_obj_kpc,"[kpc]")
    return

################################################
## FUNCTIONS: PHOTOMETRY #######################
def Mag2Lum(Mvin):
    MvSun = 4.87
    Lum   = 10.**((MvSun-Mvin)/2.5)
    print("Lv=",'%.2E' % Decimal(str(Lum)),"Lsun")
    return Lum

def mag2Mag(magvin,dist):
    MagV   = magvin - 5.*np.log10(dist*1000./10)
    print("MagV=",'%.2E' % Decimal(str(MagV)),"mag")
    return MagV

def Lum2Mag(Lvin):
    MvSun = 4.87
    Mvout = MvSun - 2.5*np.log10(Lvin) 
    print("Mv=",str(round(Mvout,2)),"mag")
    return Mvout


################################################
## FUNCTIONS: COSMOLOGY  #######################
def deltaf(zin):
    xf = Planck15.Om(zin)-1
    delta=(18*np.pi**2+82.*xf-39*xf**2)/(xf+1)
    return delta

def rvirf(Mvir_in,zin):  
#     rvir_out=258.*(340.*0.3/102)**(-1./3)*(Mvir_in/1.E12)**(1./3)
#     redshift= 0
    rvir_out = (Mvir_in)**(1./3)/(4.*np.pi/3.*Planck15.critical_density(zin).to(u.solMass/u.kpc**3).value*Planck15.Om(zin)*deltaf(zin))**(1./3)
    return rvir_out

def Rvirf(Mvir_in,zin):  # from Bland-Hawthorn & Gerhard 2016
#     rvir_out=258.*(340.*0.3/102)**(-1./3)*(Mvir_in/1.E12)**(1./3)
#     redshift= 0
    Rvir_out = rvirf(Mvir_in,zin)/af(zin)
    return Rvir_out

def af(zin):
    a = 1./(1.+zin)
    return a

def looktime2z(lookbacktime_in):
    from astropy.cosmology import Planck15, z_at_value
    zredshift_out = np.zeros((len(lookbacktime_in),))
    for i in range(0,len(lookbacktime_in)): zredshift_out[i]= z_at_value(Planck15.lookback_time, (lookbacktime_in[i]))
    return zredshift_out

def looktime2zEq(lookbacktime_in):
    # VALID ONLY FOR HIGH REDSHIFT z>2 !!
    from astropy.cosmology import Planck15, z_at_value
    AgeU = Planck15.age(0).to(u.Myr)
    H0_15, Om0_15 = ((67.7* (1.*u.km).to(u.Mpc)) / u.Mpc/ (1.*u.s)).to(1./u.Myr), 0.307
    zredshift_out = np.zeros((len(lookbacktime_in),))
    zredshift_out = (2./(3.*H0_15*Om0_15**0.5*(AgeU-lookbacktime_in)))**(2./3.)-1
    return zredshift_out

################################################
## FUNCTIONS: DYNAMICS #########################
def rvir_BG16f(Mvir_in):  # from Bland-Hawthorn & Gerhard 2016
    Omega = 0.3
    rvir_out=258.*(340.*Omega/102)**(-1./3)*(Mvir_in/1.E12)**(1./3)
    return rvir_out

def tff(Mi,Ri):
    rho = 3.*Mi/Ri**3/4./np.pi
    tff_out = (3.*np.pi/32./G4/rho)**(1./2)
    print("tff=",tff_out,"Myr")
    return tff_out

def tff_AU(Mi,Ri) :
    G_units1 = 39.478 # AU/Msun(AU/yr)^2
    rho      = 3.*Mi/Ri**3/4./np.pi
    tff_out  = (3.*np.pi/32./G_units1/rho)**(1./2)
    print("tff=",tff_out,"yr")
    return tff_out

def rtidal(mi,Mi,Di):
    rtidal = Di*(mi/Mi/3.)**(1./3)
#     print("Rt=",rtidal
    return rtidal

def rtidal_MWnfw(mi,Di):
    Mvir = 1.3E+12
    rvir = 282.    #kpc
    c    = 11.5 
    rs   = rvir/c
    Om   = 0.3
    deno = Mvir/(4.*np.pi*rs**3*(np.log(1.+rvir/rs)-rvir/rs/(1.+rvir/rs)))
    rtidal = Di*(mi/Mvir/3.)**(1./3)
    mask = Di<rvir
    Mi   = 4.*np.pi*deno*rs**3*(np.log(1.+Di[mask]/rs)-Di[mask]/rs/(1.+Di[mask]/rs))
    rtidal[mask] = Di[mask]*(mi/Mi/3.)**(1./3)
    return rtidal

def Mvc(vcin,rin):
    Mvco=rin*vcin**2/Gpc
    print("Mvco("+str(rin)+"pc)=",'%.2E' % Decimal(str(Mvco)),"M")
    return Mvco


################################################
## FUNCTIONS: MASS PROFILES ####################
def Mnfw(r_in,rs_in,rho_in):
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i] = 4.*np.pi*rho_in*(rs_in**3.)*(np.log(1.0+r_in[i]/rs_in)-(r_in[i]/(rs_in+r_in[i])))
    return Mout

def Mnfw_MW(r_in):
    Mout = np.zeros((len(r_in)))
    Mvir = 1.3E+12
    rvir = 282.    #kpc
    c    = 11.5
    rs   = rvir/c
#     Om   = 0.3
    deno = Mvir/(4.*np.pi*rs**3*(np.log(1.+rvir/rs)-rvir/rs/(1.+rvir/rs)))
    Mout   = 4.*np.pi*deno*rs**3*(np.log(1.+r_in/rs)-r_in/rs/(1.+r_in/rs))    
    return Mout

def Mnfw_rho(r_in,rho_in):
    rs_in = 1.0 #kpc
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i] = 4.*np.pi*rho_in*(rs_in**3.)*(np.log(1.0+r_in[i]/rs_in)-(r_in[i]/(rs_in+r_in[i])))
    return Mout

def LogMnfw(r_in,rs_in,rho_in):
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i] = 4.*np.pi*rho_in*(rs_in**3.)*(np.log(1.0+r_in[i]/rs_in)-(r_in[i]/(rs_in+r_in[i])))
    return np.log10(Mout)


def MpsisoN(r_in,rs_in,rho_in,nrel):
    from scipy.interpolate import interp1d
    r_min = 0.5*np.amin(r_in)
    r_max = 1.5*np.amax(r_in)
    rstp  = (r_max-r_min)/(nrel*len(r_in))
    r_i   = np.arange(r_min,r_max+rstp,rstp)
    rho_i = rho_in/(1.+(r_i/rs_in)**2)
    Mshell_i = 4.*np.pi*rho_i*r_i**2*rstp
    Mout_i = np.zeros((len(r_i)))
    counter = 0.
    for i in range(0,len(r_i)):
        Mout_i[i] = counter + Mshell_i[i]
        counter  = Mout_i[i]
    Mout = np.zeros((len(r_in)))
    Mf = interp1d(r_i, Mout_i)
    Mout = Mf(r_in)
    return Mout

def Mpsiso(r_in,rs_in,rho_in):
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i] = 4.*np.pi*rho_in*rs_in**2*(r_in[i]-rs_in*np.arctan(r_in[i]/rs_in))
    return Mout

def LogMpsiso(r_in,rs_in,rho_in):
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i] = 4.*np.pi*rho_in*rs_in**2*(r_in[i]-rs_in*np.arctan(r_in[i]/rs_in))
    return np.log10(Mout)


def McoredW09(r_in,rs_in,rho_in):
    Cte = 4.*np.pi*rho_in*rs_in**3*((rs_in*(3.*rs_in+4.*0.))/(2.*(rs_in+0.)**2)+np.log(rs_in+0.))
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i]=4.*np.pi*rho_in*rs_in**3*((rs_in*(3.*rs_in+4.*r_in[i]))/(2.*(rs_in+r_in[i])**2)+np.log(rs_in+r_in[i])) -Cte
    return Mout

def LogMcoredW09(r_in,rs_in,rho_in):
    Cte = 4.*np.pi*rho_in*rs_in**3*((rs_in*(3.*rs_in+4.*0.))/(2.*(rs_in+0.)**2)+np.log(rs_in+0.))
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i]=4.*np.pi*rho_in*rs_in**3*((rs_in*(3.*rs_in+4.*r_in[i]))/(2.*(rs_in+r_in[i])**2)+np.log(rs_in+r_in[i])) -Cte
    return np.log10(Mout)


def Mburk(r_in,rs_in,rho_in):
#     Cte = 4.*np.pi*rho_in*rs_in**3.*(1./4.*np.log(rs_in**2+0.**2)+1./2.*np.log(rs_in+0.)-1./2.*np.arctan(0./rs_in))
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i]=4.*np.pi*rho_in*rs_in**3.*(1./4.*np.log(1.0+(r_in[i]**2)/rs_in**2)+1./2.*np.log(1.0+r_in[i]/rs_in)-1./2.*np.arctan(r_in[i]/rs_in)) #- Cte
    return Mout

def Mburk_rho(r_in,rho_in):
    rs_in = 0.25 #kpc
#     Cte = 4.*np.pi*rho_in*rs_in**3.*(1./4.*np.log(rs_in**2+0.**2)+1./2.*np.log(rs_in+0.)-1./2.*np.arctan(0./rs_in))
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i]=4.*np.pi*rho_in*rs_in**3.*(1./4.*np.log(1.0+(r_in[i]**2)/rs_in**2)+1./2.*np.log(1.0+r_in[i]/rs_in)-1./2.*np.arctan(r_in[i]/rs_in)) #- Cte
    return Mout


def LogMburk(r_in,rs_in,rho_in):
#     Cte = 4.*np.pi*rho_in*rs_in**3.*(1./4.*np.log(rs_in**2+0.**2)+1./2.*np.log(rs_in+0.)-1./2.*np.arctan(0./rs_in))
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i]=4.*np.pi*rho_in*rs_in**3.*(1./4.*np.log(1.0+r_in[i]**2/rs_in**2)+1./2.*np.log(1.0+r_in[i]/rs_in)-1./2.*np.arctan(r_in[i]/rs_in)) #- Cte
    return np.log10(Mout)


####################################################################
########## OTHER POTENTIAL GRADIENTS ###############################
def Mburk(r_in,rs_in,rho_in):
#     Cte = 4.*np.pi*rho_in*rs_in**3.*(1./4.*np.log(rs_in**2+0.**2)+1./2.*np.log(rs_in+0.)-1./2.*np.arctan(0./rs_in))
    Mout = np.zeros((len(r_in)))
    for i in range(0,len(r_in)):
        Mout[i]=4.*np.pi*rho_in*rs_in**3.*(1./4.*np.log(1.0+(r_in[i]**2)/rs_in**2)+1./2.*np.log(1.0+r_in[i]/rs_in)-1./2.*np.arctan(r_in[i]/rs_in)) #- Cte
    return Mout


def plum_LeoT_ac(rvin):
    Mpl_gas  = 5.5E5
    rpl_gas  = 0.19394
    Mpl_st   = 2.0E5 
    rpl_st   = 0.145 
    rin  = (np.dot(rvin,rvin))**0.5
    acc   = -G4*Mpl_gas/(rin**2+rpl_gas**2)**1.5*rvin + -G4*Mpl_st/(rin**2+rpl_st**2)**1.5*rvin
    return acc #kpc/Myr^2


def acc_burkert(rvin):
#     rs_in  = 0.4  # Burkert Fidiucial HI Patra 2018
#     rho_in = 1.089372E+08 # Burkert Fidiucial HI Patra 2018
    
    rs_in  = 0.4
    rho_in =  2.262736E+07
    
#     rs_in  = 0.4  # Burkert Massive HI data
#     rho_in = 1.336446E+08 # Burkert Massive HI data
    
#     rs_in = 0.4 #kpc
#     rho_in= 1.240247E+08 # Msun/kpc^3

#     rs_in = 0.07945 #kpc
#     rho_in= 2.268142E+09 # Msun/kpc^3 
    
    rinf  = 0.0001 # =0.1pc
    rin  = (np.dot(rvin,rvin))**0.5
    if rin<=rinf: rin = rinf
    
    MbR   = 4.*np.pi*rho_in*rs_in**3.*1./2*(1./2.*np.log(1.+rin**2/rs_in**2) + np.log(1.+rin/rs_in) - np.arctan(rin/rs_in)) #- Cte
     
    acc   = -G4*MbR/(rin**3)*rvin + plum_LeoT_ac(rvin)
    return acc #kpc/Myr^2

def acc_burkert2(rvin):
    rs_in  = 0.4  # Burkert Fidiucial HI Patra 2018
#     rho_in = 1.089372E+08 # Burkert Fidiucial HI Patra 2018
    rho_in = 1.33644648e+08  #Model S6M3 Massive Burkert
#     rs_in  = 0.4  # Burkert Massive HI data
#     rho_in = 1.336446E+08 # Burkert Massive HI data
    
#     rs_in = 0.4 #kpc
#     rho_in= 1.240247E+08 # Msun/kpc^3

#     rs_in = 0.07945 #kpc
#     rho_in= 2.268142E+09 # Msun/kpc^3 
    
    rinf  = 0.0001 # =0.1pc
    rin  = (np.dot(rvin,rvin))**0.5
    if rin<=rinf: rin = rinf
    
    MbR   = 4.*np.pi*rho_in*rs_in**3.*1./2*(1./2.*np.log(1.+rin**2/rs_in**2) + np.log(1.+rin/rs_in) - np.arctan(rin/rs_in)) #- Cte
     
    acc   = -G4*MbR/(rin**3)*rvin + plum_LeoT_ac(rvin)
    return acc #kpc/Myr^2


def Pacc_burkert(rvin):
    n,m = rvin.shape
    acc = np.zeros((n,m))
    for i in range(0,n):
        acc[i,:] = acc_burkert(rvin[i,:])
    return acc

def acc_nfw(rvin):
    r_in   = (np.dot(rvin,rvin))**0.5
    rs_in  = 0.10991988914315304
    rho_in = 1.157707E+09
    MnfwR  = 4.*np.pi*rho_in*(rs_in**3.)*(np.log(1.0+r_in/rs_in)-(r_in/(rs_in+r_in)))
#     cin  = 8.74603379
#     rvir = 266.3    #kpc
#     rs   = rvir/cin
#     Mvir = 1.1E+12 #Ms 0.8-1.1 x10^10Msun Tamm 2012, 
#     Om   = 0.3
#     A    = -G4*Mvir/(np.log(1.+cin)-cin/(1.+cin))
#     acc  = -A*(rvin/rin*(1./(rin*(rin+rs)))-np.log(1.+rin/rs)/rin**2*rvin/rin)
#     if (rin>=rvir): acc =-G4*Mvir/(rin*rin)*(rvin/rin) #(kpc/Myr)^2/kpc
    acc  = -G4*MnfwR/r_in**3*rvin + plum_LeoT_ac(rvin)
    return acc #kpc/Myr^2


def accdynfr(dtin,vvin,rvin,Mvir,rvir,cin):
    rin  = (np.dot(rvin,rvin))**0.5
    if (rin>rvir):
        acc = np.zeros((len(rvin)))
    else:
        vin  = (np.dot(vvin,vvin))**0.5  # from Binney+Trimaine-8
        rh   =   rvir/cin
#         mass_sat = 1.E+8 #minimum to form stars (Read+06)
#         rvir_sat = 10.
        mass_sat = 1.E+10 #LMC ?
        rvir_sat = 50.         
        rh_sat   = rvir_sat/10
#         LAM      = (Mvir/mass_sat)*(rin/(rvir))
        LAM      = (Mvir/mass_sat)*(rin/(rh))
#         LAM      = rh/rh_sat
        sigma    = 200.*kms2kpcMyr
        X        = vin/(2.**0.5*sigma)
        lnLAM    = abs(np.log(LAM))
        den      = denMW(rvin,Mvir,rvir,cin)
        acc      = (-4.)*np.pi*G4**2*lnLAM*den*mass_sat/vin**3*(special.erf(X)-2.*X/np.pi**0.5*np.e**(-X**2))*vvin
#         acc      = (-4.)*np.pi*G4**2*lnLAM*den*mass_sat*vvin/vin**3
#         acc      = np.sign(dtin)*(-4.)*np.pi*G4**2*lnLAM*den*mass_sat*vvin/vin**3
    return acc #kpc/Myr^2

###################################################################################
#### To use GyrfalcON in NEMO to calculate accelerations from N-body model ########
def accnbody(rvin):
#     filesim     = "MWpart1509"
#     filesim     = "diceMW1_t2120" #"MWpart1509"
    foldernbody = "Data/Data_Nbody/"
    filesim     = foldernbody+"MW1_t2120v3"
    filegrid    = foldernbody+"sinkgrid"
    fileaccgrid = foldernbody+"accgrid"
    x_saveGrid,y_saveGrid,z_saveGrid = (rvin.T)*kpc2ud
    os.system("rm "+filegrid+".a1")
    SinkGyr(x_saveGrid,y_saveGrid,z_saveGrid,filegrid+".a1")
    os.system("rm "+filegrid+".s1")
    os.system("a2s in="+filegrid+".a1 out="+filegrid+".s1 N=1 read=x")
    os.system("rm "+fileaccgrid+".a1")
    os.system("rm "+fileaccgrid+".s1")    
    os.system("getgravity srce="+filesim+".s1 sink="+filegrid+".s1 theta=0.5 eps=0.05 out="+fileaccgrid+".s1 Ncrit=6")
    os.system("snapprint in="+fileaccgrid+".s1 options=ax,ay,az format=%g separ=0 times=0 > "+fileaccgrid+".a1")
    datacc = np.loadtxt(fileaccgrid+".a1")
    acx,acy,acz = (datacc.T)*uac2kpcMyrsq
    acc = np.array([acx,acy,acz])
    return acc


##################################################################################################################
## Densities for calculations of the dynamical friction F_dyn ~ den/v^2 and Ram Pressure Stripping (RPS) #########
##################################################################################################################

def denhalonfwnostars_MW_ac(rvin,Mvir,rvir,cin):
    rs   = rvir/cin
    Om   = 0.3
    Mdmg = (1.0 - GLfrac_Mbulge - GLfrac_Mdisk)*Mvir  #Ms    
    deno = 3.*Mdmg/(4.*np.pi*rvir**3*Om)    
    rin  = (np.dot(rvin,rvin))**0.5
    if (rin<=rvir):
        den_halo = deno/(rin/rs)/(1.+rin/rs)**2
    elif (rin>rvir):
        den_halo = 0. 
    return den_halo #Msun/kpc^3

def dendisc_MW_ac(rvin,Mviri):
    zh   = 0.3 #kpc #Ms BO15w
    Rd   = 2.5 #kpc #Ms BO15
    Md   = GLfrac_Mdisk*Mviri #Ms BO15 = stellar disk+ gas disk
    xi   = rvin[0]
    yi   = rvin[1]
    zi   = rvin[2]
    ri   = (xi**2 + yi**2 + zi**2)**0.5
    Ri   = (xi**2 + yi**2)**0.5
    C    = Md * zh**2 *(Rd*Ri**2 + (Rd + (3*(zi**2+zh**2)**0.5)*(Rd+(zi**2+zh**2)**0.5)**2))
    D    = 4.*np.pi * (Ri**2+(Rd+(zi**2+zh**2)**0.5)**2)**(5./2)*(zi**2+zh**2)**(3./2)
    den_disk = C/D  
    return den_disk #Msun/kpc^3

def denbulgeplum_MW_ac(rvin,Mviri):
    Mpl   = GLfrac_Mbulge*Mviri
    rpl   = 0.3      # BO2015 (0-25%CB)
    rin  = (np.dot(rvin,rvin))**0.5
    den_plum   = 3./4/np.pi*Mpl/rpl**3/(1.+(rin/rpl)**2)**2.5
    return den_plum #Msun/kpc^3

def denMW(rvin,Mvir,rvir,cin):
    den_out = denhalonfwnostars_MW_ac(rvin,Mvir,rvir,cin)+dendisc_MW_ac(rvin,Mvir)+denbulgeplum_MW_ac(rvin,Mvir)
    return den_out

# def MRCvirfun(time_in):
#     Mvir_out = f_Mvir_t(time_in)
#     c_out    = f_c_t(time_in)
#     rvir_out = f_Rvir_t(time_in)
#     return Mvir_out,c_out,rvir_out

#####################################################
## STILL USED FOR THE MW - M31 ORBIT PRE-COMPUTATION 
#####################################################
def MRcfun(objclass,case,t_in):
    # return Mviri,Rviri,ci
    global GLfrac_Mdisk, GLfrac_Mbulge    
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

# ######################################################################################
# ### SETUPS or CASES FOR DIFFERENT POTENTIAL GRADIENTS FOR ACCELERATIONS ##############
# ######################################################################################
# def accase(objclass,case,ti,dto,ri,vi,riPOT1,viPOT1,riPOT2,Mviri,Rviri,ci):    
#     if case=="case1" or case=="case1a" or case=="case1b" or case=="case2" or case=="case2a" or case=="case2b":
#         acc = acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
# #     elif case=="case3" or case=="case3a" or case=="case3b" or case=="case4" or case=="case4a" or case=="case4b" or case=="case3_cosmo" or case=="case4_cosmo":
# #         acc = acchalonfw_M31(ri-riPOT2)+acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
#     elif case=="case3_cosmo" or case=="case4_cosmo":
#         acc = acchalonfwnostars_MW_ac(ri-riPOT2,Mviri,Rviri,ci)+acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
#     elif case=="case3" or case=="case3a" or case=="case3b" or case=="case4" or case=="case4a" or case=="case4b":
#         acc = acchalonfw_M31(ri-riPOT2)+acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
#     elif case=="case5" or case=="case6" or case=="case5fut":
#         acc = acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri) + accdynfr(dto,vi-viPOT1,ri-riPOT1,Mviri,Rviri,ci)
#     elif case=="case1fut":
#         acc = acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
#     elif case=="caseCOtest":
#         acc = acchalonfwnostars_MW_ac(ri-riPOT1)
#     elif case=="case1nbody":
#         acc = accnbody(ri-riPOT1)
#     elif case=="case1_cosmo" or case=="case2_cosmo": 
#         acc = acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)+accdisc_MW_ac(ri-riPOT1,Mviri)+accbulgeplum_MW_ac(ri-riPOT1,Mviri)
#     elif case=="case1_kepler":
#         acc = acchalokeplernostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)
#     elif case=="case1_burk":
#         acc = acc_burkert(ri-riPOT1)
#     elif case=="case2_burk":
#         acc = acc_burkert2(ri-riPOT1)    
#     elif case=="case1_nfw":
#         acc = acc_nfw(ri-riPOT1)  
#     elif case=="Fornax_Setup1":
#         acc = acc_Fornax_nfw(ri-riPOT1)  
#     elif case=="Fornax_Setup2":
#         acc = acc_Fornax_nfw_ac(ri-riPOT1,Mviri,Rviri,ci)
#     elif case=="Malin1_Setup1":
#         acc = acc_Malin1(ri-riPOT1,Mviri,Rviri,ci)
#     elif case=="Malin1_Setup1i":
#         acc = acc_Malin1(ri-riPOT1,Mviri,Rviri,ci,i=35.,PA=90.)
#     elif case=="Malin1_Setup2":
#         acc = acc_Malin1(ri-riPOT1,Mviri,Rviri,ci)
#     elif case=="Malin1_Setup3":
#         acc = acc_Malin1(ri-riPOT1,Mviri,Rviri,ci)     #  rvin-riPOT1      
#     elif case=="bar_case1":
#         # acc = accbar(ri,ti)   +acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)  #  rvin-riPOT1
#         acc = accbarLM(ri,ti)   +acchalonfwnostars_MW_ac(ri-riPOT1,Mviri,Rviri,ci)  #  rvin-riPOT1   
#     return acc
