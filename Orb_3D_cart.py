"""This script plots and animates the orbits of globular clusters (GC).

Orbits are calculated backwards in time in a Milky Way potential model.
It uses packages from the DELOREAN_v1 file and animates trajectories
for selected clusters.
"""

import sys
import numpy as np
sys.path.append('Code')

from matplotlib import pyplot as plt
from matplotlib import animation
from DELOREAN_v1 import StarCluster, MySetup, ExploreParam_v2


###########################################################
# Globular clusters of interest: NGC 5694, NGC 5824,
# NGC 6229, NGC 7006, NGC 7492, Pal 13, Whiting 1 y Munoz 1
#
# CLUSTER        RA             DEC			l          b
#  NAME          [°]            [°]		   [°]        [°]
# NGC 5694   219.901245		-26.538776	 331.056 	 30.360
# NGC 5824   225.994156 	-33.068138 	 332.555 	 22.071
# NGC 6229   251.745249		+47.527796 	  73.639 	 40.306
# NGC 7006   315.372626		+16.187323 	  63.770	-19.407
# NGC 7492	 347.111170 	-15.611469 	  53.386 	-63.478
# Pal 13  	 346.685190 	 12.771539 	  87.103 	-42.700
# Whiting 1   30.737499 	 -3.252778 	 161.618 	-60.636
# Munoz 1	 15.0300056     66.9686944
#
# CLUSTER        R☉             R_GC             V_r
#  NAME         [kpc]           [kpc]          [km/sec]
# NGC 5694   34.84 ± 0.74	29.16 ± 0.70 	-139.55 ± 0.49
# NGC 5824   31.71 ± 0.60	25.42 ± 0.57 	 -25.24 ± 0.52
# NGC 6229   30.11 ± 0.47	29.45 ± 0.44 	-137.89 ± 0.71
# NGC 7006   39.32 ± 0.56	36.67 ± 0.53 	-383.47 ± 0.73
# NGC 7492	 24.39 ± 0.57	23.57 ± 0.52 	-176.70 ± 0.27
# Pal 13  	 23.48 ± 0.40	24.57 ± 0.37 	  25.30 ± 0.22
# Whiting 1  30.59 ± 1.17	35.15 ± 1.11 	-130.41 ± 1.79
# Munoz 1	    45 ± 5 	 	 	 	 	 	   −137 ± 4
#
# CLUSTER       μ_α*cosδ          μ_δ        ρ_{μ_α}_{μ_δ}
#  NAME         [mas/yr]        [mas/yr]
# NGC 5694   -0.476 ± 0.012	 -1.102 ± 0.011	     -0.36
# NGC 5824   -1.193 ± 0.010	 -2.228 ± 0.009	     -0.20
# NGC 6229   -1.156 ± 0.017	 -0.461 ± 0.018	      0.17
# NGC 7006   -0.115 ± 0.010	 -0.619 ± 0.009	      0.03
# NGC 7492	  0.777 ± 0.011	 -2.321 ± 0.011	      0.21
# Pal 13  	  1.740 ± 0.039	  0.116 ± 0.032	      0.11
# Whiting 1  -0.244 ± 0.051	 -2.019 ± 0.042	      0.06
# Munoz 1
#
# CLUSTER          X               Y              Z
#  NAME          [kpc]           [kpc]          [kpc]
# NGC 5694   -18.13 ± 0.56 	-14.55 ± 0.31 	 17.61 ± 0.37
# NGC 5824   -17.90 ± 0.49 	-13.55 ± 0.26 	 11.92 ± 0.23
# NGC 6229     1.71 ± 0.10 	 22.03 ± 0.34 	 19.47 ± 0.30
# NGC 7006    -8.21 ± 0.23 	 33.27 ± 0.47 	-13.06 ± 0.19
# NGC 7492	   1.68 ± 0.15 	  8.74 ± 0.21 	-21.82 ± 0.51
# Pal 13  	   7.31 ± 0.01 	 17.23 ± 0.29 	-15.92 ± 0.27
# Whiting 1   22.41 ± 0.54 	  4.73 ± 0.18 	-26.66 ± 1.02
# Munoz 1
#
# CLUSTER         V_x             V_y             V_z
#  NAME         [km/sec]        [km/sec]        [km/sec]
# NGC 5694   114.37 ± 1.32 	 144.79 ± 1.72 	-170.18 ± 1.60
# NGC 5824    97.40 ± 0.96 	 -61.59 ± 1.35 	-179.17 ± 1.35
# NGC 6229    -7.59 ± 2.53 	  30.75 ± 1.77 	  45.79 ± 1.93
# NGC 7006    62.61 ± 1.61	-134.21 ± 1.07 	  84.85 ± 1.69
# NGC 7492	  -1.29 ± 1.26 	 -75.12 ± 1.16 	  70.68 ± 0.62
# Pal 13  	 161.51 ± 4.18 	 216.80 ± 2.57 	 -77.87 ± 2.82
# Whiting 1 -246.96 ± 6.29 	  32.88 ± 6.60 	  -7.97 ± 3.57
# Munoz 1
#
# CLUSTER       R_Peri          R_Apo
#  NAME         [kpc]           [kpc]
# NGC 5694   29.14 ± 0.71	 59.26 ± 1.47
# NGC 5824   13.94 ± 1.19	 34.99 ± 1.61
# NGC 6229    1.43 ± 0.20	 30.36 ± 0.40
# NGC 7006    2.33 ± 0.18	 50.56 ± 0.78
# NGC 7492	  2.81 ± 0.38	 26.03 ± 0.64
# Pal 13  	  6.91 ± 0.34	 58.24 ± 1.76
# Whiting 1  35.15 ± 1.12	 61.94 ± 4.86
# Munoz 1
###########################################################

######################################################
# EXAMPLE 1: Calculate orbits in a Milky Way potential
######################################################

testorbit1 = StarCluster()
testorbit1.name = 'Testorbit1'
testorbit1.ncpu = 1
testorbit1.setup = MySetup(setup_name='Setup:Blana2020',
                           case='case1').setup
testorbit1.coordtype = 'Cartesian'
# testorbit1.orbprop   = True
# testorbit1.paramtime = [0, -5000, 0.1]
testorbit1.paramtime = [0, -5000, 0.1]  # t_init, t_fin, dt
x, y, z, vx, vy, vz = -18.13, -14.55, 17.61, 114.37, 144.79, -170.18
testorbit1.paramvars = [[x, y, z, vx, vy, vz]]
x, y, z, vx, vy, vz = -17.90, -13.55, 11.92, 97.40, -61.59, -179.17
testorbit1.paramvars.append([x, y, z, vx, vy, vz])
x, y, z, vx, vy, vz = 1.71, 22.03, 19.47, -7.59, 30.75, 45.79
testorbit1.paramvars.append([x, y, z, vx, vy, vz])
x, y, z, vx, vy, vz = -8.21, 33.27, -13.06, 62.61, -134.21, 84.85
testorbit1.paramvars.append([x, y, z, vx, vy, vz])
x, y, z, vx, vy, vz = 1.68, 8.74, -21.82, -1.29, -75.12, 70.68
testorbit1.paramvars.append([x, y, z, vx, vy, vz])
x, y, z, vx, vy, vz = 7.31, 17.23, -15.92, 161.51, 216.80, -77.87
testorbit1.paramvars.append([x, y, z, vx, vy, vz])
x, y, z, vx, vy, vz = 22.41, 4.73, -26.66, -246.96, 32.88, -7.97
testorbit1.paramvars.append([x, y, z, vx, vy, vz])
listobjects_init = [testorbit1]
listobjects_out = ExploreParam_v2(listobjects_init)
testorbit1 = listobjects_out[0]
# plotorbplaneXY(testorbit1,L=100,path='Data/Data_Outputs/Figs/fig_test_3NGC')

N_trajectories = len(testorbit1.orbs)

# Solve for the trajectories
t = testorbit1.paramtime[3]
x_t0 = []

for i in range(len(testorbit1.orbs)):
    x_t0.append(testorbit1.orbs[i][0][:, :3].tolist())

x_t = np.array(x_t0)

# Set up figure & 3D axis for animation
fig = plt.figure(figsize=(15, 15))
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
ax.clear()
ax.axis('on')

# Choose a different color for each trajectory
colors = plt.cm.jet(np.linspace(0, 1, N_trajectories))
cluster_name = ['NGC 5694', 'NGC 5824', 'NGC 6229', 'NGC 7006',
                'NGC 7492', 'Pal 13', 'Whiting 1']

# Set up lines and points
lines = sum([ax.plot([], [], [], '-', c=c)
             for c in colors], [])
pts = sum([ax.plot([], [], [], 'o', c=c)
           for c in colors], [])

# Set axis names and legend
ax.set_xlabel(r"$X$ [kpc]", size=22)
ax.set_ylabel(r"$Y$ [kpc]", size=22)
ax.set_zlabel(r"$Z$ [kpc]", size=22)
ax.legend(cluster_name, loc='upper right', fontsize=16)

# Prepare the axes limits
ax.set_xlim((-95, 95))
ax.set_ylim((-95, 95))
ax.set_zlim((-95, 95))

# Set point-of-view: specified by (altitude degrees, azimuth degrees)
ax.view_init(30, 0)


# Initialization function: plot the background of each frame
def init():
    """Plot the background of each frame."""
    for line, pt in zip(lines, pts):
        line.set_data([], [])
        line.set_3d_properties([])

        pt.set_data([], [])
        pt.set_3d_properties([])
    return lines + pts


# Animation function. This will be called sequentially with the frame number
def animate(i):
    """We'll step twenty time-steps per frame. This leads to nice results."""
    i = (40*i) % x_t.shape[1]

    for line, pt, xi in zip(lines, pts, x_t):
        x, y, z = xi[:i].T
        line.set_data(x, y)
        line.set_3d_properties(z)

        pt.set_data(x[-1:], y[-1:])
        pt.set_3d_properties(z[-1:])

    ax.view_init(30, 0.3*(i/20))

    # set title name
    fig.suptitle(f'Time = -{i/10:.0f} [Myr]', fontsize=30)
    fig.canvas.draw()
    return lines + pts


# Instantiate the animator.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=2500, interval=1, blit=True)

# anim.save('final_result.gif')
writervideo = animation.FFMpegWriter(fps=20)
anim.save('Anim3D_orb_cart_case1.mov', writer=writervideo)

plt.show()
