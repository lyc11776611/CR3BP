import numpy as np
import matplotlib.pyplot as plt

from orbit_elements import *
from Settings import mu

import numpy.random as rd

# Initial condition and parameters

TU=100 # Time in total

d2r = np.pi/180.0
r2d = 1.0/d2r

omega = 0.1 # The omega in the perturbing term.

ts = np.linspace(0,TU,30000) # All time

# Keplarian initial orbital elements
a = 1
e = 0.5
i = 45
o = 90
O = 0
M = 0

# a = rd.rand()*5
# e = rd.rand()
# i = rd.rand()*180-90
# o = rd.rand()*360
# O = rd.rand()*360
# M = rd.rand()*360

orbelems = [[a,e,i,o,O,M]]
rv0 = get_r_v_from_orbelem(*orbelems[0])
rvecs = [rv0[0]]
vvecs = [rv0[1]]

orbelems[0][2] = 0 # Inclination is set to 0 for z axis set parallel to omega.

orbelems_dnlp = [aeioom_2_dnlp(orbelems[0],omega)]

for it in range(1,len(ts)):

    # aeioom = orbelems[-1]
    dt = ts[it] - ts[it-1]

    Lp,Gp,Hp,Mp,op,Op = orbelems_dnlp[-1]


    # Propagate
    Mp+= (1/Lp**3 * dt)*r2d

    op = np.mod(op,360)
    Op = np.mod(Op, 360)
    # Mp = np.mod(Mp, 360)

    dnlp = np.array([Lp,Gp,Hp,Mp,op,Op])
    dnl = dnlp_2_dnl(dnlp,omega)

    orbelems_dnlp.append(dnlp)
    orbelems.append(dnl_2_aeioom(dnlp_2_dnl(dnlp,omega)))

    # print(dnlp)

    aeioom = orbelems[-1]

    aeioom[2] = i # Put back the actual inclination during coordinate transform to cartesian
    rv = get_r_v_from_orbelem(*aeioom)

    # print(aeioom)

    rvecs.append(rv[0])
    vvecs.append(rv[1])

rvecs = np.array(rvecs)
vvecs = np.array(vvecs)


#### Plot stuffs.
ax = plt.figure().add_subplot(projection='3d')

# ax.plot(r_vec_s1[:, 0], r_vec_s1[:, 1], r_vec_s1[:, 2], label="Earth")
# ax.plot(r_vec_s2[:, 0], r_vec_s2[:, 1], r_vec_s2[:, 2], label="Object")
plt.plot(rvecs[:,0],rvecs[:,1],rvecs[:,2],color='grey')
print(rvecs[:,2])

ind0 = 0
for i in range(1,len(ts)):

    if rvecs[ind0,2]*rvecs[i,2]==0:
        ind0 = ind0 + 1

    if rvecs[ind0,2]*rvecs[i,2]<0:
        if rvecs[ind0,2]<=0:
            plt.plot(rvecs[ind0:i, 0], rvecs[ind0:i, 1], rvecs[ind0:i, 2], color='blue',alpha=0.3)
        else:
            plt.plot(rvecs[ind0:i, 0], rvecs[ind0:i, 1], rvecs[ind0:i, 2], color='blue',alpha=0.8)
        ind0 = i
if rvecs[ind0,2]<=0:
    plt.plot(rvecs[ind0:, 0], rvecs[ind0:, 1], rvecs[ind0:, 2], color='blue',alpha=0.3)
else:
    plt.plot(rvecs[ind0:, 0], rvecs[ind0:, 1], rvecs[ind0:, 2], color='blue',alpha=0.8)

# inds = np.where(rvecs[:,2]>=0)[0]
# plt.plot(rvecs[inds,0],rvecs[inds,1],rvecs[inds,2],color='red')
# inds = np.where(rvecs[:,2]<=0)[0]
# plt.plot(rvecs[inds,0],rvecs[inds,1],rvecs[inds,2],color='blue')


ax.scatter(0, 0, 0, c='red', marker="*", label="Primary")

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_zlim(-1.5, 1.5)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("$\omega = %.3f$  $\Delta t=%.1d$ Q1"%(omega,TU))


xs = np.linspace(-2, 2, 100)
zs = np.linspace(-2, 2, 100)

X, Z = np.meshgrid(xs, zs)

# ax.plot_surface(X, Z, np.zeros_like(X),alpha=0.1,color='black')

# plt.plot(rvecs[:,0],rvecs[:,1])
plt.show()