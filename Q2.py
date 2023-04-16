import numpy as np
import matplotlib.pyplot as plt

from orbit_elements import get_r_v_from_orbelem
from Settings import mu
import numpy.random as rd


# Initial condition and parameters

TU=100

d2r = np.pi/180.0
r2d = 1.0/d2r

omega = 0.05

ts = np.linspace(0,TU,30000)

# Keplarian initial orbital elements
a = 1
e = 0.5
i = 45
o = 90
O = 0
M = 0

# a = rd.rand()*3
# e = rd.rand()
# i = rd.rand()*180-90
# o = rd.rand()*360
# O = rd.rand()*360
# M = rd.rand()*360

# Setup system in initial conditions.
orbelems = [[a,e,i,o,O,M]]

rv0 = get_r_v_from_orbelem(*orbelems[0])

rvecs = [rv0[0]]
vvecs = [rv0[1]]

# Propagate
for it in range(1,len(ts)):

    a,e,i,o,O,M = orbelems[-1]
    dt = ts[it] - ts[it-1]

    n = (mu/a**3)**0.5
    L = n * a ** 2

    # Propagate
    M+= (1/L**3 * dt)*r2d
    O+= (omega*dt)*r2d

    o = np.mod(o,360)
    O = np.mod(O, 360)
    M = np.mod(M, 360)

    orbelems.append([a,e,i,o,O,M])

    rv = get_r_v_from_orbelem(*orbelems[-1])

    rvecs.append(rv[0])
    vvecs.append(rv[1])

rvecs = np.array(rvecs)
vvecs = np.array(vvecs)

print(orbelems)


# Plot trajectory
ax = plt.figure().add_subplot(projection='3d')

# ax.plot(r_vec_s1[:, 0], r_vec_s1[:, 1], r_vec_s1[:, 2], label="Earth")
# ax.plot(r_vec_s2[:, 0], r_vec_s2[:, 1], r_vec_s2[:, 2], label="Object")
plt.plot(rvecs[:,0],rvecs[:,1],rvecs[:,2],color='grey')
print(rvecs[:,2])


# Use different transparency in different part of trajectory
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
ax.set_title("$\omega = %.3f$  $\Delta t=%.1d$ Q2"%(omega,TU))


xs = np.linspace(-2, 2, 100)
zs = np.linspace(-2, 2, 100)

X, Z = np.meshgrid(xs, zs)

# ax.plot_surface(X, Z, np.zeros_like(X),alpha=0.1,color='black')

# plt.plot(rvecs[:,0],rvecs[:,1])
plt.show()










