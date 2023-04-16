import numpy as np
from Settings import *
deg2rad = np.pi/180.0
rad2deg = 1/deg2rad

def get_orbital_element(r_vec,v_vec):

    r = np.sum(r_vec**2)**0.5
    v = np.sum(v_vec**2)**0.5

    vr = np.dot(r_vec,v_vec)/r

    h_vec = np.cross(r_vec,v_vec)
    h = np.sum(h_vec**2)**0.5

    i = np.arccos(h_vec[2]/h)

    N_vec = np.cross(np.array([0,0,1]),h_vec)
    N = np.sum(N_vec**2)**0.5

    if N_vec[1] >=0:
        Omega = np.arccos(N_vec[0]/N)
    else:
        Omega = 2*np.pi - np.arccos(N_vec[0]/N)

    e_vec = 1/mu*((v**2-mu/r)*r_vec-r*vr*v_vec)

    e = np.sum(e_vec**2)**0.5

    omega = np.arccos(np.dot(N_vec,e_vec)/N/e)
    if e_vec[2] <0:
        omega = 2*np.pi - omega

    theta = np.arccos(np.dot(e_vec,r_vec)/e/r)
    if vr<0:
        theta = 2*np.pi - theta

    if e<1: # Elliptical orbits
        a = h**2/mu/(1-e**2)
        return(a / au,
               i / np.pi*180.0,
               Omega/np.pi*180.0,
               e,
               omega/np.pi*180.0,
               theta/np.pi*180.0)
    if e>=1: # Hyperbolic trajectory
        q = h ** 2 / mu / (1 + e)

        # Get time of periapsis passage
        # all_thetas = np.linspace(0,theta,1000)
        # dtheta = all_thetas[1] - all_thetas[0]
        # t = np.sum(dtheta/(1+e*np.cos(all_thetas))**2)/mu**2*h**3

        t = 1/(e**2-1)*(e*np.sin(theta)/(1+e*np.cos(theta)))-1/(e**2-1)**(3/2)*np.log(((e+1)**0.5+(e-1)**0.5*np.tan(theta/2))/((e+1)**0.5-(e-1)**0.5*np.tan(theta/2)))
        t = -t/mu**2*h**3
        tp = reftime+t*np.timedelta64(1, "s")


        return (q / au,
                i / np.pi*180.0,
                Omega / np.pi * 180.0,
                e,
                omega / np.pi * 180.0,
                t/dayins)


# Convert from mean anomaly to true anomaly
def M2theta(M,e):

    # Check if M=0/360:
    M = np.mod(M,360)
    if M==0:
        return 0.0

    M = M*deg2rad

    # 1.  Solve the equation from M to E: E-esin(E)=M
    E = M
    # Newton's method
    Nmax = 1000
    for i in range(Nmax):

        E1 = E-(M-E+e*np.sin(E))/(-1+e*np.cos(E))
        if np.abs((E1-E)/E1)<1e-8:
            E = E1
            break
        else:
            E = E1

    if i == Nmax-1:
        raise BrokenPipeError("Not converge!@")

    # Convert E to theta
    theta = 2*np.arctan(np.tan(E/2.0)*((1+e)/(1-e))**0.5)

    return theta*rad2deg

# Convert from mean anomaly to Eccentric anomaly
def M2E(M,e):

    # Check if M=0/360:
    M = np.mod(M,360)
    if M==0:
        return 0.0

    M = M*deg2rad

    # 1.  Solve the equation from M to E: E-esin(E)=M
    E = M
    # Newton's method
    Nmax = 1000
    for i in range(Nmax):

        E1 = E-(M-E+e*np.sin(E))/(-1+e*np.cos(E))
        if np.abs((E1-E)/E1)<1e-8:
            E = E1
            break
        else:
            E = E1

    if i == Nmax-1:
        raise BrokenPipeError("Not converge!@")

    return E*rad2deg

# Convert true anomaly to mean anomaly
def theta2M(theta,e):

    theta = theta * deg2rad
    M = 2*np.arctan((((1-e)/(1+e))**0.5*np.tan(theta/2)))-e*(1-e**2)**0.5*np.sin(theta)/(1+e*np.cos(theta)) # Curtis 2020, Eq.3.6
    M = M*rad2deg

    return M

def get_r_v_from_orbelem(a,e,i,o,O,M):

    o = o * deg2rad
    O = O * deg2rad
    # M = M * deg2rad
    i = i * deg2rad

    # True anomaly
    theta = M2theta(M,e)*deg2rad
    # Radius
    r = a*(1-e**2)/(1+e*np.cos(theta))
    # Specific angular momentum
    h = (mu*a*(1-e**2))**0.5
    p = a*(1-e**2)
    # Positions:
    x = r * (np.cos(O)*np.cos(o+theta) - np.sin(O)*np.sin(o+theta)*np.cos(i))
    y = r * (np.sin(O) * np.cos(o + theta) + np.cos(O) * np.sin(o + theta) * np.cos(i))
    z = r * (np.sin(i)*np.sin(o+theta))
    # Velocities
    xdot = x*h*e/r/p*np.sin(theta)-h/r*(np.cos(O) * np.sin(o + theta)+np.sin(O) * np.cos(o + theta)* np.cos(i))
    ydot = y * h * e / r / p * np.sin(theta) - h / r * (
                np.sin(O) * np.sin(o + theta) - np.cos(O) * np.cos(o + theta) * np.cos(i))

    zdot = z * h * e / r / p * np.sin(theta) + h / r * (np.sin(i)*np.cos(o+theta))

    r_vec = np.array([x,y,z])
    v_vec = np.array([xdot,ydot,zdot])

    return [r_vec,v_vec]

####### Coordinate transform functions in this problem:  #######
# rv: non-rotating location and velocity.
# dnl (lhgmoo): Denaulay variable
# dnlp (lhgmoop): Denaulay variable after transformed to Kamiltonian by Lie-Hori theory
# aeioom: Keplarian orbital elements
# equ: Equinoctial element
# o (as a single variable function input): the omega of rotating frame, not the one in orbital element.

def rv_2_dnl(r_vec,v_vec):
    r = np.sum(r_vec**2)**0.5
    v = np.sum(v_vec**2)**0.5

    vr = np.dot(r_vec,v_vec)/r

    h_vec = np.cross(r_vec,v_vec)
    h = np.sum(h_vec**2)**0.5

    i = np.arccos(h_vec[2]/h)

    N_vec = np.cross(np.array([0,0,1]),h_vec)
    N = np.sum(N_vec**2)**0.5

    if N_vec[1] >=0:
        Omega = np.arccos(N_vec[0]/N)
    else:
        Omega = 2*np.pi - np.arccos(N_vec[0]/N)

    e_vec = 1/mu*((v**2-mu/r)*r_vec-r*vr*v_vec)

    e = np.sum(e_vec**2)**0.5

    omega = np.arccos(np.dot(N_vec,e_vec)/N/e)
    if e_vec[2] <0:
        omega = 2*np.pi - omega

    theta = np.arccos(np.dot(e_vec,r_vec)/e/r)
    if vr<0:
        theta = 2*np.pi - theta


    a = h**2/mu/(1-e**2)
    n = (mu / a ** 3) ** 0.5
    L = n*a**2
    G = L*(1-e**2)**0.5
    H = G*np.cos(i)

    M = theta2M(theta*rad2deg,e)

    return np.array([L,G,H,M,omega,Omega])

def dnl_2_dnlp(lhgmoo,o):
    L, G, H, M, omega, Omega = lhgmoo

    M *= deg2rad
    omega *= deg2rad
    Omega *= deg2rad

    Lp = L - o*H*L**3
    Gp = G
    Hp = H
    Mp = M - 3*o*H*M*L**2
    omegap = omega
    Omegap = Omega - o*M*L**3

    Mp *= rad2deg
    omegap *= rad2deg
    Omegap *= rad2deg
    # Mp = np.mod(Mp,360)
    omegap = np.mod(omegap,360)
    Omegap = np.mod(Omegap,360)

    return np.array([Lp, Gp, Hp, Mp, omegap, Omegap])

def dnlp_2_dnl(lhgmoop, o):
    Lp, Gp, Hp, Mp, omegap, Omegap = lhgmoop

    Mp *= deg2rad
    omegap *= deg2rad
    Omegap *= deg2rad

    L = Lp + o * Hp * Lp ** 3
    G = Gp
    H = Hp
    M = Mp + 3 * o * Hp * Mp * Lp ** 2
    omega = omegap
    Omega = Omegap + o * Mp * Lp ** 3

    M *= rad2deg
    omega *= rad2deg
    Omega *= rad2deg
    # M = np.mod(M,360)
    omega = np.mod(omega,360)
    Omega = np.mod(Omega,360)

    return np.array([L, G, H, M, omega, Omega])

def dnl_2_aeioom(lhgmoo):
    L, G, H, M, omega, Omega = lhgmoo

    i = np.arccos(H/G)*rad2deg
    e = (1-(G/L)**2)**0.5
    a = L**2/mu

    return [a,e,i,omega,Omega,M]

def dnl_2_rv(lhgmoo):

    return get_r_v_from_orbelem(*dnl_2_aeioom(lhgmoo))

def aeioom_2_dnl(aeioom):
    a, e, i, o, O, M = aeioom
    n = (mu/a**3)**0.5
    L = n * a ** 2

    G = L*(1-e**2)**0.5
    H = G*np.cos(i*deg2rad)

    return np.array([L,G,H,M,o,O])

def aeioom_2_dnlp(aeioom,o):

    return dnl_2_dnlp(aeioom_2_dnl(aeioom),o)

def aeioom_2_equ(aeioom):
    a, e, i, o, O, M = aeioom

    o *= deg2rad
    O *= deg2rad
    M *= deg2rad
    i *= deg2rad

    h = e * np.sin(o + O)
    k = e * np.cos(o + O)
    p = np.tan(i/2)*np.sin(O)
    q = np.tan(i/2)*np.cos(O)

    L = o+O+M
    L *= rad2deg

    return np.array([a,h,k,p,q,L])
