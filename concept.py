import numpy as np
from typing import Callable

"""
This script intends to demonstrate a comprehensive understanding on the
computational approximation of differential equations through the
equations of motion in the circular restricted three body problem.
Specifically, this is a proof of concept of the computational construct
of the model created by our equations of motion.

The script is annotated to describe each step.
"""

# Gravitational constant
G = 6.67408e-11

# Mass of the Sun
m1 = 1.9885e30 # kg

# Mass of the Earth
m2 = 5.972e24 # kg

# Distance of the Earth to the Sun
R = 1.49598e11 # m

# Radial velocity of the Earth
radial_velocity = np.sqrt(G*(m1+m2)/R**3)
        

def LN(r):
    """
    Simplified accleration in the x-component to find the Lagrange
    Points.
    """
    f1 = -G*m1/((r)*abs(r))
    f2 = -G*m2/((r-R)*abs(r-R))
    ang_v = G*r*(m1+m2)/(R**3)
    return ang_v + f1 + f2


def get_LagrangePoints():
    """
    Finding the Lagrange points will require a way to computationally

    """
    ...


def r1(x,y,z):
    """The position vector of the object relative to the Sun."""
    return np.sqrt(x**2 + y**2 + z**2)


def r2(x,y,z):
    """The position vector of the object relative to the Earth."""
    return np.sqrt((x-R)**2 + y**2 + z**2)


def xddot(t, x, y, z, ydot):
    """Acceleration in the x-component."""
    R1 = r1(x,y,z)
    R2 = r2(x,y,z)
    xddot = (x*radial_velocity**2 + 2*radial_velocity*ydot - G*m1*x/(R1**3)
        - G*m2*(x-R)/(R2**3))
    return xddot


def yddot(t, x, y, z, xdot):
    """Acceleration in the y-component."""
    R1 = r1(x,y,z)
    R2 = r2(x,y,z)
    yddot = (y*radial_velocity**2 - 2*radial_velocity*xdot - G*m1*y/(R1**3)
        - G*m2*y/(R2**3))
    return yddot


def zddot(t, x, y, z):
    """Acceleration in the z-component."""
    R1 = r1(x,y,z)
    R2 = r2(x,y,z)
    zddot = - G*m1*z/(R1**3) - G*m2*z/(R2**3)
    return zddot


def POSITION(x_0, y_0, z_0, xdot_0, ydot_0, zdot_0, t, h=1e-6):
    ...


def main():
    ...

if __name__ == '__main__':
    main()