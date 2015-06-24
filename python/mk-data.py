################################################################################
# drag-force-mla.py:
#   runs the MCMC simulation to make figure 6 in McCourt & Madigan 2015.
#
# Copyright 2015 Mike McCourt and Ann-Marie Madigan
#
# notes:
#
#   - this program consists of two parts: there's a physical model
#     which integrates trajectories for G2 using the equation of
#     motion described in the paper.  this part is done by the
#     function integrate_model.  next, there's a likelihood
#     calculation (i.e., a "chi^2") determined by comparing with the
#     observational data from Gillessen et al. (2013b) and from Pfuhl
#     et al. (2015).  the likelihood is really all we need... we feed
#     that into emcee and let it sample the parameter space.
#
#   - this is written to be run over MPI on a computer cluster.
#     comments including "MPI" describe how to modify the code to run
#     on a single processor
#
# usage: mpirun python drag-force-mla.py
#
#   - see stampede.pbs for an example batch script to run this on the
#     stampede supercomputer



################################################################################
# start by importing some things we need
#
from math import acos, cos, pi, sin, sqrt, exp, log

import numpy as np
from numpy import array, cross, dot, zeros, sign

from scipy.integrate import odeint


################################################################################
# define global constants (units of pc and years)
#
gm = 1.9393e-8



################################################################################
# math functions which should be in python, but aren't
#
def norm(vec):
    return sqrt(dot(vec,vec))

# rotation matrix which translates one vector to lie along another
# also returns the inverse
#
# from http://svn.gna.org/svn/relax/tags/1.3.4/maths_fns/rotation_matrix.py
#
def rotation_matrix(vector_orig, vector_fin):
    # normalize to unit vectors
    vector_orig = vector_orig / norm(vector_orig)
    vector_fin  = vector_fin  / norm(vector_fin)

    # unit vector along the rotation axis
    axis = cross(vector_orig, vector_fin)
    axis_len = norm(axis)
    if axis_len != 0.0:
        axis = axis / axis_len

    # alias the axis coordinates
    x = axis[0]
    y = axis[1]
    z = axis[2]

    # the rotation angle
    angle = acos(dot(vector_orig, vector_fin))

    # trig functions (only need to do this maths once)
    ca = cos(angle)
    sa = sin(angle)

    # calculate the rotation matrix elements.
    R = zeros((3,3))
    R[0,0] =   1.0 + (1.0-ca) * (x**2-1.0)
    R[0,1] = -z*sa + (1.0-ca) * x*y
    R[0,2] =  y*sa + (1.0-ca) * x*z
    R[1,0] =  z*sa + (1.0-ca) * x*y
    R[1,1] =   1.0 + (1.0-ca) * (y**2-1.0)
    R[1,2] = -x*sa + (1.0-ca) * y*z
    R[2,0] = -y*sa + (1.0-ca) * x*z
    R[2,1] =  x*sa + (1.0-ca) * y*z
    R[2,2] =   1.0 + (1.0-ca) * (z**2-1.0)

    # now calculate the inverse matrix
    sa = -sa

    # calculate the rotation matrix elements.
    Rinv = zeros((3,3))
    Rinv[0,0] =   1.0 + (1.0-ca) * (x**2-1.0)
    Rinv[0,1] = -z*sa + (1.0-ca) * x*y
    Rinv[0,2] =  y*sa + (1.0-ca) * x*z
    Rinv[1,0] =  z*sa + (1.0-ca) * x*y
    Rinv[1,1] =   1.0 + (1.0-ca) * (y**2-1.0)
    Rinv[1,2] = -x*sa + (1.0-ca) * y*z
    Rinv[2,0] = -y*sa + (1.0-ca) * x*z
    Rinv[2,1] =  x*sa + (1.0-ca) * y*z
    Rinv[2,2] =   1.0 + (1.0-ca) * (z**2-1.0)

    return (R, Rinv)



################################################################################
# integrate the equation of motion
#
# helper function to compute the drag force
#
def drag_force(r,v,  alpha, logbeta, fkep, theta, phi, loggf):
    gf = exp(loggf)             # gf is L/R (= "geometric factor")
    Sigma = 1015691.12          # 3 earth masses / (100 AU)^2

    # get the relative velocity
    #
    J = array([sin(theta) * cos(phi),
               sin(theta) * sin(phi),
               cos(theta)])

    vbg  = (cross(J,r)/norm(J)/norm(r)) * fkep * sqrt(gm/norm(r))
    vrel = v - vbg

    # get the prefactor for the drag force
    #
    rho     = (1.65e5) * (norm(r)/0.04)**(-alpha)
    machsq  = (norm(v - vbg))**2 / ((5.0/3.0)*gm/norm(r))
    prefact = rho/Sigma * (1.0 + 2.0/(exp(logbeta)*machsq))

    # rotate to the coordinate system aligned with the cloud, find the
    # drag force in that frame, and rotate back to the original
    # coordinate system
    #
    (c, cinv) = rotation_matrix(v, array([1.0, 0.0, 0.0]))
    tmp = dot(c, vrel)

    area = array([[1.0,    0.0,    0.0],
                  [0.0, 1.0+gf,    0.0],
                  [0.0,    0.0, 1.0+gf]])

    ram = sign(tmp) * dot(area, tmp**2);
    return dot(cinv, ram) * prefact

# derivative vector
# - Y is the state vector [r, v]
# - return the derivative vector [v, a]
#
def f(Y,t,  alpha, logbeta, fkep, theta, phi, loggf):
    x, y, z, vx, vy, vz = Y
    d3 = (x**2 + y**2 + z**2)**1.5

    fx, fy, fz = drag_force(array([ x,  y,  z]),
                            array([vx, vy, vz]),
                            alpha, logbeta, fkep, theta, phi, loggf)

    return np.array([ vx, vy, vz,
                      -gm*x/d3 - fx,
                      -gm*y/d3 - fy,
                      -gm*z/d3 - fz])

# integrate the trajectory
#
def integrate_model(alpha, logbeta, fkep, theta, phi, loggf,
                    raShift, decShift, vraShift, vdecShift):
    # initial condition
    #
    t0   = 2013.33
    tend = 2040.0
    tbeg = 1990.0

    y0 = np.array([ 0.00353743   + raShift,
                   -0.0000267992 + decShift,
                   -0.00108785,
                   #
                   -0.002075    + vraShift,
                    0.000850869 + vdecShift,
                    0.00218925])


    # times to store the solution
    #
    dt   = 5.0e-2
    tfwd = np.arange(t0,  tend,  dt)
    tbak = np.arange(t0,  tbeg, -dt)


    # integrate forwards and backwards
    #
    ofwd = odeint(f, y0, tfwd, args=(alpha, logbeta, fkep, theta, phi, loggf))
    obak = odeint(f, y0, tbak, args=(alpha, logbeta, fkep, theta, phi, loggf))


    # join the solutions
    #
    sol = np.concatenate((np.flipud(obak), ofwd))
    t   = np.concatenate((np.flipud(tbak), tfwd))

    return (t, sol)



alpha     =  0.6
logbeta   =  0.3
fkep      =  0.5
theta     =  1.754460524831643609e+00
phi       = -6.648098517051298506e-01
loggf     =  log(14.0)
raShift   =  1.114870617771124360e-04
decShift  =  1.661446928084823094e-05
vraShift  =  9.256909974821245766e-05
vdecShift =  3.162990726641473789e-05

t, sol = integrate_model(alpha, logbeta, fkep, theta, phi, loggf,
                         raShift, decShift, vraShift, vdecShift)

x  = sol[:,0]
y  = sol[:,1]
vz = sol[:,5]


# interpolate solutions to the times matching telescope observations
#
g2vtime = array([2003.2711, 2004.5356, 2005.2153, 2008.2624, 2009.3851,
                 2010.3539, 2011.3169, 2011.5665, 2012.2126, 2012.3612,
                 2012.4942, 2012.7036, 2013.2650])

g2vr = np.interp(g2vtime, t, vz)


g1vtime = array([2004.50, 2006.21, 2008.25]) + 13.0

g1vr = np.interp(g1vtime, t, vz)



g2time = array([2004.5356, 2005.2153, 2006.2040, 2008.2624, 2009.3851,
                2010.3539, 2011.3169, 2011.5665, 2012.2126, 2012.4942,
                2013.2650])

g2xpos = np.interp(g2time, t, x)
g2ypos = np.interp(g2time, t, y)

g1time = array([2003.36, 2004.50, 2005.36, 2006.43, 2007.25, 2009.21, 2010.50,
                #
                2004.54, 2006.21, 2008.26]) + 13.0

g1xpos = np.interp(g1time, t, x)
g1ypos = np.interp(g1time, t, y)



print "g2vr"
print g2vr

print "g1vr"
print g1vr

print "g2xpos"
print g2xpos

print "g2ypos"
print g2ypos

print "g1xpos"
print g1xpos

print "g1ypos"
print g1ypos




# fin

# Local Variables:
# coding: utf-8-unix
# End:
