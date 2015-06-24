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
from scipy.optimize  import minimize_scalar

import emcee

# comment this out if not using MPI
from emcee.utils import MPIPool
import sys
# end MPI



################################################################################
# define global constants (units of pc and years)
#
gm = 1.9393e-8



################################################################################
# G2 data from the wiki: https://wiki.mpe.mpg.de/gascloud/PlotsNData
#
# x = RA, y = dec, z increases along the line of sight
#

# radial velocity data, Br gamma
#
g2vtime = array([2003.2711, 2004.5356, 2005.2153, 2008.2624, 2009.3851,
                 2010.3539, 2011.3169, 2011.5665, 2012.2126, 2012.3612,
                 2012.4942, 2012.7036, 2013.2650])

g2vr = array([ 990.0, 1052.0, 1010.0, 1256.0, 1404.0, 1522.0, 1692.0,
              1640.0, 1934.0, 1930.0, 2030.0, 2046.0, 2180.0]) * 1.02269032e-6

g2vrerr = array([50.0, 30.0, 50.0, 20.0, 30.0, 25.0, 60.0, 60.0, 50.0,
                 50.0, 50.0, 50.0, 50.0]) * 1.02269032e-6


# radial velocity for G1, Br gamma.  from Pfuhl et al. (2015)
#
g1vtime = array([2004.50, 2006.21, 2008.25])

g1vr = array([-2051.526718, -1593.511450, -1135.496183]) * 1.02269032e-6

g1vrerr = array([143.129771, 143.129771, 143.129771]) * 1.02269032e-6



# astrometry data for G2, Br gamma.  from the wiki
#
g2time = array([2004.5356, 2005.2153, 2006.2040, 2008.2624, 2009.3851,
                2010.3539, 2011.3169, 2011.5665, 2012.2126, 2012.4942,
                2013.2650])

g2xpos = array([0.264754, 0.233448, 0.243039, 0.204608, 0.212171, 0.160914,
                0.146314, 0.147089, 0.132222, 0.122072, 0.084641]) * 0.04

g2xerr = array([0.00247242, 0.00478122, 0.00489195, 0.00349524, 0.00338314,
                0.00352225, 0.00268501, 0.00779815, 0.00289233, 0.00158594,
                0.00151250]) * 0.04

g2ypos = array([-0.1455190, -0.1358090, -0.1239330, -0.0849295, -0.0812734,
                -0.0510249, -0.0381918, -0.0275061, -0.0201099, -0.0223565,
                 0.0010505]) * 0.04

g2yerr = array([0.00264724, 0.00337244, 0.00441347, 0.00335883, 0.00265295,
                0.00361687, 0.00266301, 0.00633855, 0.00307700, 0.00294338,
                0.00160000]) * 0.04


# astrometry data for G1, from Pfuhl et al. (2015)
# first 7 points are L-band; last 3 are Br-gamma
#
g1time = array([2003.36, 2004.50, 2005.36, 2006.43, 2007.25, 2009.21, 2010.50,
                #
                2004.54, 2006.21, 2008.26])

g1xpos = array([-0.085751, -0.078626, -0.080662, -0.077354, -0.068448,
                -0.044529, -0.030025,
                #
                -0.108052, -0.072468, -0.039740]) * 0.04

g1xerr = array([0.005089, 0.004071, 0.003053, 0.003308, 0.003053, 0.003053,
                0.003308,
                #
                0.007013, 0.008312, 0.007792]) * 0.04

g1ypos = array([-0.063436, -0.071806, -0.081498, -0.096256, -0.109692,
                -0.120264, -0.122687,
                #
                -0.058559, -0.109459, -0.126577]) * 0.04

g1yerr = array([0.004846, 0.003965, 0.003084, 0.002863, 0.003084, 0.003084,
                0.003304,
                #
                0.007658, 0.007207, 0.007207]) * 0.04



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
    vbg  = vbg + r/norm(r) * sqrt(gm/norm(r))
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



################################################################################
# chi^2 error calculation
#
# chi^2 for a single variable (ra, dec, vr)
#
def chi_sq_term(t, sol, obstime, obsdata, obserr, lnf):
    # evaluate the model at the times of the observations
    model = np.interp(obstime, t, sol)

    # "extra" error s = f * rms(data)
    syssq = np.exp(2*lnf) * np.mean(obsdata**2)

    # assume P = (2 pi s^2)^(-1/2) * exp(-delta^2 / (2 s^2))
    #
    # analogue of chi^2 here is -ln P
    #
    err = map(lambda a, b, c:
                  0.5*((a-b)**2/(c**2+syssq) + np.log(2*pi*(c**2+syssq))),
              model, obsdata, obserr)

    return sum(err)


# chi^2 for the g2 data
#
def chi_sq_g2(t, sol, lnfp, lnfv):
    x  = sol[:,0]
    y  = sol[:,1]
    vz = sol[:,5]

    xerr = chi_sq_term(t,  x, g2time,  g2xpos, g2xerr,  lnfp)
    yerr = chi_sq_term(t,  y, g2time,  g2ypos, g2yerr,  lnfp)
    verr = chi_sq_term(t, vz, g2vtime, g2vr,   g2vrerr, lnfv)

    return (xerr + yerr + verr)


# chi^2 for the g1 data. (requires minimizing over dtg1, the time-lag
# between G2 and G1.)
#
def chi_sq_g1_intermediate(dtg1, t, sol, lnfp, lnfv):
    x  = sol[:,0]
    y  = sol[:,1]
    vz = sol[:,5]

    xerr = chi_sq_term(t,  x, g1time  + dtg1,  g1xpos, g1xerr,  lnfp)
    yerr = chi_sq_term(t,  y, g1time  + dtg1,  g1ypos, g1yerr,  lnfp)
    verr = chi_sq_term(t, vz, g1vtime + dtg1,  g1vr,   g1vrerr, lnfv)

    return (xerr + yerr + verr)

def chi_sq_g1(t, sol, lnfp, lnfv):
    res = minimize_scalar(chi_sq_g1_intermediate,
                          [0.0, 20.0],
                          args=(t, sol, lnfp, lnfv))

    return (res.x, res.fun)


# total chi^2 for the model
#
def chi_sq(alpha, logbeta, fkep, theta, phi, loggf,
           raShift, decShift, vraShift, vdecShift, lnfp, lnfv):

    t, sol = integrate_model(alpha, logbeta, fkep, theta, phi, loggf,
                             raShift, decShift, vraShift, vdecShift)

    dtg1, chisqg1 = chi_sq_g1(t, sol, lnfp, lnfv)

    return (chi_sq_g2(t, sol, lnfp, lnfv) + chisqg1)



################################################################################
# now, we've done our part.  all that remains to be done is to feed
# the chi_sq function to the emcee optimizer, let it explore the
# parameter space, and keep track of the steps it takes



################################################################################
# MCMC optimization
#
# initial conditions for the MCMC walkers.  follow the advice from the
# emcee documentation and initialize walkers in a narrow ball around
# the best-fit solution.
#
# I used Mathematica's simulated annealing optimizer to find this
# best-fit solution.  that also indicated that the probability
# distribution isn't multi-modal, so it's ok to start the walkers off
# around just a single point.
#
def random_ic():
    # start off at the best-fit solution
    alpha     =  8.598847643691621689e-01
    logbeta   =  3.199742509393102008e-01
    fkep      =  6.851834170728353657e-01
    theta     =  1.754460524831643609e+00
    phi       = -6.648098517051298506e-01
    loggf     =  2.695946924174811077e+00
    raShift   =  1.114870617771124360e-04
    decShift  =  1.661446928084823094e-05
    vraShift  =  9.256909974821245766e-05
    vdecShift =  3.162990726641473789e-05
    lnfp      = -3.633534064585596912e+00
    lnfv      = -7.304906199422774193e+00

    # seed them with a small, ~1% spread
    alpha     = np.random.normal(alpha,     np.abs(0.01 * alpha))
    logbeta   = np.random.normal(logbeta,   np.abs(0.01 * logbeta))
    fkep      = np.random.normal(fkep,      np.abs(0.01 * fkep))
    theta     = np.random.normal(theta,     np.abs(0.01 * theta))
    phi       = np.random.normal(phi,       np.abs(0.01 * phi))
    loggf     = np.random.normal(loggf,     np.abs(0.01 * loggf))
    raShift   = np.random.normal(raShift,   np.abs(0.01 * raShift))
    decShift  = np.random.normal(decShift,  np.abs(0.01 * decShift))
    vraShift  = np.random.normal(vraShift,  np.abs(0.01 * vraShift))
    vdecShift = np.random.normal(vdecShift, np.abs(0.01 * vdecShift))
    lnfp      = np.random.normal(lnfp,      np.abs(0.01 * lnfp))
    lnfv      = np.random.normal(lnfv,      np.abs(0.01 * lnfv))

    return [alpha, logbeta, fkep, theta, phi, loggf,
            raShift, decShift, vraShift, vdecShift, lnfp, lnfv]


# assume likelihood = exp(-chi^2)
#
def lnprob(p):
    alpha     = p[0]
    logbeta   = p[1]
    fkep      = p[2]
    theta     = p[3]
    phi       = p[4]
    loggf     = p[5]
    raShift   = p[6]
    decShift  = p[7]
    vraShift  = p[8]
    vdecShift = p[9]
    lnfp      = p[10]
    lnfv      = p[11]

    # apply parameter constraints here.  priors should go here, as well.
    #
    if alpha<=0.0 or alpha>1.2:
        return -np.inf

    if logbeta<=np.log(1.0e-2) or logbeta>np.log(1.0e2):
        return -np.inf

    if fkep<=0.0 or fkep>1.0:
        return -np.inf

    if theta<0.0 or theta>pi:
        return -np.inf

    if phi<-pi or phi>pi:
        return -np.inf

    if loggf<np.log(1e-10) or loggf>np.log(1.0e2):
        return -np.inf

    if lnfp < -10.0 or lnfp > 2.0:
        return -np.inf

    if lnfv < -10.0 or lnfv > 2.0:
        return -np.inf

    # constraints on rashift, etc?

    # integration sometimes fails if the cloud spirals into the black
    # hole.  I couldn't find a way to tell the integrator to stop, so
    # we have to let it fail and return a large chi^2 when it does.

    try:
        ret = chi_sq(alpha, logbeta, fkep, theta, phi, loggf,
                     raShift, decShift, vraShift, vdecShift, lnfp, lnfv)
    except:
        ret = np.inf

    return -1.0 * ret



################################################################################
# set up and run the MCMC simulation
#
# use 4096 walkers to explore the 12-dimensional parameter space
#
ndim     = 12
nwalkers = 4096


# comment this out if not using MPI
pool = MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)
# end MPI



print "starting burn-in..."

# empirically, we need a very long burn-in phase.  something in the
# ballpark of 1000 steps!
#
# if not using MPI, remove pool=pool from the next line.  can replace
# it with threads=N to run on N threads
#
p0 = [random_ic() for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[], pool=pool)
pos, prob, state = sampler.run_mcmc(p0, 1000)

sampler.reset()

print "...finished burn-in"
print " "



# now, run the simulation
print "starting simulation..."

iter = 0
for result in sampler.sample(pos, iterations=10000, thin=100, storechain=True):
    print iter

    # write data to the chain file as we go
    if iter % 100 == 0:
        position = result[0]
        f = open("chain-mla.dat", "a")

        for k in range(position.shape[0]):
            tab = '   '.join([str(i) for i in position[k]])
            f.write("{0:4d} {1:s}\n".format(k, tab))

        f.close()

    # occasionally write diagnostic information about the random walks
    if iter % 500 == 0:
        print("Mean acceptance fraction: {0:.3f}"
              .format(np.mean(sampler.acceptance_fraction)))

        print("Mean autocorrelation time: {0:.5f}"
              .format(np.mean(sampler.acor)))


    iter = iter+1


# write out probabilities at the end, so we can find the best fit
#
np.savetxt("lnprob-mla-thin.dat", sampler.lnprobability.reshape(-1))

# finally, print out some diagnostics and quit
#
print("Mean acceptance fraction: {0:.3f}"
      .format(np.mean(sampler.acceptance_fraction)))

print("Mean autocorrelation time: {0:.5f}"
      .format(np.mean(sampler.acor)))

print "...finished simulation"

# comment out if not using MPI
pool.close()
# end MPI

# fin

# Local Variables:
# coding: utf-8-unix
# End:
