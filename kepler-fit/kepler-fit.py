################################################################################
# kepler-fit.py:
#   fit G2's astrometry and velocity to a Keplerian orbit.  used in
#   the appendix of McCourt & Madigan 2015
#
# Copyright 2015 Mike McCourt and Ann-Marie Madigan
#
# notes:
#
#   - this is written to be run over MPI on a computer cluster.
#     comments including "MPI" describe how to modify the code to run
#     on a single processor
#
# usage: mpirun python kepler-fit.py
#
#   - see stampede.pbs for an example batch script to run this on the
#     stampede supercomputer



################################################################################
# start by importing things we need
#
import math
from math import cos, sin, pi, tan, atan, sqrt, exp, log

import numpy as np
import emcee

# use brentq method to solve for eccentric anomaly
from scipy.optimize import brentq

# comment this out if not using MPI
from emcee.utils import MPIPool
import sys
# end MPI



################################################################################
# define constants and data.  use units of pc and years.  data from
# the wiki: https://wiki.mpe.mpg.de/gascloud/PlotsNData
#
# x = RA, y = dec, z increases along the line if sight
#
gm = 1.9393e-8


# radial velocity data, Br gamma
#
vtime = np.array([2003.2711, 2004.5356, 2005.2153, 2008.2624, 2009.3851,
                  2010.3539, 2011.3169, 2011.5665, 2012.2126, 2012.3612,
                  2012.4942, 2012.7036, 2013.2650])

vr = np.array([ 990.0, 1052.0, 1010.0, 1256.0, 1404.0, 1522.0, 1692.0,
               1640.0, 1934.0, 1930.0, 2030.0, 2046.0, 2180.0]) * 1.02269032e-6

vrerr = np.array([50.0, 30.0, 50.0, 20.0, 30.0, 25.0, 60.0, 60.0, 50.0, 50.0,
                  50.0, 50.0, 50.0]) * 1.02269032e-6


# astrometry data, Br gamma
#
time = np.array([2004.5356, 2005.2153, 2006.2040, 2008.2624, 2009.3851,
                 2010.3539, 2011.3169, 2011.5665, 2012.2126, 2012.4942,
                 2013.2650])

xpos = np.array([0.264754, 0.233448, 0.243039, 0.204608, 0.212171, 0.160914,
                 0.146314, 0.147089, 0.132222, 0.122072, 0.084641]) * 0.04

xerr = np.array([0.00247242, 0.00478122, 0.00489195, 0.00349524, 0.00338314,
                 0.00352225, 0.00268501, 0.00779815, 0.00289233, 0.00158594,
                 0.00151250]) * 0.04

ypos = np.array([-0.1455190, -0.1358090, -0.1239330, -0.0849295, -0.0812734,
                 -0.0510249, -0.0381918, -0.0275061, -0.0201099, -0.0223565,
                  0.0010505]) * 0.04

yerr = np.array([0.00264724, 0.00337244, 0.00441347, 0.00335883, 0.00265295,
                 0.00361687, 0.00266301, 0.00633855, 0.00307700, 0.00294338,
                 0.00160000]) * 0.04



################################################################################
# code to convert kepler elements to cartesian coordinates
#
# unit vector along angular momentum vector (normal to orbital plane)
def jhat(i, Omega, omega):
    return np.array([ sin(i) * sin(Omega),
                     -sin(i) * cos(Omega),
                      cos(i) ])

# unit vector along eccentricity vector (pointing towards periapse)
def ehat(i, Omega, omega):
    return np.array(
        [ cos(omega) * cos(Omega) - sin(omega) * sin(Omega) * cos(i),
          cos(omega) * sin(Omega) + sin(omega) * cos(Omega) * cos(i),
          sin(omega) * sin(i) ])

# the "other" vector (along the semi-minor axis)
def qhat(i, Omega, omega):
    return np.array(
        [ -sin(omega) * cos(Omega) - cos(omega) * sin(Omega) * cos(i),
          -sin(omega) * sin(Omega) + cos(omega) * cos(Omega) * cos(i),
           sin(i) * cos(omega) ])

# eccentric anomaly in terms of mean anomaly and eccentricity
# http://en.wikipedia.org/wiki/Mean_anomaly
def eccanomaly(m, e):
    # put m in the range [-pi, pi]
    mym = 2.0 * atan(tan(m / 2.0))

    # search for E in the range [-1.5 pi, 1.5 pi]
    return brentq(lambda E: E - e * sin(E) - mym,
                  -1.5*pi, 1.5*pi)



def cartesian_pos(a, e, ev, qv, meananomaly):
    ecc = eccanomaly(meananomaly, e)
    return a*(cos(ecc) - e) * ev + a*sqrt(1-e*e) * sin(ecc) * qv

def cartesian_vel(a, e, ev, qv, meananomaly):
    ecc = eccanomaly(meananomaly, e)

    prefact = 1.0/(1.0 - e*cos(ecc)) * sqrt(gm/a)
    term1 = sqrt(1.0-e*e) * cos(ecc) * qv
    term2 = e * sin(ecc) * ev

    return prefact * (term1 - term2)



################################################################################
# calculate a chi^2 against the data
#
#
def chisq_term(model, data, err, lnf):
    sigma2 = (err**2 + model**2 * exp(2*lnf))

    return 0.5 * ((model - data)**2 / sigma2 + log(2*pi*sigma2))


def chisq_pos(a, e, i, Omega, omega, tperi, lnf):
    # orientation vectors for the orbit
    p = ehat(i, Omega, omega)
    q = qhat(i, Omega, omega)

    # orbital phase (mean anomaly) at each observed time
    m = map(lambda t: sqrt(gm/(a**3)) * (t - tperi), time)

    pos = map(lambda mm: cartesian_pos(a, e, p, q, mm), m)
    pos = np.asarray(pos)
    xmodel = pos[:,0]
    ymodel = pos[:,1]

    errx = map(lambda a, b, c: chisq_term(a, b, c, lnf), xmodel, xpos, xerr)
    erry = map(lambda a, b, c: chisq_term(a, b, c, lnf), ymodel, ypos, yerr)

    return (sum(errx) + sum(erry))


def chisq_vel(a, e, i, Omega, omega, tperi, lnf):
    # orientation vectors for the orbit
    p = ehat(i, Omega, omega)
    q = qhat(i, Omega, omega)

    # orbital phase (mean anomaly) at each observed time
    m = map(lambda t: sqrt(gm/(a**3)) * (t - tperi), vtime)

    model = map(lambda mm: cartesian_vel(a, e, p, q, mm), m)
    model = np.asarray(model)
    model = model[:,2]

    err = map(lambda a, b, c: chisq_term(a, b, c, lnf), model, vr, vrerr)

    return (sum(err))



################################################################################
# initial conditions for the MCMC walkers
#
def random_ic():
    a     = np.random.normal(1.25,   0.0247)
    e     = np.random.normal(0.9762, 0.00074)

    i     = np.random.normal(61.9, 0.20)
    Omega = np.random.normal( 8.1, 0.43)
    omega = np.random.normal(97.2, 0.22)

    tperi = np.random.normal(2014.25, 0.006)

    lnfp = np.random.normal(-2.0, 0.5)
    lnfv = np.random.normal(-2.0, 0.5)

    return [a, e, i, Omega, omega, tperi, lnfp, lnfv]



################################################################################
# assume likelihood = exp(-chi^2)
#
def lnprob(p):
    a     = p[0] * 0.04
    e     = p[1]

    i     = p[2] * pi / 180.0
    Omega = p[3] * pi / 180.0
    omega = p[4] * pi / 180.0

    tperi = p[5]

    lnfp  = p[6]
    lnfv  = p[7]

    # apply parameter constraints here
    #
    if a<=0.0 or e>=1.0 or e<=0.0:
        return -np.inf

    if tperi<2013.0 or tperi>2015.5:
        return -np.inf

    if lnfp > 2.0 or lnfp < -10.0:
        return -np.inf

    if lnfv > 4.0 or lnfv < -10.0:
        return -np.inf

    chisqp = chisq_pos(a, e, i, Omega, omega, tperi, lnfp)
    chisqv = chisq_vel(a, e, i, Omega, omega, tperi, lnfv)

    return -1.0 * (chisqp + chisqv)



################################################################################
# set up and run the MCMC simulation
#
ndim     = 8
nwalkers = 4096


# comment this out if not using MPI
pool = MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)
# end MPI


# generate an initial set of walkers and burn them in
p0 = [random_ic() for i in range(nwalkers)]

print "starting burn-in..."

# if not using MPI, delete pool=pool in the line below
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[], pool=pool)
pos, prob, state = sampler.run_mcmc(p0, 200)

print("Mean acceptance fraction: {0:.3f}"
      .format(np.mean(sampler.acceptance_fraction)))

print("Mean autocorrelation time: {0:.5f}"
      .format(np.mean(sampler.acor)))

sampler.reset()

print "...finished burn-in"
print " "



################################################################################
# now, run the simulation

iter = 0
for result in sampler.sample(pos, iterations=1000, storechain=True, thin=100):
    print iter
    iter = iter+1

print sampler.flatchain.shape
print sampler.lnprobability.reshape(-1).shape

np.savetxt("chain-g2.dat",  sampler.flatchain)
np.savetxt("lnprob-g2.dat", sampler.lnprobability.reshape(-1))

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
