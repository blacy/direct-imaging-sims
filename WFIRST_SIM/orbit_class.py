import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import brentq
# =====================================================================
# orbit object containing keplerian elements
# able to compute useful things, like phase angle and 
# planet-star separation for a given time
# taking into account geometry of the orbit
# NOTE this is missing some logic for edge cases
# =====================================================================
class Orbit:
    def __init__(self, param_array):

        # Keplerian orbital elements:
        ecc,inc,ome,tp,a,argperi = param_array[:6]	
        self.ecc = ecc # eccentricity
        self.inc = inc * np.pi/180.# inclination 
        self.ome = ome * np.pi/180.# longitude of ascending node
        self.tp = tp # time of pericenter passage
        self.a = a # semi-major axis in AU
        self.argperi = argperi * np.pi/180.# argument of periastron

        # planet-star system params:
        starmass, rp, dist = param_array[6:]
        self.starmass = starmass # stellar mass in solar masses
        self.rp = rp # radius of planet in jupiter radii
        self.dist = dist # distance to system in parsecs

        # pre-calculate a few things:
        self.per = self.a**1.5 /self.starmass # orbital period of planet around star in years


    def ta_to_alpha(self, ta):
        # given a true anomaly in radians
        # compute the phase angle also in radians
        cos_alpha = np.sin(ta + self.argperi)*np.sin(self.inc)*np.sin(self.ome) \
                    - np.cos(self.ome)*np.cos(ta+self.argperi)
        alpha = np.arccos(cos_alpha) # arccos maps to values 0 - pi
        return alpha

    def ta_to_dt(self, ta):
        # given a true anomaly 
        # compute the time since tp 
        dt = (self.per/(2.*np.pi)) * \
            (2.*np.arctan(np.sqrt((1.-self.ecc)/(1.+self.ecc))*np.tan(ta/2.)) \
            - (self.ecc*np.sin(ta)*np.sqrt(1.-self.ecc**2.))/(1. + self.ecc*np.cos(ta)))
        return dt

    def dt_ta_diff(self,ta,t):
        # function used in finding the 
        # true anomaly for a given time
        return t - self.ta_to_dt(ta)

    def dt_to_ta(self, t):
        # find ta such that: 0 = t - ta_to_t(ta)
        ta =  brentq(self.dt_ta_diff,-np.pi,np.pi,args=(t))
        return ta

    def ophase_to_dt(self, ophase): 
        # convert from orbital phase to time (in julian days)
        dt = ophase*self.per + self.tp
        return dt

    def t_to_ophase(self, t):
        # convert from a time (in julian days) to orbital phase
        return (t-self.tp) / self.per

    def ophase_to_sep(self, ophase):
        # convert from an orbital phase
        # to the instantaneous planet-star separation 
        # in units of au
        ta = self.dt_to_ta(self.ophase_to_dt(ophase))
        sep = self.a*(1.-self.ecc**2.)/(1.+self.ecc*np.cos(ta)) 
        return sep

    def ophase_to_alpha(self, ophase):
        # convert from an orbital phase
        # to the phase angle (in radians)
        ta = self.dt_to_ta(self.ophase_to_dt(ophase))
        alpha = self.ta_to_alpha(ta)
        return alpha

    def alpha_to_ta(self,alpha):
        # convert from phase angle to true anomaly
        # both in radians
        # note that arcsin maps to values -pi/2 to pi/2
        ta =  np.arcsin(np.cos(alpha)/np.sin(self.inc)) - self.argperi
        return ta 

    def alpha_to_ophase(self,alpha):
        # convert from phase angle in radians
        # to orbital phase (a fraction of the period)
        ophase = self.t_to_ophase(self.ta_to_dt(self.alpha_to_ta(alpha)))
        return ophase

    def alpha_to_sep(self, alpha):
        # given the phase angle in radians
        # return the instantaneous planet-star separation
        # in au
        sep = self.ophase_to_sep(self.alpha_to_ophase(alpha))
        return sep

