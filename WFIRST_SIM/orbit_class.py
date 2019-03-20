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
        self.per_days = self.per*365.0
        
    def ta_to_alpha(self, ta):
        # given a true anomaly in radians
        # compute the phase angle also in radians
        # cos_arg = np.sin(ta + self.argperi)*np.sin(self.inc) 
        # seem to need a shift of pi radians for this to behave as 
        # expected...
        cos_arg = np.sin(ta + self.argperi + np.pi)*np.sin(self.inc)
        alpha = np.arccos(cos_arg) # arccos maps to values 0 - pi
        return alpha 

    def ta_to_dt(self, ta):
        # given a true anomaly 
        # compute the time since most
        # recent periastron passage
        E = 2.*np.arctan(np.sqrt((1.-self.ecc)/(1.+self.ecc))*np.tan(ta/2.))
        dt = (self.per_days/(2.*np.pi))*( E - self.ecc*np.sin(E))
        if dt < 0:
            return dt + self.per_days
        return dt 

    def ta_to_sep(self,ta):
        # given a true anomaly, compute the instantaneous
        # separation in au
        sep = self.a*(1.-self.ecc**2.)/(1.+self.ecc*np.cos(ta))    
        return sep
    
    def dt_ta_diff(self,ta,t):
        # function used in finding the 
        # true anomaly for a given time
        return t - self.ta_to_dt(ta)

    def dt_to_ta(self, dt):
        # find ta such that: 0 = t - ta_to_dt(ta)
        # need to be careful with choosing
        # initial bounds for the root finder
        # so at this point
        # this function isn't vectorized... can do that with masking
        if dt < self.per_days:
            if dt == 0.0:
                return 0.0
            else:
                ta =  brentq(self.dt_ta_diff,0.0,np.pi*2.0,args=(dt))
                return ta 
        else:
            dt = self.jd_to_t(dt)
            if dt == 0.0:
                return 0.0
            else:
                ta =  brentq(self.dt_ta_diff,0.0,np.pi*2.0,args=(dt))
                return ta 
                        
    def jd_to_t(self,jd):
        # convert from jd to time past last
        # periastron passage 
        return (jd-self.tp)%self.per_days     
    
    def jdlist_to_alpha(self,jdlist):   
        # convert from a list of julian dates
        # to the phase angles the planet
        # will be at at those dates
        t_list = self.jd_to_t(jdlist)
        ta = [self.dt_to_ta(t) for t in t_list]
        alpha = self.ta_to_alpha(np.array(ta))   
        return alpha

    def jdlist_to_sep(self,jdlist):
        # get the instantaneous separation 
        # between star and planet at a given
        # julian date
        t_list = self.jd_to_t(jdlist)
        ta = [self.dt_to_ta(t) for t in t_list]
        sep = self.ta_to_sep(np.array(ta))   
        return sep

    def jdlist_to_alpha_sep(self,jdlist):   
        # convert from a list of julian dates
        # to the phase angles the planet
        # will be at on those dates
        t_list = self.jd_to_t(jdlist)
        ta = [self.dt_to_ta(t) for t in t_list]
        alpha = self.ta_to_alpha(np.array(ta))   
        sep = self.ta_to_sep(np.array(ta)) 
        return alpha, sep
    
    def ophase_to_dt(self, ophase): 
        # convert from orbital phase to time past
        # most recent periastron (in julian days)
        dt = ophase*self.per_days 
        return dt 
    
    def jd_to_ophase(self,jd):
        # convert from jd to orbital phase
        # as a fraction of the period
        return ((jd-self.tp)%self.per_days) / self.per_days

    def dt_to_ophase(self, dt):
        # convert from dt (a t less than a period past tp (in julian days) 
        # to orbital phase)
        return dt / self.per_days    
    
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
    
    def ophase_to_wa(self,ophase):
        pc_to_m = 30855152366503100.0  # convert parsecs to meters
        au_to_m = 149590000000.0   # convert AU to meters 149590000000
        alpha = self.ophase_to_alpha(ophase)
        sep = self.ophase_to_sep(ophase)
        wa = np.arcsin(np.sin(alpha)*sep*au_to_m/(self.dist*pc_to_m))
        return wa 