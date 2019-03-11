import numpy as np

def OrbDist(TA, a, ecc):
    """
    Calculate orbital distance given true anamoly (TA)
    semi-major axis (a) and eccentricity (ecc).
    """
    return a*(1.0 - ecc*ecc)/(1.0 + ecc*np.cos(TA))

def FluxRatio(phi, dist, Ag, Rp, ecc, factor=1.0e9):
    return Ag*np.power(Rp/dist, 2)*phi*factor

class Orbit:
    """
    This class contains the necessary functions 
    to calculate the phase curves as a function 
    of time given arbitrary orbital parameters.
    """

    def __init__(self, ecc, incl, argperi, Omega):
        deg2rad = (np.pi/180.0)

        self.ecc = ecc
        self.incl = incl*deg2rad
        self.argperi = argperi*deg2rad
        self.Omega = Omega*deg2rad

    def get_EA(self, MA):
        """
        Get the eccentric anomaly (EA) 
        given the mean anomaly (MA) and 
        eccentricity (ecc).
        """
        EA = MA 
        count = 0
        delta = 1.0
        if MA%np.pi<1e-8:
            EA = MA
        else:
            while delta > 1.0e-10:
                EA_prev = EA
                EA = MA + self.ecc*np.sin(EA)
                delta = np.abs((EA - EA_prev)/EA_prev)
        return EA

    def get_TA(self, MA, EA):
        """
        Get the true anomaly (TA) given the eccentric
        anomaly (EA) and eccentricity (ecc).
        """
        TA = np.arccos((np.cos(EA) - self.ecc)/(1.0 - self.ecc*np.cos(EA)))
        if MA >= np.pi: TA = 2*np.pi - TA
        return TA

    def get_t(self, MA, TA):
        """
        Get t/P given the true anomaly (TA) and eccentricity (ecc).
        """
        term_1 = np.arctan2(np.tan(TA/2.0)*np.sqrt(1.0-self.ecc), np.sqrt(1.0+self.ecc))
        if MA > np.pi: term_1 = np.pi + term_1
        term_2 = (self.ecc*np.sin(TA)*np.sqrt(1.0-self.ecc*self.ecc))/(1+self.ecc*np.cos(TA))
        factor = 1.0/(2.0*np.pi)
        t = factor*(2.0*term_1 - term_2)
        return t

    def get_phase(self, TA):
        """
        Get phase angle.
        """
        f = TA + self.argperi
        phase = np.arccos(np.sin(f)*np.sin(self.incl)*np.sin(self.Omega) -\
                          np.cos(self.Omega)*np.cos(f))
        return phase

    def get_arr(self, MA_arr, arr=['alpha']):
        """
        Returns a numpy array of the indicated parameter(s), 
        where the elements correspond to the input MA values.  
        """
        alpha_arr = []
        EA_arr = []
        TA_arr = []
        t_arr = []
        return_arr = []
        for MA in MA_arr:
            EA = self.get_EA(MA)
            TA = self.get_TA(MA, EA)
            t = self.get_t(MA, TA)
            alpha = self.get_phase(TA)
            alpha_arr.append(alpha)
            EA_arr.append(EA)
            TA_arr.append(TA)
            t_arr.append(t)
        for arr_ in arr:
            if arr_=='alpha':
                return_arr.append(np.array(alpha_arr))
            elif arr_=='EA':
                return_arr.append(np.array(EA_arr))
            elif arr_=='TA':
                return_arr.append(np.array(TA_arr))
            elif arr_=='t':
                return_arr.append(np.array(t_arr))
        return return_arr if len(arr)>1 else return_arr[0]
