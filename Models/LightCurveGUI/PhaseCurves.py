#!/usr/bin/env python 

import numpy as np

class PhaseCurves:

    def __init__(self, omega_s):
        self.omega_s = omega_s
        self.Agdict = {'lam':(2./3.)*omega_s, 'iso':0.0, 'ani':0.0, 'ray':0.0, 'vec':0.0}

    def get_Ag(self, model):
        return self.Agdict[model]
       
    def lam(self, alpha):
        """
        Analytic Lambert phase function
        """
        phi = (np.sin(alpha) + (np.pi - alpha)*np.cos(alpha))/np.pi
        return phi

    def phi_gui(self, orb_alpha, delta):
        import os
        from GlobalVars import phaseDIR, colnames, ani_Ag_cols
        from scipy.interpolate import interp1d

        files = [dir for dir in os.listdir(phaseDIR) if dir[:3]=='phi']

        omega_arr = np.array([float(f[4:]) for f in files])
        arg = np.argmin(np.abs(omega_arr - self.omega_s))

        Phi = np.loadtxt(phaseDIR+files[arg], dtype=colnames, skiprows=1)

        phi_iso = interp1d(Phi['alpha'], Phi['phi_iso'], kind='cubic')
        phi_ani = interp1d(Phi['alpha'], Phi['g0='+str(delta)], kind='cubic')
        phi_ray = interp1d(Phi['alpha'], Phi['phi_ray'], kind='cubic')
        phi_vec = interp1d(Phi['alpha'], Phi['phi_vec'], kind='cubic')

        ani_Ag = np.loadtxt(phaseDIR+'anisotropic_Ag', dtype=ani_Ag_cols, skiprows=1)

        self.Agdict['ani'] = ani_Ag[(np.round(ani_Ag['omega'],2)==round(self.omega_s,2))\
                                  & (np.round(ani_Ag['g0'],1)==round(delta,1))]['Ag']

        Phi_iso = phi_iso(orb_alpha)
        Phi_ani = phi_ani(orb_alpha)
        Phi_ray = phi_ray(orb_alpha)
        Phi_vec = phi_vec(orb_alpha)

        return Phi_iso, Phi_ani, Phi_ray, Phi_vec
