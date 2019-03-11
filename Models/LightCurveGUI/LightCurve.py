
import numpy as np
import matplotlib.pyplot as plt
from PhaseCurves import PhaseCurves
from Orbit import Orbit, OrbDist, FluxRatio
from AgFits import isofit, rayfit, vecfit
R_jup = 71492000.0 # m
au = 149597870700.0 # m

def make_lc(a1, omega_s, ecc, incl, argperi, Omega, Rp, a, models, delta):

    Rp = Rp*R_jup
    a = a*au

    MA = np.linspace(0,2*np.pi,361)
    Phi = PhaseCurves(omega_s)
    orbit = Orbit(ecc, incl, argperi, Omega)
    time, TA, alpha = orbit.get_arr(MA, ['t', 'TA', 'alpha'])
    dist = OrbDist(TA, a, ecc)

    iso_phi, ani_phi, ray_phi, vec_phi = Phi.phi_gui(alpha, delta)

    lw = 2

    if models[0]:
        lam_phi = Phi.lam(alpha)
        lam_flux = FluxRatio(lam_phi, dist, Phi.get_Ag('lam'), Rp, ecc)
        a1.plot(time, lam_flux, 'k:', label='Lambert', lw=lw)
    if models[1]:
        iso_flux = FluxRatio(iso_phi, dist, isofit(omega_s), Rp, ecc)
        a1.plot(time, iso_flux, 'g-', label='isotropic', lw=lw)
    if models[4]:
        ani_flux = FluxRatio(ani_phi, dist, Phi.get_Ag('ani'), Rp, ecc)
        a1.plot(time, ani_flux, 'r-', label='anisotropic', lw=lw)
    if models[2]:
        ray_flux = FluxRatio(ray_phi, dist, rayfit(omega_s), Rp, ecc)
        a1.plot(time, ray_flux, 'b--', label='scalar', lw=lw)
    if models[3]:
        vec_flux = FluxRatio(vec_phi, dist, vecfit(omega_s), Rp, ecc)
        a1.plot(time, vec_flux, 'b-', label='vector', lw=lw)

    a1.set_ylim(ymin=0)
    a1.minorticks_on()
    a1.set_xlabel(r'$(t - t_p)/P$')
    a1.set_ylabel(r'$\left(F_p/F_{\!\star}\right) \times\ 10^{-9}$')
    a1.legend(loc=0)
