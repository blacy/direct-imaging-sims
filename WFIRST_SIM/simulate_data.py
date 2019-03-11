import matplotlib.pyplot as plt
import numpy as np
from orbit_class import Orbit
import misc_utils as mu   
import ifs_noise_model as snr
import imaging_noise_model as snr_im
import juneper_model as hpm  
import coolTLUSTY_model as cpm
     
def generate_noisey_spectra(TEXP, bandcenter, width, params):
    wavelengths = mu.create_wl_range(bandcenter*(1.0 - width/200.0), bandcenter*(1.0 + width/200.0),params['R']) 
    noiseless_signal, noisy_signal, noise_level = snr.generate_noisey_data(TEXP, wavelengths, params)
    noisy_signal[np.where(noisy_signal<0.0)] = 0.0
    return wavelengths, noiseless_signal, noisy_signal, noise_level

def generate_noisey_imaging(TEXP, bandcenter, width, params):
    wavelengths = mu.create_wl_range(bandcenter*(1.0 - width/200.0), bandcenter*(1.0 + width/200.0),params['R']) 
    noiseless_signal, noisy_signal, noise_level = snr_im.generate_noisey_image_data(TEXP, wavelengths, params)
    if noisy_signal <0.0:
        noisy_signal = 0.0
    return bandcenter, noiseless_signal, noisy_signal, noise_level

def change_alpha(orbit, params, newalpha, pfunc):
    # NEWALPHA MUST BE IN RADIANS
    sep = orbit.alpha_to_sep(newalpha)
    phi = pfunc(newalpha)
    params.update({'phi': phi})
    params.update({'sep':sep})
    wa = mu.calc_wa(params['sep'],newalpha*180.0/np.pi,params['d_obs'])
    params.update({'wa':wa})

def create_header_string(TEXP, params):
    header = ''
    if 'met' in params.keys():
        for key in ['sep', 'phi', 'rp','met']:
            header = header +'%s = %f\n'%(key,params[key])
    else:
        for key in ['sep', 'phi', 'rp','meth','chromo']:
            header = header +'%s = %f\n'%(key,params[key])
    header = header + 'TEXP = %f\n'%TEXP 
    header = header + ' 1 : wavelength (microns)\n 2 : clean planet-star flux ratio\n 3 : noisy planet-star flux ratio\n 4 : error'
    return header


if __name__ == '__main__':

    # Example of setting up a star+planet system 
    # These are semi-arbitrary choices 
    # or you could make them consistent with a known RV planet
    partag = 'Fiducial'            # name your planet
    params = {}                    # create a dictionary for your simulation's parameters
    params.update({'d_obs':10.0})  # distance to the observer in parsecs
    params.update({'Mstar':5.0})   # absolute stellar V-band magnitude
    params.update({'Nez':1.0})     # exozodi level
    stellartype = 'g0v'            # stellar type
    stellar_mass = 1.0             # stellar mass in units of solar mass
    stellar_spec = 'AuxiliaryData/'+stellartype+'.dat' 
    ref_wl, ref_flambda = np.loadtxt(stellar_spec,unpack=True,usecols=(0,1))
    fstar_lambda = interp1d(ref_wl,ref_flambda)  # specific flux density W / m^2 / micron, for zero mag star
    params.update({'fstar_lambda':fstar_lambda}) # a function which returns 
                                                 # specific flux density for any wavelength
    params.update({'rp':1.0})  # planet radius in Jupiter radii
    params.update({'chromo':0.75,'meth':0.25}) # chromophore and methane level 
    params.update({'Ag':hpm.juneper_Agfunc(params['chromo'],params['meth'])}) 
    a = 3.8           # semimajor axis (in au)
    ecc = 0.000000001  # eccentricity
    inc = 90.0         # inclination (degrees)
    ome = 90.0         # longitude of ascending node (degrees)
    tp = 0.0           # epoch of perihelion passage (julian date)
    argperi = 90.0     # argument of perihelion (degrees)
    orbit_pars = np.array([ecc,inc,ome,tp,a,argperi,
                           stellar_mass,params['rp'],params['d_obs']])
    orbit = Orbit(orbit_pars) 

    # Compute self-consistent instantaneous planet-star separation, 
    # Phi, and working angle for desired phase angle 
    # phase function, and orbit parameters 
    # add them to the parameter dictionary with the 
    # change_alpha() function
    alpha = 65.0*np.pi/180.0 # phase angle in radians
    pfunc = mu.avg_empirical # type of phase function to use
    change_alpha(orbit, params, alpha, pfunc) # updates params['sep'], params['phi'], params['wa']

    # Now need to add the desired coronagraph versions
    # this also adds needed astronomical/physical constants 
    # and some things related to assumed calibration method
    # the coronagraph/detector options are: 
    # org_hlc_pars,org_spc_pars  (from before WIETER)
    # cbe_hlc_pars, cbe_spc_pars ("current best estimates" as of August 2018)
    # req_spc_pars, req_hlc_pars (mission requirements as of August 2018)

    params.update(mu.cbe_spc_pars)

    # Now generate some mock IFS data 
    TEXP = 500.0        # hours
    bandcenter = 0.7    # microns
    width = 18.0        # the full width divided by bandcenter * 100.0
    wl, truth, random, noise = generate_noisey_spectra(TEXP, bandcenter, width, params)
    plt.plot(wl,truth,'ko')
    plt.plot(wl,random,'ro')
    plt.errorbar(wl,random,yerr=noise)

    plt.show()

    # to do imaging data, one should change the coronagraph to the HLC
    params.update(mu.cbe_hlc_pars)
    # NOTE*** For use with retrievals, the following imaging band centers and widths 
    # are currently hardcoded into juneper_model.py and coolTLUSTY_model.py 
    # centers = np.array([0.506,0.575,0.661,0.721,0.883,0.940])
    # widths = np.array([0.103,0.101,0.10,0.050,0.052,0.060])

    TEXP = 100.0
    bandcenter = 0.506
    width = 10.3
    wl, truth, random, noise = generate_noisey_imaging(TEXP, bandcenter, width, params)
    plt.errorbar(wl,truth,xerr=bandcenter*width/200.0,color='k')
    plt.plot(wl,random,'rs')
    plt.plot(wl,truth,'ks')
    plt.errorbar(wl,random,yerr=noise)

    plt.show()

