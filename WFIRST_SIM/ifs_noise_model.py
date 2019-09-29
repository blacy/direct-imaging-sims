import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
import misc_utils as mu
from matplotlib import colors, ticker, cm
from orbit_class import Orbit
# -----------------------------------------------------------------
# Note that you must have 5 files in correct paths/formats/units for this
# to work correctly
# 1) contrast table for the coronagraph design
# 2) solar spectrum
# 3) geometric albedo spectrum for modeled planet
# 4) stellar spectrum for host star
# 5) quantum efficiency for e2_v detector
# may want to add something with phase function as function of wavlength...
# -----------------------------------------------------------------
def planet_star_fluxratio(wavelength, params):
    # input wavelength unit: microns 
    return params['Ag'](wavelength)*params['phi']*(params['rp']*params['rjup'] \
            / (params['sep']*params['au_to_m']))**2.0

def stellar_photon_flux(wavelength, params):
    # input wavelength unit: microns
    # units of flux returned are photons/m^2/second in bins centered at input wavelength
    # with resolving power R specified in params dictionary
    mstar = params['Mstar'] +5.0*np.log10(params['d_obs']/10.0) # distance modulus to get apparent V band magnitude
    flux = (wavelength*10.0**-6.0)/(params['h']*params['c']) \
           * params['fstar_lambda'](wavelength) * wavelength/params['R']\
           * 10.0**(-mstar/2.5) # modify flux from apparent mag of 0 to target's apparent mag
    return flux

def zodi_photon_flux(wavelength, params):
    # input wavelength unit: microns
    # units of flux returned are photons/m^2/second in bins centered at input wavelength
    # with resolving power R specified in params dictionary
    ratio = params['fstar_lambda'](wavelength)*10.0**(-params['Mstar']/2.5) \
            *(10.0*params['pc_to_m']/au_to_m)**2.0 / params['FsolV']
    exozodi_photon_flux_density = params['Nez']*ratio*params['F0V']*wavelength/params['R']  \
                               *10.0**(-params['MezV']/2.5)*params['sep']**-2.0 \
                               * (wavelength*10.0**-6.0)/(params['h']*params['c']) #\
                               #* (10.0/params['d_obs'])**2.0  
    ratio = params['fsun_lambda'](wavelength) / params['FsolV']
    zodi_photon_flux_density =  ratio*params['F0V']* wavelength/params['R']  \
                               *10.0**(-params['MzV']/2.5) \
                               * (wavelength*10.0**-6.0)/(params['h']*params['c']) 
    return (exozodi_photon_flux_density + zodi_photon_flux_density)*calc_psf_area(wavelength,params)

def calc_mpix(wavelength, params):
    # input wavelength unit: microns
    mpix = params['mpix']*(wavelength/params['wl_c'])**2.0 
    return mpix

def calc_mpixTable(wavelength, params):
    mpix = params['mpixTable']#*(wavelength/params['wl_d'])**2.0  
    return mpix

def calc_psf_area(wavelength,params):
    # wavelength in microns
    # output area in square arcseconds
    return np.median(params['area_sqarcsec'])*(wavelength/params['wl_d'])**2.0

# -----------------------------------------------------------------
def calc_rplanet(wavelength, params):
    # input wavelength unit: microns
    # returns photons/second in bins centered at input wavelength
    # with resolving power R specified in params dictionary
    wa_lamD = params['wa']/(wavelength*10.0**-6.0/params['D'])
    tau_pla = params['tau_ref']*params['tau_fil']*params['tau_core'](wa_lamD)*params['tau_pol']
    rpl = stellar_photon_flux(wavelength, params)*planet_star_fluxratio(wavelength, params) \
         *params['Apm']*tau_pla*params['eta'](wavelength)*params['cte_cr_hot']*params['phc']
    return rpl 

def calc_rzodi(wavelength, params):
    # input wavelength unit: microns
    # returns photons/second in bins centered at input wavelength
    # with resolving power R specified in params dictionary 
    wa_lamD = params['wa']/(wavelength*10.0**-6.0/params['D'])
    tau_zod = params['tau_occ'](wa_lamD)*params['tau_ref']*params['tau_fil']*params['tau_pol']
    rzodi = zodi_photon_flux(wavelength, params)*params['Apm']*tau_zod* \
            params['eta'](wavelength)*params['cte_cr_hot']*params['phc']
    return rzodi

def calc_rspeckle(wavelength, params): 
    # input wavelength unit: microns
    # returns photons/second in bins centered at input wavelength
    # with resolving power R specified in params dictionary
    wa_lamD = params['wa']/(wavelength*10.0**-6.0/params['D'])
    tau_spe = params['tau_ref']*params['tau_fil']*params['tau_pol']
    rsp = stellar_photon_flux(wavelength, params) \
          *params['contrast_func'](wa_lamD)*params['Ipk_func'](wa_lamD)\
          *calc_mpixTable(wavelength,params)\
          *params['Apm']*tau_spe*params['eta'](wavelength)*params['cte_cr_hot']*params['phc'] 
    return rsp

def calc_rdark(wavelength, params):
    # input wavelength unit: microns
    # returns photons/second in bins centered at input wavelength
    # with resolving power R specified in params dictionary
    return params['idark']*calc_mpix(wavelength,params)*params['cte_cr_hot']*params['phc']

def calc_rcic(wavelength, params):
    # input wavelength unit: microns
    # returns photons/second in bins centered at input wavelength
    # with resolving power R specified in params dictionary
    return params['qcic']*calc_mpix(wavelength,params)/params['tfr']*params['cte_cr_hot']*params['phc']

def calc_rread(wavelength, params):
    # input wavelength unit: microns
    # returns photons/second in bins centered at input wavelength
    # with resolving power R specified in params dictionary
    return params['sig_read']*calc_mpix(wavelength,params)/params['tfr']*params['cte_cr_hot']

# -----------------------------------------------------------------
def calc_rates(wavelength,params):
    rp = calc_rplanet(wavelength, params)
    rz = calc_rzodi(wavelength, params)
    rs = calc_rspeckle(wavelength,params)
    rd = calc_rdark(wavelength, params)
    rcic = calc_rcic(wavelength, params)
    rr = calc_rread(wavelength, params)
    rn = rp + (rz + rs)*params['zodi_multiplier'] + (rd + rcic + rr)*params['detector_multiplier']
    return rp,rz,rs,rd,rcic,rr,rn
 
def calc_snr(texp,wavelength,params):
    # input wavelength unit: microns
    # input texp unit: seconds
    rp,rz,rs,rd,rcic,rr,rn = calc_rates(wavelength,params)
    signal = rp*texp
    variance = np.sqrt(rn*texp + (params['fpp']*rs*texp)**2.0)
    return signal/variance

def calc_texp(SNR,wavelength,params):
    # input wavelength unit: microns
    # time returned in seconds
    rp,rz,rs,rd,rcic,rr,rn = calc_rates(wavelength,params)
    return SNR**2.0*rn / (rp**2.0 - (SNR*params['fpp']*rs)**2.0)

def calc_tSNRmax(wavelength,params,percent=0.9):
    # returns time in seconds to reach percent*SNRmax (recall time to 100% SNRmax is infinite)
    rp,rz,rs,rd,rcic,rr,rn = calc_rates(wavelength,params)
    SNRmax = rp/(params['fpp']*rs)
    return rn*(percent*SNRmax)**2.0/(rp**2.0*(1.0-percent**2.0))

def generate_noisey_data(TEXP, wavelengths, params):
    SNR = calc_snr(TEXP*60.0*60.0,wavelengths,params)
    noiseless_signal = planet_star_fluxratio(wavelengths,params)
    noise = noiseless_signal/SNR 
    noisy_signal = [np.random.normal(noiseless_signal[k],noise[k]) for k in range(len(wavelengths))]
    noisy_signal = np.array(noisy_signal)
    return noiseless_signal, noisy_signal, noise
    
def planet_star_fluxratio_rayleigh(wavelength, params):   
    ratio = params['Ag'](wavelength)*params['phi_func'](wavelength)*(params['rp']\
            *params['rjup']/(params['sep']*params['au_to_m']))**2.0
    return ratio

def generate_noisey_data_rayleigh(TEXP, wavelengths, params):
    SNR = calc_snr(TEXP*60.0*60.0,wavelengths,params)
    noiseless_signal = planet_star_fluxratio_rayleigh(wavelengths,params)
    noise = noiseless_signal/SNR 
    noisy_signal = [np.random.normal(noiseless_signal[k],noise[k]) for k in range(len(wavelengths))]
    noisy_signal = np.array(noisy_signal)
    return noiseless_signal, noisy_signal, noise
# -----------------------------------------------------------------
# ASTROPHYSICAL CONSTANTS AND CONVERSIONS
rad_to_arcsec = mu.rad_to_arcsec 
pc_to_m = mu.pc_to_m
au_to_m = mu.au_to_m
rsun = mu.rsun          # radius of the sun in meters
rjup = mu.rjup          # radius of jupiter in meters
h = mu.h                # planck constant m^2 kg s^-1
kb = mu.kb              # Boltzman constant m^2 kg s^-2 K^-1
sb = mu.sb              # steffan boltzman constant W m^-2 K^-4
c = mu.c                # speed of light m s^-1
F0V = mu.F0V            # 3.6*10.0**-8.0 # zero mag V band flux W/m^2/um
solar_spec = mu.solar_spec # specific flux density W / m^2 / micron, at 1AU
ref_wl, ref_flambda = np.loadtxt(solar_spec,unpack=True,usecols=(0,1))
fsun_lambda = interp1d(ref_wl,ref_flambda) # units still specific flux density
MsunV = mu.MsunV        # 4.83
FsolV = mu.FsolV        # solar V band flux W/m^2/um at 1AU
MzV = mu.MzV            # 23.0
MezV = mu.MezV          #22.0

# -----------------------------------------------------------------
if __name__ == '__main__':

    # ------------------------ EXAMPLE  ---------------------------------
    # SET UP TARGET SYSTEM
    d_obs = 10.0
    mstar = 5.0  #APPARENT vband magnitude
    Mstar = mstar - 5.0*np.log10(d_obs/10.0) # distance modulus to get ABSOLUTE V band magnitude
    stellartype = 'g0v'
    stellar_spec = 'AuxiliaryData/'+stellartype+'.dat' # specific flux density W / m^2 / micron, for zero mag star
    ref_wl, ref_flambda = np.loadtxt(stellar_spec,unpack=True,usecols=(0,1))
    fstar_lambda = interp1d(ref_wl,ref_flambda) # units still specific flux density
    sep = 3.8
    rp = 1.0
    Nez = 1.0
    ## set phi*Ag = 0.25
    Ag = mu.fixed_albedo # this function just returns 0.25 for every wavelength
    alpha = 65.0  # will only be used for wa, not Phi(alpha)
    phi = 1.0 
    wa = mu.calc_wa(sep,alpha,d_obs) # working angle in radians, consistent with alpha, sep, and circular orbit
    # PUT IT ALL IN DICTIONARY          
    params = {'fstar_lambda':fstar_lambda,'d_obs':d_obs,'Nez':Nez,
              'rp':rp,'sep':sep,'Ag':Ag,'phi':phi,'wa':wa,'Mstar':Mstar}
    # ADD CORONAGRAPH, DETECTOR AND TELESCOPE TO DICTIONARY
    # available coronagraph versions are:
    # org_hlc_pars, org_spc_pars, cbe_hlc_pars, 
    # cbe_spc_pars, req_spc_pars, req_hlc_pars
    cversion_dict = mu.cbe_spc_pars.copy()
    params.update(cversion_dict)


    # COMPUTE + PRINT SOME RATES AND TIMES 
    wavelength = np.linspace(0.575,1.0,1000) 
    times = np.arange(10.0,420.0,20.0)
    rp, rz, rs, rd, rcic, rr, rn = calc_rates(wavelength, params)
    times = calc_texp(5.0,wavelength,params)
    # index 435 is 0.76 microns, lets look at the center of that band
    SNRmax = rp[435]/(params['fpp']*rs[435])
    con = planet_star_fluxratio(0.76,params) 
    print('Planet-Star flux ratio: ', con)
    print('SNR max: ', SNRmax)
    print('wa (arcsec): ', wa*rad_to_arcsec)
    print('exposure time to SNR 5.0: %f hours' % (times[435]/60.0/60.0))
    print('count rates rp: %f rz: %f rs: %f rd: %f rclk: %f rr: %f  rn: %f all e-/sec/signalregion' %\
        (rp[435], rz[435], rs[435], rd[435], rcic[435], rr[435],rn[435]))



