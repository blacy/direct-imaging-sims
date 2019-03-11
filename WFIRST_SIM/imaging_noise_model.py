import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import misc_utils as mu
import ifs_noise_model as spec_noise_model 
from matplotlib import colors, ticker, cm
# -----------------------------------------------------------------
# Altered spectra noise model for imaging by changing R to something high (500000.0), 
# generating wavelength grid over the band pass with that high R, summing rates over band 
# also need to vary number of pixels in PSF region
# and to incorporate ENF since won't be in photon counting mode
# -----------------------------------------------------------------
def calc_rates(wavelength,params):
    rp = np.sum(spec_noise_model.calc_rplanet(wavelength, params))
    rz = np.sum(spec_noise_model.calc_rzodi(wavelength, params))
    rs = np.sum(spec_noise_model.calc_rspeckle(wavelength,params))
    rd = np.mean(spec_noise_model.calc_rdark(wavelength, params))
    rcic = np.mean(spec_noise_model.calc_rcic(wavelength, params))
    rr = np.mean(spec_noise_model.calc_rread(wavelength, params))
    rn = params['ENF']*(rp + (rz + rs)*params['zodi_multiplier'] + (rd + rcic + rr/params['ENF'])*params['detector_multiplier'])
    return rp,rz,rs,rd,rcic,rr,rn
 
def calc_snr(texp,wavelength,params):
    # input wavelength unit: microns
    # input texp unit: seconds
    rp, rz,rs,rd,rcic,rr,rn = calc_rates(wavelength,params)
    signal = rp*texp
    variance = np.sqrt(rn*texp + (params['fpp']*rs*texp)**2.0)
    return signal/variance

def calc_texp(SNR,wavelength,params):
    # input wavelength unit: microns, needs to be correct resolution list (the R of image params) 
    # time returned in seconds
    rp,rz,rs,rd,rcic,rr,rn = calc_rates(wavelength,params)
    return SNR**2.0*rn / (rp**2.0 - (SNR*params['fpp']*rs)**2.0)

def calc_tSNRmax(wavelength,params,percent=0.9):
    # input wavelength unit: microns, needs to be correct resolution list (the R of image params) 
    # returns time in seconds to reach percent*SNRmax (recall time to 100% SNRmax is infinite)
    rp,rz,rs,rd,rcic,rr,rn = calc_rates(wavelength,params)
    SNRmax = rp/(params['fpp']*rs)
    return rn*(percent*SNRmax)**2.0/(rp**2.0*(1.0-percent**2.0))

def generate_noisey_image_data(TEXP, wavelengths, params):
    # input TEXP in hours
    # input wavelength in microns
    SNR = calc_snr(TEXP*60.0*60.0,wavelengths,params)
    noiseless_signal = np.mean(spec_noise_model.planet_star_fluxratio(wavelengths,params))
    noise = noiseless_signal/SNR 
    noisy_signal = np.random.normal(noiseless_signal,noise)
    return noiseless_signal, noisy_signal, noise

def generate_noisey_image_data_rayleigh(TEXP, wavelengths, params):
    SNR = calc_snr(TEXP*60.0*60.0,wavelengths,params)
    noiseless_signal = np.mean(spec_noise_model.planet_star_fluxratio_rayleigh(wavelengths,params))
    noise = noiseless_signal/SNR 
    noisy_signal = np.random.normal(noiseless_signal,noise)
    return noiseless_signal, noisy_signal, noise

# -----------------------------------------------------------------
if __name__ == '__main__':

    # parameter dictionary
    params = spec_noise_model.params.copy() # take in same example target as ifs_noise_model.py
    # change coronagraph to imaging coronagraph
    params.update(mu.cbe_hlc_pars)
    # -----------------------------------------------------------------
    bandcenter, width = 0.506, 10.3
    wavelengths_in_band1 = mu.create_wl_range(bandcenter*(1.0 - width/200.0), bandcenter*(1.0 + width/200.0), params['R'])   
    bandcenter, width = 0.575, 10.1
    wavelengths_in_band2 = mu.create_wl_range(bandcenter*(1.0 - width/200.0), bandcenter*(1.0 + width/200.0), params['R'])  
    bandcenter, width = 0.661, 10.0
    wavelengths_in_band6 = mu.create_wl_range(bandcenter*(1.0 - width/200.0), bandcenter*(1.0 + width/200.0), params['R'])  
    bandcenter, width = 0.883, 5.2
    wavelengths_in_band8 = mu.create_wl_range(bandcenter*(1.0 - width/200.0), bandcenter*(1.0 + width/200.0), params['R'])  
    bandcenter, width = 0.721, 5.0
    wavelengths_in_band7 = mu.create_wl_range(bandcenter*(1.0 - width/200.0), bandcenter*(1.0 + width/200.0), params['R'])  
    bandcenter, width = 0.940, 6.4
    wavelengths_in_band9 = mu.create_wl_range(bandcenter*(1.0 - width/200.0), bandcenter*(1.0 + width/200.0), params['R'])  
    bandcenter, width = 0.760, 18.0
    wavelengths_in_ifsband = mu.create_wl_range(bandcenter*(1.0 - width/200.0), bandcenter*(1.0 + width/200.0), params['R'])
    rates = []
    for wavelengths in [wavelengths_in_band1, wavelengths_in_band2, wavelengths_in_band6,  wavelengths_in_band8,wavelengths_in_band7,wavelengths_in_band9 ]:
        rates.append(calc_rates(wavelengths, params))
    rp,rz,rs,rd,rcic,rr,rn = np.array(rates).T
    print('Count Rates (pla, zod+exz, spe, dark, cic, read): ', rp[1], rz[1], rs[1], rd[1], rcic[1], rr[1]) ## the 0.575 micron band
    max1 = rp[1]/(rs[1]*params['fpp'])
    time1 = calc_texp(5,wavelengths_in_band2,params)/60.0/60.0
    print('Max_SNR: %.2f'% max1, 'ExpTime: %.2f'% time1, 'Bijan Max_SNR: %.2f hrs'%100.0 , 'Bijan ExpTime: %.2f hrs'% 0.4)

