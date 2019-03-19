import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
import numpy as np
from scipy.interpolate import interp1d
from orbit_class import Orbit
import misc_utils as mu   
import ifs_noise_model as snr_ifs
import imaging_noise_model as snr_im
import juneper_model as hpm  
import coolTLUSTY_model as cpm
import simulate_data as simdata

def alpha_tsnr_plot(params,orbit,wavelengths,snr,cbar=True):
    sep,phi,wa=params['sep'], params['phi'], params['wa'] 
    texps = []
    alphas = []
    for alpha in np.linspace(0.0,180.0,1000):
        alphas.append(alpha)    
        alpha_r = alpha*np.pi/180.0 # phase angle in radians
        pfunc = mu.avg_empirical # type of phase function to use
        simdata.change_alpha(orbit, params, alpha_r, pfunc) # updates params['sep'], params['phi'], params['wa']    
        wl_texps = snr_ifs.calc_texp(snr,wavelengths,params)
        wl_texps = wl_texps/60.0/60.0
        mask = np.where(wl_texps < 0.0)[0]
        subvalue = 10000000.0
        for k in mask:
            wl_texps[k] = subvalue
        texps.append(wl_texps)
    #mask out the negative texps as super high
    texps = np.array(texps)
    zz = np.array(texps) 
    xx, yy = np.meshgrid(wavelengths, alphas)
    minorLocator = ticker.AutoMinorLocator()
    pcm = plt.pcolor(xx, yy, zz,
                       norm=colors.LogNorm(vmin=1.0,vmax=10000000.0),
                       cmap=mu.CMAP2, rasterized=True)
    plt.xlabel('Wavelength ($\mu$m)',fontsize=17)
    plt.ylabel('Phase Angle, $\\alpha$',fontsize=17)
    if cbar:
        cb = plt.colorbar(pcm, extend='max',orientation='horizontal')
        cb.set_label('Exposure Time to reach SNR = 5 (hours)',fontsize=17)

    params['sep'], params['phi'], params['wa'] = sep, phi, wa        
        
        
def rp_tsnr_plot(params,wavelengths,snr,cbar=True):
    rplanet = params['rp']
    texps = []
    radii = np.linspace(0.1,2.0,1000)
    for rp in radii:
        params.update({'rp':rp})
        wl_texps = snr_ifs.calc_texp(snr,wavelengths,params)
        wl_texps = wl_texps/60.0/60.0
        mask = np.where(wl_texps < 0.0)[0]
        subvalue = 10000000.0
        for k in mask:
            wl_texps[k] = subvalue
        texps.append(wl_texps)
    #mask out the negative texps as super high
    texps = np.array(texps)
    zz = np.array(texps) 
    xx, yy = np.meshgrid(wavelengths, radii)
    minorLocator = ticker.AutoMinorLocator()
    pcm = plt.pcolor(xx, yy, zz,
                       norm=colors.LogNorm(vmin=1.0,vmax=10000000.0),
                       cmap=mu.CMAP2, rasterized=True )
    if cbar:
        cb = plt.colorbar(pcm, extend='max',orientation='horizontal')
        cb.set_label('Exposure Time to reach SNR = 5 (hours)',fontsize=17)
    plt.xlabel('Wavelength ($\mu$m)',fontsize=17)
    plt.ylabel('Planet Radius (R$_{jup}$)',fontsize=17)

    params['rp'] = rplanet
    
    
def max_snr_plot(param_list):
    pass

def snr_vs_time_imaging(param_list):
    pass

def snr_vs_time_ifs(param_list):
    pass

