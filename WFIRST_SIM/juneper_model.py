 import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
import misc_utils as mu

def jovian_Ag(wl):
    # input - wavelength in microns
    # output - measured Jupiter geometric albedo
    #          at that wavelength
    ref_wl, ref_ag = np.loadtxt('AuxiliaryData/jovian.albedo.dat',unpack=True, 
                                 usecols=(0,3),skiprows=11)
    ref_wl = ref_wl/1000.
    func = interp1d(ref_wl,ref_ag)
    return func(wl) 

def neptune_Ag(wl):
    # input - wavelength in microns
    # output - measured Neptune geometric albedo
    #          at that wavelength
    ref_wl, ref_ag = np.loadtxt('AuxiliaryData/jovian.albedo.dat',unpack=True, 
                                usecols=(0,6),skiprows=11)
    ref_wl = ref_wl/1000.
    func = interp1d(ref_wl,ref_ag,)
    return func(wl) 

def juneper_Ag(wavelengths,chromo=1.0, meth=1.0):
    # input:
    #   wavelengths - units of microns
    #   chromo - 1.0 means that short wavelength 
    #            behavior is entirely like Jupiter,
    #            0.0 means it is entirely like Neptune
    #   meth -   1.0 means that long wavelength 
    #            behavior is entirely like Jupiter,
    #            0.0 means it is entirely like Neptune  
    # ouptut:
    #   geometric albedo at the specified wavelengths
    jov = jovian_Ag(wavelengths) 
    nep = neptune_Ag(wavelengths)

    transition_lb = 0.55
    transition_ub = 0.65

    mask1 = np.where(wavelengths<transition_lb)[0]

    mask3 = np.where(wavelengths>transition_ub)[0]

    mask2 = [k for k in range(len(wavelengths)) if k not in mask1 and k not in mask3]
    mask2 = np.array(mask2)

    if len(mask2)>0:
        ptrans = np.linspace(chromo,meth,len(mask2))
        spec2 = ptrans*jov[mask2] + (1-ptrans)*nep[mask2]
    else:
        spec2=np.array([])

    spec1 = (chromo*jov + (1-chromo)*nep)[mask1] 
    spec3 = (meth*jov + (1-meth)*nep)[mask3] 
    
    Ag = np.concatenate([spec1,spec2,spec3])

    return Ag

def juneper_Agfunc(chromo, meth):
    # input:
    #   wavelengths - units of microns
    #   chromo - 1.0 means that short wavelength 
    #            behavior is entirely like Jupiter,
    #            0.0 means it is entirely like Neptune
    #   meth -   1.0 means that long wavelength 
    #            behavior is entirely like Jupiter,
    #            0.0 means it is entirely like Neptune  
    # ouptut:
    #   function for geometric albedo as a function of wavelength
    wavelengths = np.linspace(0.4,1.0,10000.0)
    jov = jovian_Ag(wavelengths) 
    nep = neptune_Ag(wavelengths)

    transition_lb = 0.55
    transition_ub = 0.65

    mask1 = np.where(wavelengths<transition_lb)[0]
    mask3 = np.where(wavelengths>transition_ub)[0]
    mask2 = [k for k in range(len(wavelengths)) if k not in mask1 and k not in mask3]
    mask2 = np.array(mask2)

    ptrans = np.linspace(chromo,meth,len(mask2))
    spec1 = (chromo*jov + (1-chromo)*nep)[mask1] 
    spec2 = ptrans*jov[mask2] + (1-ptrans)*nep[mask2]
    spec3 = (meth*jov + (1-meth)*nep)[mask3] 
    
    Ag = np.concatenate([spec1,spec2,spec3])

    return interp1d(wavelengths,Ag)

def compute_fluxratio(chromo, meth, rad, ast_dict, wavelengths):
    sep, phi = ast_dict['sep'], ast_dict['phi']
    planet_radius_meters = rad * mu.rjup
    star_planet_sep_meters = sep * mu.au_to_m
    Ag = juneper_Ag(wavelengths,chromo=chromo, meth=meth)
    # check for imaging band centers
    # and integrate appropriately
    centers = np.array([0.506,0.575,0.661,0.721,0.883,0.940])
    widths = np.array([0.103,0.101,0.10,0.050,0.052,0.060])
    imagemask,include = [],[]
    for k in range(len(centers)):
        ind = np.where(wavelengths==centers[k])
        if len(ind[0])>0:
            imagemask.append(ind[0][0])
            include.append(k)
    imagemask,include = np.array(imagemask),np.array(include)
    if len(imagemask)>0:
        for k in range(len(include)):
            ind = include[k]
            imageAgs = juneper_Ag(np.linspace(centers[ind]*(1-widths[ind]/2.0),centers[ind]*(1+widths[ind]/2.0),500),chromo=chromo, meth=meth)
            Ag[imagemask[k]] = np.mean(imageAgs)

    return Ag*phi*(planet_radius_meters/star_planet_sep_meters)**2.0 

def lnlike(p, wavelength, fluxratio, errors, ast_dict):
    chromo, meth, rad = p
    model = compute_fluxratio(chromo, meth, rad, ast_dict,wavelength)
    inv_sigma2 = 1.0/(errors**2) 
    return -0.5*(np.sum((fluxratio-model)**2*inv_sigma2 - np.log(inv_sigma2)))


def lnprior(p):
    chromo, meth, rad = p
    if not 0.0 <= chromo <= 1.0:
        return -np.inf
    if not 0.0 <= meth <= 1.0:
        return -np.inf
    if not 0.0 <= rad <= 10.0:
        return -np.inf
    return 0.0

def lnprob(p, wavelength, fluxratio, errors, ast_dict):
    lp = lnprior(p)
    if not np.isfinite(lp):
        return -np.inf
    lnprob = lp + lnlike(p, wavelength, fluxratio, errors, ast_dict)

    return lnprob