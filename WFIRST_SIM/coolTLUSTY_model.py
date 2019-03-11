import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d
import misc_utils as mu

def ctlustyAg(metallicity):
    # input - metallicity of atmosphere
    #         in multiples of solar metallicity
    #         can range from 1 to 30 times
    # output - returns a fuction for 
    #          the geometric albedo as a 
    #          function of wavelength in microns
    metallicities = np.array([1.0,3.0,10.0,30.0])
    filenames = ['AuxiliaryData/fort.19.jupiter.solar.tholin4d_7.0.05.power2.5.0.85mic',
                 'AuxiliaryData/fort.19.jupiter.3solar.tholin4d_7.0.05.power2.5.0.85mic',
                 'AuxiliaryData/fort.19.jupiter.10solar.tholin4d_7.0.05.power2.5.0.85mic',
                 'AuxiliaryData/fort.19.jupiter.30solar.tholin4d_7.0.05.power2.5.0.85mic']
    wls, ags = [],[]
    for filename in filenames:
        wl, ag = np.loadtxt(filename,unpack=True,usecols=(0,2),)
        wls.append(wl[:650])
        ags.append(ag[:650])
    wls, ags = np.array(wls), np.array(ags)
    # case that metallicity not in range:
    if metallicity < metallicities[0] or metallicity > metallicities[3]:
        print('metallicity outside interpolation range, make sure it is between 1.0 and 30.0')
        return None
    # case that metallicity exactly on the grid:
    elif metallicity == metallicities[0] or metallicity == metallicities[1] or metallicity == metallicities[2]  or metallicity == metallicities[3] :
        index = np.where(metallicity == metallicities)[0]
        agcombo = ags[index].reshape(len(wls[0]))
    else:
        lb = np.max(np.where(metallicity > metallicities))
        ub = np.min(np.where(metallicity < metallicities))
        p1 = 1.0 - (metallicity - metallicities[lb])/(metallicities[ub]-metallicities[lb])
        p2 = 1.0 - (metallicities[ub] - metallicity)/(metallicities[ub]-metallicities[lb])
        agcombo = p1*ags[lb] + p2*ags[ub]
    func = interp1d(wls[0],agcombo)
    return func

def compute_fluxratio(met, rad, ast_dict, wavelengths):
    sep, phi = ast_dict['sep'], ast_dict['phi']
    planet_radius_meters = rad * mu.rjup
    star_planet_sep_meters = sep * mu.au_to_m
    Agfunc = ctlustyAg(met)
    Ag = Agfunc(wavelengths)
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
            imageAgs = Agfunc(np.linspace(centers[ind]*(1-widths[ind]/2.0),centers[ind]*(1+widths[ind]/2.0),500))
            Ag[imagemask[k]] = np.mean(imageAgs)

    return Ag*phi*(planet_radius_meters/star_planet_sep_meters)**2.0


def lnlike(p, wavelength, fluxratio, errors, ast_dict):
    met, rad = p
    model = compute_fluxratio(met, rad, ast_dict,wavelength)
    inv_sigma2 = 1.0/(errors**2) 
    return -0.5*(np.sum((fluxratio-model)**2*inv_sigma2 - np.log(inv_sigma2)))


def lnprior(p):
    met, rad = p
    if not 1.0 <= met <= 30.0:
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