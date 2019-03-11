import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
import numpy as np
from scipy.interpolate import interp1d, CubicSpline
from orbit_class import Orbit

def calc_wa(sep, alpha, d_obs):
    # sep should be in au
    # alpha should be in degrees
    # d_obs should be in parsecs
    # returns wa in radians
    sep = sep*au_to_m
    d_obs = d_obs*pc_to_m
    alpha = alpha*np.pi/180.0
    # dplanet = np.sqrt(sep**2.0+d_obs**2.0-2.0*sep*d_obs*np.cos(wa))
    # wa = np.arcsin(np.sin(alpha)*sep/dplanet)  
    wa = np.arcsin(np.sin(alpha)*sep/d_obs)  
    return wa 

def lambert_phase(alpha, degrees=False):
    if degrees:
        alpha =  alpha*np.pi/180.0
    return (np.sin(alpha) + (np.pi - alpha)*np.cos(alpha))/np.pi

def quantum_eff(wl):
  if type(wl) == np.ndarray:
    wlcopy = [wl[k] for k in range(len(wl))]
    wlcopy = np.array([wlcopy]).reshape(len(wl))
    for k in np.where(wlcopy<0.7)[0]:
        wlcopy[k] = 0.7
  elif wl < 0.7:
      wlcopy = 0.7
  else:
      wlcopy = wl
  return 0.9*(1.0-(wlcopy-0.7)/0.3)*0.9*0.98*0.865
  
def old_quantum_efficiency(wavelength):
    #input wavelength unit: microns
    if type(wavelength) == np.ndarray:
        wlcopy = [wavelength[k] for k in range(len(wavelength))]
        wlcopy = np.array([wlcopy]).reshape(len(wavelength))
        for k in np.where(wlcopy<0.7)[0]:
            wlcopy[k] = 0.7
    elif wavelength < 0.7:
        wlcopy = 0.7
    else:
        wlcopy = wavelength
    return (0.87*(1.0-(wlcopy-0.7)/0.32))#*0.9*0.98*0.865 # factors

def quantum_efficiency(wavelength):
    #input wavelength unit: microns
    wls, eff = np.loadtxt('AuxiliaryData/ev2_qe_curve.dat',unpack=True)
    qefunc = interp1d(wls/1000.0, eff/100.0)
    return qefunc(wavelength)

def create_wl_range(start,end,R):
    wlrange = [start]
    wl = start
    while wl < end:
        wlnew = wl + wl/R
        if wlnew < end:
            wlrange.append(wlnew)
        wl = wlnew
    wlrange =  np.array(wlrange) 
    return wlrange

def fixed_albedo(wavelength):
    if type(wavelength) == np.ndarray:
        return 0.25 + np.zeros(len(wavelength))
    else:
        return 0.25
        
def VIO(x):
    return 1.0 + x*(-1.815) + x**2.0*(0.940) + x**3.0*(-2.399) + x**4.0*(4.990) + x**5.0*(-2.715)
def BL1(x):
    return 1.0 + x*(-1.311) + x**2.0*(-2.382) + x**3.0*(5.893) + x**4.0*(-4.046) + x**5.0*(0.846)
def GRN(x):
    return 1.0 + x*(-1.507) + x**2.0*(-0.363) + x**3.0*(-0.062) + x**4.0*(2.809) + x**5.0*(-1.876)
def RED(x):
    return 1.0 + x*(-0.882) + x**2.0*(-3.923) + x**3.0*(8.142) + x**4.0*(-5.776) + x**5.0*(1.439)
def CB2(x):
    return 1.0 + x*(-1.121) + x**2.0*(-1.720) + x**3.0*(1.776) + x**4.0*(1.757) + x**5.0*(-1.691)
def CB3(x):
    return 1.0 + x*(-0.413) + x**2.0*(-6.932) + x**3.0*(11.388) + x**4.0*(-3.261) + x**5.0*(-1.783)
def avgMayorga(x):
    return np.mean(np.array([VIO(x),BL1(x),GRN(x),RED(x),CB2(x),CB3(x)]))

def IFSavgMayorga(x):
    return np.mean(np.array([GRN(x),RED(x),CB2(x),CB3(x)]))


def piecewise_empirical(wavelength, alpha, degrees=True):

    # currently piece-wise, perhaps should make some gradual connections?

    def VIO(x):
        return 1.0 + x*(-1.815) + x**2.0*(0.940) + x**3.0*(-2.399) + x**4.0*(4.990) + x**5.0*(-2.715)
    def BL1(x):
        return 1.0 + x*(-1.311) + x**2.0*(-2.382) + x**3.0*(5.893) + x**4.0*(-4.046) + x**5.0*(0.846)
    def GRN(x):
        return 1.0 + x*(-1.507) + x**2.0*(-0.363) + x**3.0*(-0.062) + x**4.0*(2.809) + x**5.0*(-1.876)
    def RED(x):
        return 1.0 + x*(-0.882) + x**2.0*(-3.923) + x**3.0*(8.142) + x**4.0*(-5.776) + x**5.0*(1.439)
    def CB2(x):
        return 1.0 + x*(-1.121) + x**2.0*(-1.720) + x**3.0*(1.776) + x**4.0*(1.757) + x**5.0*(-1.691)
    def CB3(x):
        return 1.0 + x*(-0.413) + x**2.0*(-6.932) + x**3.0*(11.388) + x**4.0*(-3.261) + x**5.0*(-1.783)

    x = alpha/180.0
    oldwl = wavelength
    phases = []
    if type(oldwl) == float or type(oldwl)==np.float64:
        wavelength = np.array([wavelength])

    for k in range(len(wavelength)):

        if wavelength[k] < 0.440:
            phi = VIO(x) 
        if wavelength[k] > 0.440 and wavelength[k] < 0.500:
            phi = BL1(x)
        if wavelength[k] > 0.500 and wavelength[k] < 0.575:
            phi = GRN(x)  
        if wavelength[k] > 0.575 and wavelength[k] < 0.650:
            phi = 0.5*GRN(x) + 0.5*RED(x)  
        if wavelength[k] > 0.650 and wavelength[k] < 0.740:
            phi = RED(x)  
        if wavelength[k] > 0.740 and wavelength[k] < 0.9:
            phi = CB2(x)  
        if wavelength[k] > 0.9:
            phi = CB3(x)  

        phases.append(phi)

    return np.array(phases)

def avg_empirical(alpha, degrees=False):
    def VIO(x):
        return 1.0 + x*(-1.815) + x**2.0*(0.940) + x**3.0*(-2.399) + x**4.0*(4.990) + x**5.0*(-2.715)
    def BL1(x):
        return 1.0 + x*(-1.311) + x**2.0*(-2.382) + x**3.0*(5.893) + x**4.0*(-4.046) + x**5.0*(0.846)
    def GRN(x):
        return 1.0 + x*(-1.507) + x**2.0*(-0.363) + x**3.0*(-0.062) + x**4.0*(2.809) + x**5.0*(-1.876)
    def RED(x):
        return 1.0 + x*(-0.882) + x**2.0*(-3.923) + x**3.0*(8.142) + x**4.0*(-5.776) + x**5.0*(1.439)
    def CB2(x):
        return 1.0 + x*(-1.121) + x**2.0*(-1.720) + x**3.0*(1.776) + x**4.0*(1.757) + x**5.0*(-1.691)
    def CB3(x):
        return 1.0 + x*(-0.413) + x**2.0*(-6.932) + x**3.0*(11.388) + x**4.0*(-3.261) + x**5.0*(-1.783)
    if degrees:
        alpha = alpha*np.pi/180.0
    alpha = alpha*180.0/np.pi
    x = alpha/180.0
    phi = (VIO(x)+BL1(x)+GRN(x)+RED(x)+CB2(x)+CB3(x))/6.0
    return phi


def update_c_version(new_param_dict,old_param_dict,cperformance_table,tablesampling):
    old_param_dict.update(new_param_dict)
    rlamD,r_arcsec,Intensity,Contrast,coreThruput,PSFpeak,area_sqarcsec,occTrans = cperformance_table
    tau_occ = interp1d(rlamD,occTrans,fill_value=0.0,bounds_error=False) # working angle dependent occulter transmission
    tau_core = interp1d(rlamD,coreThruput,fill_value=0.0,bounds_error=False) # working angle dependendt 
    tau_PSF = interp1d(rlamD,coreThruput/occTrans,fill_value=0.0,bounds_error=False)
    contrast_func = interp1d(rlamD,Contrast,fill_value=0.0,bounds_error=False)
    Ipk_func = interp1d(rlamD,PSFpeak,fill_value=0.0,bounds_error=False)
    wl_d = r_arcsec[0]*D/(rlamD[0]*rad_to_arcsec) *10.0**6.0 # the wavelength in microns for which the contrast table file was computed
    mpixTable = np.mean(area_sqarcsec)*(np.pi/180.0/3600.0)**2.0/(tablesampling*wl_d*10.0**-6.0/D)**2 # a fixed value!!!
    coronagraph_dict = {'rlamD':rlamD,'r_arcsec':r_arcsec,'contrast_func':contrast_func, 'Ipk_func':Ipk_func,
          'area_sqarcsec':area_sqarcsec,'wl_d':wl_d,'mpixTable':mpixTable,'tau_occ':tau_occ,'tau_core':tau_core,}
    old_param_dict.update(coronagraph_dict)


COLORMAP = cm.nipy_spectral
CMAP2 = cm.nipy_spectral_r #cm.hot_r
plt.rcParams['font.family'] = 'serif'
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction']= 'in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['ytick.major.size'] = 9
plt.rcParams['xtick.major.size'] = 9
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.minor.size'] = 3.5
plt.rcParams['xtick.minor.size'] = 3.5
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 14

# astrophysical constants and standards
rad_to_arcsec = 206265.0  # convert radians to arcsec
pc_to_m = 30855152366503100.0  # convert parsecs to meters
au_to_m = 149590000000.0   # convert AU to meters 149590000000
rsun = 1392000000.0    # radius of the sun in meters
rjup = 69911000.0    # radius of jupiter in meters 69911000
h = 6.62607004 * 10.0**-34.0 # planck constant m^2 kg s^-1
kb = 1.38064852 *10.0**-23.0 # Boltzman constant m^2 kg s^-2 K^-1
sb = 5.67 *10.**-8.0 # steffan boltzman constant W m^-2 K^-4
c = 299792458.0              # speed of light m s^-1
F0V = 3.6*10.0**-8.0 # zero mag V band flux W/m^2/um

solar_spec = 'AuxiliaryData/kurucz_solar_spec_1AU.dat' # specific flux density W / m^2 / micron, at 1AU
ref_wl, ref_flambda = np.loadtxt(solar_spec,unpack=True,usecols=(0,1))
fsun_lambda = interp1d(ref_wl,ref_flambda) # units still specific flux density
MsunV = 4.83 # absolute V magnitude of the sun
MzV = 23.0 # surface brightness of zodiacal light at 1 AU
FsolV = 1.86*10.0**3.0 # solar V band flux W/m^2/um at 1AU
MezV = 22.0 # magnitudes arcsec^-2
Lsun = 3.828*10.0**26.0# W 
# things all versions of coronagraph have in common:
tau_fil = 0.9   # filter
tau_pol = 1.0   # polarizer
ast_constants_dict = { 'MsunV':MsunV,'MzV':MzV,'MezV':MezV, 'rjup':rjup,'rsun':rsun,
                       'fsun_lambda':fsun_lambda,'F0V':F0V,'FsolV':FsolV, 
                       'h':h,'c':c,'au_to_m':au_to_m,'pc_to_m':pc_to_m,
                       'rad_to_arcsec':rad_to_arcsec,'tau_fil':tau_fil,'tau_pol':tau_pol}

# relating to assumptions about calibration
# and image subtraction with other stars
delta_compstar = 3.0 # difference in apparent magnitudes between target and comparison star
fratio  = 10.0**(delta_compstar/2.5) # 
t_comp = 0.2 # 20% of the time on the target star
vratio = 1.0/(fratio*t_comp) # variance ratio target/comparison
zodi_multiplier = 1 + vratio 
detector_multiplier = 1 + t_comp # detector noise will increase with the additional time

ast_constants_dict.update({'zodi_multiplier':zodi_multiplier,'detector_multiplier':detector_multiplier})



## All the differenct versions of the coronagarph parameters... 
## org_hlc_pars,org_spc_pars,cbe_hlc_pars, cbe_spc_pars, req_spc_pars, req_hlc_pars
#############################Original SPC
# TELESCOPE
D = 2.38                        
tau_obs = 0.835                  
Apm = np.pi*(D/2.0)**2.0*tau_obs   
# DETECTOR 
eta = old_quantum_efficiency      
cte_cr_hot = 0.865*0.9*0.98  
phc = 0.90                     
tfr =  30.0                   
qcic = 0.01                  
idark =  2.1*10.0**-4.0            
sig_read = 0.00000001          
# CORONAGRAPH
cperformance_table = np.loadtxt("AuxiliaryData/SPC_20170714_660_1masStr_1masJit_CBEFIT.csv",unpack=True,delimiter=',')
tablesampling = 0.2 
tau_ref = 0.474
OWA = 2.7 
IWA = 9.0 
# IFS
wl_c = 0.6 
mpix = 20.0 
R = 50.0
# POST-PROCESSING
fpp = 1.0/10.0
# CREATE DICT
original_spc_dict = {'Apm':Apm,'eta':eta,'tfr':tfr,'qcic':qcic,'idark':idark,'phc':phc,'cte_cr_hot':cte_cr_hot,'R':R,
      'sig_read':sig_read,'wl_c':wl_c,'fpp':fpp,'OWA':OWA,'IWA':IWA,'D':D,'mpix':mpix,'tau_obs':tau_obs,'tau_ref':tau_ref}

org_spc_pars = ast_constants_dict.copy() 
update_c_version(original_spc_dict, org_spc_pars,cperformance_table,tablesampling)

########################### Original HLC
# TELESCOPE
D = 2.38                        
tau_obs = 0.835                  
Apm = np.pi*(D/2.0)**2.0*tau_obs   
# DETECTOR 
eta = old_quantum_efficiency      
cte_cr_hot = 0.865*0.9*0.98    
phc = 1.0                     
tfr =  30.0                   
qcic = 0.01                  
idark =  2.1*10.0**-4.0            
sig_read = 0.00000001  
ENF = np.sqrt(2.0)        
# CORONAGRAPH
cperformance_table = np.loadtxt("AuxiliaryData/hlc_20161228_9_polall_0.4mas_jitter_results.txt",unpack=True)
tablesampling = 0.2 
tau_ref = 0.474 
OWA = 2.7 
IWA = 9.0 
# focal plane
wl_c = 0.6 
mpix = 7.0 
R = 50000.0
# POST-PROCESSING
fpp = 1.0/10.0
# CREATE DICT
original_hlc_dict = {'Apm':Apm,'eta':eta,'tfr':tfr,'qcic':qcic,'idark':idark,'phc':phc,'cte_cr_hot':cte_cr_hot,'R':R,
      'sig_read':sig_read,'wl_c':wl_c,'fpp':fpp,'OWA':OWA,'IWA':IWA,'D':D,'mpix':mpix,'tau_obs':tau_obs,'tau_ref':tau_ref,'ENF':ENF}
org_hlc_pars = ast_constants_dict.copy() 
update_c_version(original_hlc_dict, org_hlc_pars,cperformance_table,tablesampling)

############################ SPC CBE
# TELESCOPE
D = 2.37                        
tau_obs = 0.835                  
Apm = np.pi*(D/2.0)**2.0*tau_obs   
# DETECTOR 
eta = quantum_efficiency      
cte_cr_hot = 0.934*0.874*0.983   # (CBE goes with 33% lifetime, 21 months at L2)
phc = 0.90                     
tfr =  80.0                   
qcic = 0.016                  
idark =  4.6*10.0**-4.0            
sig_read = 0.00000001          
# CORONAGRAPH
cperformance_table = np.loadtxt("AuxiliaryData/SPC_20170714_660_1masStr_1masJit_CBEFIT.csv",unpack=True,delimiter=',')
tablesampling = 0.2 
tau_ref = 0.383 
OWA = 2.7 
IWA = 9.0 
# IFS
wl_c = 0.66 
mpix = 45.0 #54.0 # 26.5
R = 50.0
# POST-PROCESSING
fpp = 1.0/12.0
# CREATE DICT
cbe_spc_dict = {'Apm':Apm,'eta':eta,'tfr':tfr,'qcic':qcic,'idark':idark,'phc':phc,'cte_cr_hot':cte_cr_hot,'R':R,
      'sig_read':sig_read,'wl_c':wl_c,'fpp':fpp,'OWA':OWA,'IWA':IWA,'D':D,'mpix':mpix,'tau_obs':tau_obs,'tau_ref':tau_ref}
cbe_spc_pars = ast_constants_dict.copy() 
update_c_version(cbe_spc_dict, cbe_spc_pars,cperformance_table,tablesampling)

########################## HLC CBE
# TELESCOPE
D = 2.37                        
tau_obs = 0.835                  
Apm = np.pi*(D/2.0)**2.0*tau_obs   
# DETECTOR 
eta = quantum_efficiency      # (CBE goes with 33% lifetime, 21 months at L2)
cte_cr_hot = 0.934*0.874*0.983   
phc = 0.90  
ENF=1.0                   
tfr =  10.0                   
qcic = 0.016                  
idark =  4.6*10.0**-4.0            
sig_read = 0.00000001          
# CORONAGRAPH
cperformance_table = np.loadtxt("AuxiliaryData/hlc_20161228_9_polall_0.4mas_jitter_results.txt",unpack=True)
tablesampling = 0.3 
tau_ref = 0.573 
OWA = 2.7 
IWA = 9.0 
# focal plane
wl_c = 0.508 
mpix = 9.0 #4.0 
R = 50000.0
# POST-PROCESSING
fpp = 1.0/12.0
# CREATE DICT
cbe_hlc_dict = {'Apm':Apm,'eta':eta,'tfr':tfr,'qcic':qcic,'idark':idark,'phc':phc,'cte_cr_hot':cte_cr_hot,'R':R,
      'sig_read':sig_read,'wl_c':wl_c,'fpp':fpp,'OWA':OWA,'IWA':IWA,'D':D,'mpix':mpix,'tau_obs':tau_obs,'tau_ref':tau_ref,'ENF':ENF}
cbe_hlc_pars = ast_constants_dict.copy() 
update_c_version(cbe_hlc_dict, cbe_hlc_pars,cperformance_table,tablesampling)

############################ SPC REQ
# TELESCOPE
D = 2.37                        
tau_obs = 0.835                  
Apm = np.pi*(D/2.0)**2.0*tau_obs   
# DETECTOR 
eta = quantum_efficiency      
cte_cr_hot = 0.596*0.838*0.95   # (REQ goes with 63 months at L2)
phc = 0.90                     
tfr =  80.0                   
qcic = 0.0232                  
idark =  5.56*10.0**-4.0            
sig_read = 0.00000001          
# CORONAGRAPH
cperformance_table = np.loadtxt("AuxiliaryData/newSPC.dat",unpack=True)
tablesampling = 0.2 
tau_ref = 0.337 
OWA = 2.8 
IWA = 8.6 
# IFS
wl_c = 0.66 # where the det is nyquist sampled
mpix = 53.0 
R = 50.0
# POST-PROCESSING
fpp = 1.0/12.0

req_spc_dict = {'Apm':Apm,'eta':eta,'tfr':tfr,'qcic':qcic,'idark':idark,'phc':phc,'cte_cr_hot':cte_cr_hot,'R':R,
      'sig_read':sig_read,'wl_c':wl_c,'fpp':fpp,'OWA':OWA,'IWA':IWA,'D':D,'mpix':mpix,'tau_obs':tau_obs,'tau_ref':tau_ref}

req_spc_pars = ast_constants_dict.copy() 
update_c_version(req_spc_dict, req_spc_pars,cperformance_table,tablesampling)

############################ HLC REQ
# TELESCOPE
D = 2.37                        
tau_obs = 0.835                  
Apm = np.pi*(D/2.0)**2.0*tau_obs   
# DETECTOR 
eta = quantum_efficiency      # (REQ goes with 63 months at L2)
cte_cr_hot = 0.798*0.95*0.988   
phc = 0.90  
ENF = 1.0                   
tfr =  6.0                   
qcic = 0.0232                  
idark =  5.56*10.0**-4.0             
sig_read = 0.00000001          
# CORONAGRAPH
cperformance_table = np.loadtxt("AuxiliaryData/newHLC.dat",unpack=True)
tablesampling = 0.3 
tau_ref = 0.573 
OWA = 2.8 
IWA = 8.6 
# focal plane
wl_c = 0.508  # where the det is nyquist sampled
mpix = 11.0 
R = 50000.0  
# POST-PROCESSING
fpp = 1.0/12.0
# CREATE DICT
req_hlc_dict = {'Apm':Apm,'eta':eta,'tfr':tfr,'qcic':qcic,'idark':idark,'phc':phc,'cte_cr_hot':cte_cr_hot, 'R':R,
      'sig_read':sig_read,'wl_c':wl_c,'fpp':fpp,'OWA':OWA,'IWA':IWA,'D':D,'mpix':mpix,'tau_obs':tau_obs,'tau_ref':tau_ref,'ENF':ENF}
req_hlc_pars = ast_constants_dict.copy() 
update_c_version(req_hlc_dict, req_hlc_pars,cperformance_table,tablesampling)

allversions = [org_hlc_pars,org_spc_pars,cbe_hlc_pars, cbe_spc_pars, req_spc_pars, req_hlc_pars]

