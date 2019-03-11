import matplotlib.pyplot as plt
import numpy as np
import emcee 
import time
from scipy.interpolate import interp1d       
import coolTLUSTY_model as pm
import misc_utils as mu             

def run_mcmc(spectra_list):
    nthreads = 4
    nwalkers = 8                                    # number of walkers to use (3 times ndim generally ok)
    burnin_steps = 100                                # 10% of produciton recommended, or until spread out
    production_run_steps = 1000                       # enough to make nice contours...
    tag = spectra_list + '.'

    # for use with CoolTlusty model
    header = open(spectra_list,'r').readlines()[:5]
    planet_pars = {'sep': float(header[0].split(' = ')[1]),
                   'phi': float(header[1].split(' = ')[1]),
                   'rp':  float(header[2].split(' = ')[1]),
                   'met': float(header[3].split(' = ')[1])}
    ndim = 2                                                                       # dimensions of model pm
    labels = ['metallicity','r$_p$']                                                  # names of model parameters
    pos0 = np.array([planet_pars['met'], planet_pars['rp']]) # initial guesses for walkers
    
    wavelength, fluxratio, noisyfluxratio, errors = np.loadtxt(spectra_list, unpack=True)    
 
    mask = np.where(errors!=np.inf)
    if len(mask[0])==0:
        print('Infinite errors only in data from '+spectra_list+'\n')
       	np.savetxt(tag+'failed._chain.dat',np.zeros(100)+999.0)
        return 1
    wavelength,	noisyfluxratio,	errors = wavelength[mask[0]], noisyfluxratio[mask[0]], errors[mask[0]]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, pm.lnprob, args=(wavelength, noisyfluxratio, errors, planet_pars), threads=nthreads)   
    # starting position:
    pos = [pos0 + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

    # burnin
    if burnin_steps>0:
        time0 = time.time()
        # burnin phase
        pos, prob, state  = sampler.run_mcmc(pos, burnin_steps)
        sampler.reset()
        time1=time.time()
        print("burnin time: %f" %(time1-time0))
    elif burnin_steps==0:
        print("no burnin requested")
    else:
        print("Warning: incorrect input for burnin steps, setting burnin to 0")

    time0 = time.time()
    # perform MCMC
    pos, prob, state  = sampler.run_mcmc(pos, production_run_steps)
    time1=time.time()
    print("production run time: %f"%(time1-time0))

    # the results:
    samples = sampler.flatchain
    np.savetxt(tag+'_chain.dat',samples)  
    return samples

