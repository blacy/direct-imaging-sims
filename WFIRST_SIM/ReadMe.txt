Contents:
coolTLUSTY_model.py
coolT_mcmc.py
juneper_model.py
hybrid_mcmc.py
ifs_noise_model.py
imaging_noise_model.py
simulate_data.py
orbit_class.py
misc_utils.py
AuxiliaryData (a directory)
	versions of coronagraph
    stellar spectra
    geometric albdeo for different metallicities - made with CoolTLUSTY
    geometric albedo for gaseous solar system bodies - observed 
    select list of known RV planets

"Pipeline" Usage (demonstrated in example_sim_and_fit.ipynb):

1.) initiate a dictionary with the desired planet parameters and star parameters. 
    The Orbit object defined in orbit_class.py could be useful for this, along with
    possible phase functions and geometric albedo options in misc_utils.py and 
    coolTLUSTY_model.py or juneper_model.py
2.) detector + coronagraph parameters matching with different iterations of
    mission design are listed in miscellaneous_utils.py, select one of these and 
    add it to your dictionary
3.) select the imaging bands and/or IFS bands you want to include in your mock
    observations, along with an exposure time (in hours) and input these into
    generate_noisey_spectra() and/or generate_noisey_imaging() functions in simulat_data.py
    to obtain the true underlying signal, the expected noise level, and a randomly generated
    noisy iteration of your mock observation
3.) optionally, save your mock observation to a file in the correct format and
    perform a MCMC retrieval on your data using the functions in hybrid_mcmc.py
    or coolT_mcmc.py 

Other useful functionalities (demonstrated in example_SNR_calcs.ipynb):

- compute exposure times necessary to reach a specified SNR for a given system
- compute maximum possible SNR for a given system 
- compute SNR as a function of exposure time for a given system
- compare results across different coronagraph/detector specs 




