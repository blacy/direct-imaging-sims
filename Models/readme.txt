Reflected Light Curve GUI: 

This GUI allows the user to specify the orbital parameters, radius,
and single scattering albedo of an observed planet, then generates 
planet-star flux ratio light curves for different scattering models. 
This interface makes use of the wxpython package (included in tar ball), 
so the user must be certain this is installed prior to using the GUI. 
Once pre-requisite packages are in place, to use the GUI, simply type: 

$ pythonw mygui.py

or 

$ ./mygui.py 

The analytic scattering models implemented here are described in 
Madhusudhan & Burrows 2012: "Analytic Models for Albedos, 
Phase Curves, and Polarization of Reflected Light from Exoplanets". 

---------------------------------------------------------------------
CoolTLusy Geometric Albedo Models: 

This class of models uses the atmosphere and spectral code 
CoolTLusty, which solves self-consistent atmospheres under 
stellar irradiation, using a detailed suite of thermochemical 
and opacity data, augmented to incorporate the ExoMol 
methane opacities. Given the uncertainty in the nature of 
haze and cloud condensate species in giant exoplanet 
atmospheres, for this class of models, we merely keep the 
temperatures below ~200 K, assume a uniform scattering cloud 
and a uniform distribution of an absorber. The scattering cloud 
has a scattering opacity set above a wavelength 0.84 microns at 
a constant 0.002 cm^2 g^-1 and below a wavelength of 0.84 microns 
it assumes a 1/lambda^2.5 behavior. The uniform haze absorber is 
taken to be Titan tholins with an assumed atomic weight of 100, 
a model particle size of 0.05 microns, and a number fraction of 
3.3 x 10^-10. Even with such a low abundance, the tholin haze can 
markedly affect the albedo at short wavelengths and serves as our 
chromophore. These specific numbers and constituents were chosen 
to fit Jupiter's albedo spectrum for an atmosphere with a metallicity 
of 0.5 dex, ~3.16 x solar elemental abundances, insolated with a 
blackbody solar spectrum at 5777 K. With this background model, we 
then varied only the metallicity to include solar, 10 x solar, and 
30 x solar. Intermediate metallicities are interpolations between 
these four computed CoolTLusty models. In this way, we have generated 
a simple model suite that crudely captures the possible metallicity 
(read methane) dependence of such exoplanet albedos. 

---------------------------------------------------------------------
Jupiter-Neptune Hybrid Geometric Albedo Models:

In this class of model the geometric albedo spectrum is formed as a 
hybrid of the observed geometric albedo spectra of Neptune and Jupiter. 
With this approach, we can construct a geometric albedo spectrum for a 
planet whose atmospheric properties lie between these extremes. While 
both planets' reflective properties are dominated by the presence of 
ammonia clouds, the differences in the albedo spectra of Jupiter and 
Neptune are largely determined redward of 0.6 microns by the abundance 
of gaseous methane in the atmosphere and blueward of 0.6 microns by the 
presence or absence of the chromophore. Because there is no 
well-established correlation between methane and Jupiter's chromophore, 
interpolation between Jovian and Neptunian properties is here done 
independently for these two regions. Specifically, for a chromophore-region 
Jovian character Pc and methane-region Jovian character Pm (i.e. 
Pc or Pm = 1 is Jovian and Pc or Pm = 0 is Neptunian in the region 
c or m, respectively). We interpolate the geometric albedo 
spectra of Jupiter and Neptune as follows: 

Ag(lambda) =  Ag_jup x Pc + Ag_nep x (1-Pc), for lambda < 0.55 microns

              Ag_jup x Pm + Ag_nep x (1-Pm), for lambda > 0.65 microns

We arbitrarily select the region between 0.55 microns and 0.65 microns 
to transition from the blue chromophore-dominated region to the red 
methane-dominated region using a linearly scaling weighted average of 
the two formulae provided in the above equation. By separating the 
Jovian character of the red and blue regions of the spectrum, we allow 
Pc and Pm to function loosely as metrics of the chromophore and methane 
content of an atmosphere. All provided models make use of the geometric 
albedo spectra measured by Karkoschka 1994.

