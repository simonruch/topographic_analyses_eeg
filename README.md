# Perform topographic analyses of EEG data in Matlab
## About
Topographic analyses of EEG data provide multivariate tests that add to mass univariate statistics as implemented in fieldtrip.

The functions provided here allow performing the following tests:
 * TANOVA: topographic analyses of variance; These analyses allow testing whether the electric field on the scalp (and hence the cortical source/oscillator) systematically different between conditions. Dependent and independent samples tests are implemented.
 *  TCT: topographic consistency tests; These tests allow assessing the consistency of the electric field on the scalp across subjects within a condition, or across trials within a subject & condition
 *  GFP contrasts:


Functions are intended to fit into the fieldtrip ecosystem (https://www.fieldtriptoolbox.org/),


## Install
 * Download repository
 * Make the folder available in the Matlab search path
 * Make sure to have the fieldtrip toolbox installed if fieldtrip should be used along with these functions

The functions were tested on Matlab 2019b.

## References
General description of TCT and TANOVA: 

*Koenig, T., & Melie-García, L. (2009). Statistical analysis of multichannel scalp field data. In C. M. Michel, D. Brandeis, J. Wackermann, L. R. R. Gianotti, & T. Koenig (Eds.), Electrical Neuroimaging (pp. 169–190). Cambridge University Press. https://doi.org/10.1017/CBO9780511596889.009*

Detailed description of randomisation (with Matlab code): 

*Koenig, T., & Melie-García, L. (2010). A method to determine the presence of averaged event-related fields using randomization tests. Brain Topography, 23(3), 233–242. https://doi.org/10.1007/s10548-010-0142-1*

Detailed description of TANOVA: 

*Koenig, T., Kottlow, M., Stein, M., & Melie-García, L. (2011). Ragu: A free tool for the analysis of EEG and MEG event-related scalp field data using global randomization statistics. Computational Intelligence and Neuroscience, 2011, 1–14. https://doi.org/10.1155/2011/938925*
