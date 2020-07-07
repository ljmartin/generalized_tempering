# generalized_tempering
Generalized serial simulated tempering

Enhanced sampling for a flexible range of non-specific collective variables / reaction coordinates / configurational parameters.

# na_cl_pmf
In this dir, find two calculations of the PMF vs. separation distance for a sodium and chloride ion pair in explicit solvent. First the PMF is generated using the (metropolized) Bennet Acceptance Ratio with the PyMBAR library. This PMF is compared to one determined by generalized serial simulated tempering + the Wang-Landau algorithm. 

Observe that the weights learned closely track the free energy as estimated by MBAR:

![nacl](./na_cl_pmf/na_cl_pmf.png "NaClPMF")



