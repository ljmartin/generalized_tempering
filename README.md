# generalized_tempering
Generalized serial simulated tempering

Enhanced sampling for a flexible range of specific or non-specific collective variables / reaction coordinates / configurational parameters.

# na_cl_pmf
In this dir, find two calculations of the PMF vs. separation distance for a sodium and chloride ion pair in explicit solvent. 

- The PMF is first generated using the (metropolized) Bennet Acceptance Ratio with the PyMBAR library. This is a typical workflow for calculating free energy across a reaction coordinate, and serves as an anchor point against which we can measure the accuracy of the generalized tempering weights.
- A second PMF is estimated using generalized serial simulated tempering + the Wang-Landau algorithm. Observe that the weights learned closely track the free energy as estimated by MBAR:

![nacl](./na_cl_pmf/na_cl_pmf.png "NaClPMF")

_Note the zero point is arbitrary, since free energies are compared using relative differences, not absolute values. Here the zero point is set to the arithmetic mean._ 

# chemical_potential a.k.a. alchemical annihilation

In this dir, find a similar comparison as above. This time, instead of separation distance, we temper 'lambda', which in this context refers to a value that ranges [0,1] and scales the nonbonded interactions of some particle or molecule. See the readme in that directory for some more details. 

![chemical_potential](./chemical_potential/chemical_potential.png, "Methane chemical potential")


