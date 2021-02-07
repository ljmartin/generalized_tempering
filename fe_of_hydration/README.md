# Chemical_potential of a united atom methane

Generalized tempering can also be used to temper the nonbonded interactions between a ligand/particle of interest and its surroundings. 
Instead of tempering a harmonic restraint that restricts distance (like in the na_cl_pmf example), this time we temper a 'lambda' parameter,
which ranges between 0 and 1 and controls the scaling of the nonbonded forces between a particle of interest and a waterbox. To keep things simple,
we make the particle a single united atom methane. 

Like the ion distance example, tempering the nonbonded interactions also alters the configurations available
to the particle. Particularly at values close to zero, the particle gains a lot of entropy and can slip through water molecules without resistance. 
For this reason, we choose lots of points near zero, and fewer points near 1 where things don't change so quickly. This also helps smooth out sampling
across lambda space. When the particle is in the non-interacting state, it can overlap coordinates with water atoms. Monte-carlo moves that bring lambda 
back into non-zero territory will often be rejected if they cause the methane to pop into existence at coordinates that overlap with water due to the 
huge energy associated with the repulsive part of the LJ potential. To reduce this, we keep the jump from 'non-interacting' to 'slightly-interacting' as
small as possible. 


As you can see below, again the generalized tempering approach has good agreement with the MBAR estimate. Both are close to the estimate in 
J. Chem. Phys. 122, 134508 (2005); https://doi.org/10.1063/1.1877132 of 2.24 kcal (-4.3*0.59 = 2.5 kcal/mol, where 0.59 is kT in kcal/mol), but
do not aim for correspondence with other simulation setups, only requiring that MBAR and generalized tempering are self-consistent within the context
of our particular simulation setup. 

![chemical_potential](chemical_potential.png "Chemical potential of united atom Methane")
