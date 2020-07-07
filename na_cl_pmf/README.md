# PMFs of sodium and chloride ion separation

Generalized serial simulated tempering can be used to calculate PMFs. In the same way as umbrella sampling 
requires adjacent biases to have configurational overlap, here to get good Monte Carlo sampling between states we 
also require that states have overlap. Because the weights learned correspond to free energy, when they are converged
you can see that the weights closely track the MBAR free energy calculation:

![nacl](./na_cl_pmf.png "Title")
