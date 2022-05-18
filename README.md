# CRISPEY_SIMULATION

This repository contains a Jupyter/iPython Notebook that can be used to simulate a simplified CRISPEY screen 
or in principle any pooled screen with varying fitness values and a sequencing based readout. It assumes multiple
Poisson processes for culture bottlenecks, sequencing, and DNA extraction,  a binomial distribution on the probability
of success for each PCR cycle duplicating a readout molecule from extracted DNA, and a two-sided exponential distribution
for the fitness values and initial counts of each strain/edit being tracked in the librry. For now, fitness is calculated from simulations
using simple linear regression of logfreq~generations for each simulated strain. Replicate competitions can be added as well.
This model is used to get a rough estimate of power and precision for pooled screens given different parameter values for 
seuqencing depth, the coverage of the library at the genome extraction step, the number of cells passaged at each culture dilution,
the number of generations, the number of PCR cycles, PCR efficiency, and the scale of fitness effects. Helper functions for plotting and 
performing various statistical tests are in there as well. 
