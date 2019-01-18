# Summary of scripts for in Fig1

1. *model.m* and *model\_f0.m*: 
	Simulation scripts of the ODE models (*model.m* is used for integration over time, *model\_f0.m* is used for the steady state solution)

2. *scan\_mods.m*:
Simulates the 1 and the 2 regulator ODE models with randomly sampled parameter sets. For each parameter set the following simulations are performed:
	- Simulation of female cells starting from monoallelic initial conditions
	- Simulation of male cells starting from an Xist low state
	- Simulation of female cells starting from biallelic initial conditions
	- Simulation of male cells starting from an Xist high state
	For each initial condition the ode23tb solver is used to integrate the system from 0h to 100h. The final state of this simulation is then used to 
	solve the system for the steady state using the function fsolve.
	The steady state solutions are written into the output file *model\_scan.txt*.

3. *lhsu.m*
	The function which performs Latin Hypercube sampling to randomly sample the parameter space, is used in *scan\_mods.m*.


