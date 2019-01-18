# Summary of scripts for Fig5,7

1. *construct\_par\_sets\_sys\_par\_var\_bifurc.m*:
Constructs the parameter sets for simulation of the Full Xist/Tsix model and all reduced models.

2. *function\_sim\_bifurc\_bas\_trans\_sysparvar\_160907.m*, *reaction\_bifurc\_bas\_trans\_160907.cpp*:
Scripts for simulation of the 1C Xist/Tsix model with the parameter sets constructed in 1.
*function\_sim\_bifurc\_bas\_trans\_sysparvar\_160907.m* loops through all cells and calls *reaction\_bifurc\_bas\_trans\_160907.cpp*,
which performs the actual Gillespie simulation. To allow matlab to call the C++ function, a C++ MEX file needs to be constructed in matlab with:
```
mex filename.cpp
```
where filename is the name of the source cpp file. This is also the case for all other *.cpp* files.

3. *construct\_par\_sets\_2C\_sim.m*
Selects the sets with stable XaXi from the simulations in 2. and combines them with randomly sampled silencing delays and reactivation rates for
simulation of the 2C model.

4. *function\_sim\_2C\_wo\_trans_XAtot2\_k1XAtdep\_161202.m*, *reaction\_2C\_wo\_trans\_XAtot2\_k1XAtdep\_161202.cpp*:
Scripts for simulation of the 2C Xist/Tsix model with the parameter sets constructed in 3. 
*function\_sim\_2C\_wo\_trans\_XAtot2\_k1XAtdep\_161202.m* loops through all cells and calls *reaction\_2C\_wo\_trans\_XAtot2\_k1XAtdep\_161202.cpp*,
which performs the actual Gillespie simulation.

5. *function\_sim\_2C\_wo\_trans\_k1XAtdep\_hetTsix\_161214.m*, *function\_sim\_2C\_wo\_trans\_k1XAtdep\_hetXist\_161214.m*, *function\_sim\_2C\_wo\_trans\_k1XAtdep\_homTsix\_161214.m* 
and *reaction\_2C\_wo\_trans\_k1XAtdep\_MUT\_161214\_minimal.cpp*:
Scripts for simulation of the Xist and Tsix mutants. The matlab function loops through all cells and calls *reaction\_2C\_wo\_trans\_k1XAtdep\_MUT\_161214\_minimal.cpp*,
which performs the actual Gillespie simulation.

6. *figure\_tsix\_model.m*:
Plots simulation of monoallelic Xist upregulation with Xist/Tsix model (Figure 5B-D).

7. *figure\_xist\_tsix\_mutants.m*:
Plots simulation of WT, het Tsix, hom Tsix and het Xist mutants (Figure 7A-F).

8. *construct\_par\_sets\_Xist\_Tsix\_1C\_red\_overlap.m*:
Constructs parameter sets for simulation of 1C Model with reduced overlap

9. *function\_1C\_Xist\_Tsix\_8kb\_overlap\_Oct2018.m*, *reaction\_bifurc\_bas\_trans\_XAtot2\_161130.cpp*:
Scripts for simulation of 1C model with reduced overlap.
*function\_1C\_Xist\_Tsix\_8kb\_overlap\_Oct2018.m loops* through all cells and calls *reaction\_bifurc\_bas\_trans\_XAtot2\_161130.cpp*,
which performs the actual Gillespie simulation.

10. *construct\_par\_sets\_Xist\_Tsix\_2C\_red\_overlap.m*:
Selects the sets with stable XaXi from the 1C simulation with reduced overlap in 9. and combines them with randomly sampled silencing delays and 
reactivation rates for simulation of the 2C model.

11. *function\_sim\_2C\_trunc\_Tsix\_8kb\_171005\_Oct2018.m*, *reaction\_trunc\_Tsix\_161202\_171005.cpp*
Scripts for simulation of the 2C model with reduced overlap.
*function\_sim\_2C\_trunc\_Tsix\_8kb\_171005\_Oct2018.m* loops through all cells and calls *reaction\_trunc\_Tsix\_161202\_171005.cpp*,
which performs the actual Gillespie simulation.

12. *figure\_reduced\_overlap.m*:
Plots simulation of XaXi maintenance and monoallelic Xist upregulation with reduced overlap (Figure S6B,D,E,F).

13. *figure\_model\_wo\_TI.m*:
Calculates and plots the activation thresholds for the full model and the model w/o TI using the function *find\_thresholds.m* (Figure S4B,C).


