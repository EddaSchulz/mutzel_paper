# Summary of scripts for Fig2,3,4

1. *scan\_tXA\_cXR\_ODE\_mod\_181024.m*: 
Simulates more sets with the tXA-cXR ODE model to get enough sets which fulfill necessary XCI conditions (XaXi stable, XiXi unstable, Xi unstable, Xa stable..)

2. *construct\_par\_sets\_stoch\_mod\_181024.m*:
Selects the sets from the ODE simulation and adapts them for the stochastic simulation (remove unnecessary parameters, add transition rate, 
adapt steady state level to new halflife, add scaling factors and silencing delays). Each of the selected ODE parameter sets is combined with 
10 random combinations of tXA and cXR silencing delays and Xist, cXR and tXA scaling factors.

3. *sim\_functions\_gil\_Oct2018.jl*, *stoch\_mod\_Oct2018.jl*:
Simulation scripts for the stochastic model (*Julia*), were run with the parameter sets saved in 'par\_sets\_stoch\_sim\_Oct2018.txt'. The results of these simulations are saved in 
the folder *stoch\_sim\_files/* and are called *TEMP\_stoch\_sim\_Oct2018\_0.txt* - *TEMP\_stoch\_sim\_Oct2018\_114.txt*.
They have the following structure: 
   - C1-33: Parameters; 
   - C34-274: Xist on X1 0-240h 0-240h; 
   - C275-515: Xist on X2 0-240h 

4. *construct\_summary\_files\_stoch\_sim.m*:
Summarizes the result files from 3. to one line per simulated parameter set using the function *calc\_ma\_ba\_swon.m* and writes this summary into the file *SUMMARY\_stoch\_sim\_Oct2018.txt*. This summary result file 
has the following structure: 
   - C1-33: Parameters; C34-274: Fraction of monoallelic cells 0-240h; 
   - C275-515: Fraction of biallelic cells 0-240h;
   - C516-756: Fraction of triallelic cells 0-240h;
   - C757-997: Fraction of tetraallelic cells 0-240h;
   - C998: Mean switch on time


   Additionally, it saves the simulation with all monoallelic parameter sets into the file *MA\_sets\_stoch\_sim_Oct2018.txt*, which has the same structure as the result files in 3.

5. *plots\_species.m*:
Generates the plots of Figure 3D-E, Figure S2 and Figure S1A. Determines the 100 parameter sets that fit the mouse in vivo data and the mESC in vitro differentiation data
the best and saves them in the files *par_best\_fits\_mouse\_in\_vivo\_for\_mutant\_sim\_Oct2018.txt* and *par\_best\_fits\_mESC\_for\_mutant\_sim\_Oct2018.txt*. It saves the temporal offset 
of the best fit with each respective parameter set in C10 and the sum of squared residuals with this offset in C15. *par\_best\_fits\_mouse\_in\_vivo\_for\_mutant\_sim\_Oct2018.txt*
are the parameter sets that are used for the simulation of aneuploid, polyploid and human cells and *par\_best\_fits\_mESC\_for\_mutant\_sim\_Oct2018.txt* are the sets that are used for
simulation of experiments (see 7.).

6. *construct\_sets\_mut\_sim\_Oct2018.m*:
Constructs the parameter sets for simulation of the experiments and the simulation of aneu-, polyploid and human cells.

7. *sim\_functions\_gil\_Oct2018.jl*, *stoch\_mod\_aneu\_Oct2018.jl*, *stoch\_mod\_poly\_XA_dil\_Oct2018.jl*, *stoch\_mod\_poly\_XA\_rep\_Oct2018.jl*, *stoch\_mod\_damp\_XR\_damp\_Oct2018.jl*, 
*stoch\_mod\_damp\_XR\_upreg\_Oct2018.jl*, *stoch\_mod\_BA\_exp\_Oct2018.jl*, *stoch\_mod\_dox\_ind\_Oct2018.jl*:
These are the simulation scripts of the polyploid, aneuploid, human cells and the simulation of experiments (doxycycline induction prior to differentiation and artificial 
induction of biallelic expression) written in *Julia*.

8. *plot\_stoch\_traj\_aneuploidies.m*:
Plots Simulation of WT and Aneuploid cells with stochastic model: Figure 2C-D,H-I,J.

9. 

   - *plot\_exp.m*: Plots the experimental data of the doxycycline induction and the artificial biallelic induction experiments: Figure 4C,4F,4G,4I and Figure S3D. Also plots the RNA-FISH counts in the embryo (Fig3C).
   - *plot\_exp\_pred.m*:
Plots Simulation of above experiments: Figure4B,G; S3A-B. 

10. *figure\_block\_feedbacks.m*:
Plots steady states on cell and allele level for full model and model with blocked positive or negative feedback (Fig2E-G).
The scripts *model\_chr.m*, *model\_mod.m* are the simulation scripts for the models with blocked feedback.
