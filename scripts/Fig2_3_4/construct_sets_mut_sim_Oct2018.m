 %% Use the 100 best RSS fits of the model w/o PP kinetic to the mouse in vivo data to simulate the mutants
 clear all;
 addpath('../Fig1');
 %Read in File with best fits to mouse in vivo data for simulation of
 %human, aneu- & polyploid cells
 par1 = dlmread('../../simulations/Fig2_3_4/par_best_fits_lik_mouse_in_vivo_for_mutant_sim_Oct2018.txt');
  %Read in File with best fits to mESC differentiation data for simulation
  %of Dox induction and BA experiment
 par2 = dlmread('../../simulations/Fig2_3_4/par_best_fits_lik_mESC_for_mutant_sim_Oct2018.txt');
 
 %Uncomment to construct file with sets for mutant simulation
%   par1(:,26)=1;
%   dlmwrite('../../simulations/Fig2_3_4/sets_2n1X_Oct2018.txt', par1(:,1:33));
%   par1(:,26)=2;
%   dlmwrite('../../simulations/Fig2_3_4/sets_WT_mouse_in_vivo_fits.txt', par1(:,1:33));
%   par1(:,26)=3;
%   dlmwrite('../../simulations/Fig2_3_4/sets_2n3X_Oct2018.txt', par1(:,1:33));
%   par1(:,26)=4;
%   dlmwrite('../../simulations/Fig2_3_4/sets_2n4X_Oct2018.txt', par1(:,1:33));
%   par1(:,33)=4;
%   dlmwrite('../../simulations/Fig2_3_4/sets_4n4X_Oct2018.txt', par1(:,1:33));
%   par1(:,26)=3;
%   par1(:,33)=3;
%   dlmwrite('../../simulations/Fig2_3_4/sets_3n3X_Oct2018.txt', par1(:,1:33));
%   par1(:,26)=0;
%   par1(:,33)=0;
%   p1_range = [0.01 0.99];
%   range_min = [p1_range(1)];
%   range_max = [p1_range(2)];
%   
%   for a = 1:size(par1,1)
%       p1 = lhsu(range_min,range_max,1);
%       par1(a,1) = p1;
%   end
%   
%   dlmwrite('../../simulations/Fig2_3_4/sets_cXR_damp_rand_p1_Oct2018.txt', par1(:,1:33));
%   
%   dlmwrite('../../simulations/Fig2_3_4/sets_mESC_exp_Oct2018.txt', par2(:,1:33));

