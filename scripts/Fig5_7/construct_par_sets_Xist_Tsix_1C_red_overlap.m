%construct par sets to simulate MA onset in cells with shorter overlap


%Read in 1C simulation of full overlap in model w/o X promoter repression,
%adapt par to simulation of 1C model with reduced overlap (k1 = k1/2)
par = dlmread('../../simulations/Fig5_7/sim_sys_bifurc_no_X_rep_170508.txt')
par = par(:,1:15);

par(:,1) = par(:,1)./2;

%Uncomment to construct parameter set
%dlmwrite('../../simulations/Fig5_7/par_1C_Xist_Tsix_8kb_Oct2018.txt', par);
