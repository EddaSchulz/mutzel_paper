addpath('../Fig1/');
%% Find the sets with stable XaXi, unstable XiXi, unstable Xi, stable Xa in the ODE model 

% read in simulation results for model comparison
sim = dlmread('../../simulations/Fig2_3_4/tXA_cXR_model_scan_for_stoch_sim.txt');
%% classify each parameter set as mono-allelic, not bi-allelic, 
% expression in males starting from and off state and expression in males starting from and Xist on state
mono=(sim(:,34)./sim(:,35))<0.1|(sim(:,34)./sim(:,35))>10;
no_bi=(sim(:,38)<0.1*max(sim(:,34),sim(:,35)))|(sim(:,39)<0.1*max(sim(:,34),sim(:,35)));
male_off = sim(:,36)<0.1*max(sim(:,34),sim(:,35));
male_on = sim(:,40)<0.1*max(sim(:,34),sim(:,35));
%% analyze maintenance of XaXi state in  tXA-cXR model
ind=zeros(1,8);
m = 2; n = 3;
ind([m n])=1;
sets=find(ismember(sim(:,25:32),ind,'rows') & mono==1 & no_bi==1 & male_off==1 & male_on==1);

par = sim(sets,1:33);
% Remove all unnecessary parameters that play no role in stochastic model
par(:,[1,2,7:10,15:33])=0;
%Write high Xist ss from ODE simulation into p2 for definition of MA expression in stochastic
%model (devide by 0.1733 to account for the longer halflife in the
%stochastic simulation)
par(:,2) = max(sim(sets,34),sim(sets,35))./0.1733;
%Transition between silencing states occurs with rate 1/h
par(:,18)=1;

% Assign to each ODE parameterset 10 randomly sampled combinations of
% p7,p8,p21,p22,p23
nr_par = 10;
p = zeros(10,33);
for i = 1:size(par,1)
par_n = floor(lhsu(1*ones(2,1),20*ones(2,1),nr_par));
par_k = floor(10.^lhsu(log10(50)*ones(3,1),log10(500)*ones(3,1),nr_par));

p = repmat(par(i,:),size(p,1),1);
p(:,7:8) = par_n;
p(:,21:23)=par_k;
%Uncomment to construct parameter set
%dlmwrite('../../simulations/Fig2_3_4/par_sets_stoch_sim_Oct2018.txt', p, '-append');
end


