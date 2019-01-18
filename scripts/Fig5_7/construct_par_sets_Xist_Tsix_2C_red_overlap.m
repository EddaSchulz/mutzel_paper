%construct par sets to simulate MA onset in cells with shorter overlap

%Read in 1C simulation of shorter overlap
sim = dlmread('../../simulations/Fig5_7/sim_1C_Xist_Tsix_8kb_Oct2018.txt');

rel = find(sim(:,22)>=0.99);

%Combine bistable parameter sets with random sil delays & react. rates 
addpath('../Fig1');

pars = sim(rel,1:15);
range_min = [log10(1) log10(1) log10(0.1) log10(0.1)];

range_max = [log10(49) log10(49) log10(100) log10(100)];

nr_par = 500;
p = zeros(nr_par,15);
for i = 1:size(pars,1)
	p12 = 10.^(lhsu(range_min(1),range_max(1),nr_par));
    p13 = 10.^(lhsu(range_min(2),range_max(2),nr_par));
    p14 = 10.^(lhsu(range_min(3),range_max(3),nr_par));
    p15 = 10.^(lhsu(range_min(4),range_max(4),nr_par));
	p = repmat(pars(i,:),size(p,1),1);
	p(:,12) = floor(p12);
    p(:,13) = floor(p13);
    p(:,14) = p14;
    p(:,15) = p15;
    %Uncomment to construct parameter set
    %dlmwrite('../../simulations/Fig5_7/par_2C_Xist_Tsix_8kb_Oct2018.txt', p, '-append');
end
