function [] = function_sim_2C_wo_trans_k1XAtdep_hetXist_161214(inputFilename,outputFilename)
fprintf('scan_par_queue %s -> %s\n',inputFilename,outputFilename)

%inputFilename = '160922_transfer_bistability_to_monoallelic_upreg/interesting_Sets_sw_off1_sw_off_2_161102.txt';
par=dlmread(inputFilename);

output = [];
v = 1/1440;
t_start=0;
t_diff=100;
t_before = 0;
output_time_step = 1;
sil_threshold = 10;
second_chr = 1;
k_adv_sil = 1;
%Tsix mutants or Xce alleles: 1= complete KO, 0= no reduction
k2_red_X1 = 0;
k2_red_X2 = 0; 
%Strength of Dox induction of Xist: k1 = k1 + k1*k1_ind_X1 => 1=2*fold increase
%or for Xist KO: -1
k1_ind_X1 = -1;
k1_ind_X2 = 0;
%Time point of Dox induction of Xist
tp_dox = t_before + 0;
const_par = [t_start, t_before, t_diff, sil_threshold, output_time_step, second_chr, k_adv_sil, k2_red_X1, k2_red_X2, k1_ind_X1, k1_ind_X2,tp_dox];
nr_cells = 100;

for loop=1:size(par,1)
%p3 and p10 are meaningless here because no basal Xist promoter transitions exist, reversal rate of Tsix induced repression is solely determined by k8
p =par(loop,:);

ma = zeros((t_diff+t_before)/output_time_step+1,nr_cells);
ba = zeros((t_diff+t_before)/output_time_step+1,nr_cells);
switch_on = zeros((t_diff+t_before)/output_time_step+1,nr_cells);
mean_xist = zeros(1,nr_cells);
switch_off = zeros(6,nr_cells);
% First time point with 1 stable Xa (Xist < thresh, XA, Tsix active) and with 1 stable Xi (Xist > thresh, XA, Tsix silenced)
switch_on_stable = zeros(1,nr_cells);


for cells = 1:nr_cells

% initial conditions
Xist_C1 = zeros(229,1);
Xist_C2 = zeros(229,1);
Tsix_C1 = zeros(357,1);
Tsix_C2 = zeros(357,1);
p_Xist_C1 = 1;
p_Xist_C2 = 1;

r = floor(length(Tsix_C1)*v*p(2));
q = randsample(length(Tsix_C1),r);
Tsix_C1(q) = 1;
r = floor(length(Tsix_C2)*v*p(2));
q = randsample(length(Tsix_C2),r);
Tsix_C2(q) = 1;

Xist_RNA_C1 = 0;
Tsix_RNA_C1 = floor(p(2)/p(5));
Xist_RNA_C2 = 0;
Tsix_RNA_C2 = floor(p(2)/p(5));


[t,xp1,tp1,xr1,tr1,xp2,tp2,xr2,tr2,test, sw_on_stable] = reaction_2C_wo_trans_k1XAtdep_MUT_161214_minimal(p, const_par, ...
    Xist_C1, Tsix_C1, Xist_RNA_C1, Tsix_RNA_C1, p_Xist_C1,...
    Xist_C2, Tsix_C2, Xist_RNA_C2, Tsix_RNA_C2, p_Xist_C2);

%cells in which no switch off happens: 


output = [output; [p, xr1', xr2', sw_on_stable]];


end

end

dlmwrite(outputFilename,output)
fprintf('scan_par_queue %s -> %s done\n',inputFilename,outputFilename)
end
