function [] = function_sim_2C_trunc_Tsix_8kb_171005_Oct2018(inputFilename,outputFilename)
fprintf('scan_par_queue %s -> %s\n',inputFilename,outputFilename)

%inputFilename = '160922_transfer_bistability_to_monoallelic_upreg/interesting_Sets_sw_off1_sw_off_2_161102.txt';
par=dlmread(inputFilename);

output = [];
v = 1/1440;
t_start=0;
t_diff=100;
t_before = 10;
output_time_step = 1;
sil_threshold = 10;
second_chr = 1;
k_adv_sil = 1;
overlap_X_T_1 = 80;
overlap_X_T_2 = 80;
%Does Tsix gene overlap with Xist promoter? 1= yes, 0=no
T_over_Xp = 0;
const_par = [t_start, t_before, t_diff, sil_threshold, output_time_step, second_chr, k_adv_sil, overlap_X_T_1, overlap_X_T_2];
nr_cells = 100;

for loop=1:size(par,1)
%p3 and p10 are meaningless here because no basal Xist promoter transitions exist, reversal rate of Tsix induced repression is solely determined by k8
p =par(loop,:);

if T_over_Xp==0
    p(7)=1; %Avoids rep of Xist promoter by Tsix
end

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
Tsix_C1 = zeros(357-(229-overlap_X_T_1),1);
Tsix_C2 = zeros(357-(229-overlap_X_T_2),1);
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


[t,xp1,tp1,xr1,tr1,xp2,tp2,xr2,tr2,test, swo, sw_on_stable] = reaction_trunc_Tsix_161202_171005(p, const_par, ...
    Xist_C1, Tsix_C1, Xist_RNA_C1, Tsix_RNA_C1, p_Xist_C1,...
    Xist_C2, Tsix_C2, Xist_RNA_C2, Tsix_RNA_C2, p_Xist_C2);

ma(:,cells) = ((xr1>10)&(xr2<10))|((xr2>10)&(xr1<10));
ba(:,cells) = (xr1>10)&(xr2>10);
for i=1:6
	switch_off(i,cells)=swo(i);
end
switch_on_stable(cells) = sw_on_stable;
temp = find(xr1>10);
temp2 = find(xr2>10);
if (~isempty(temp))
    if (~isempty(temp2))
        switch_on(min(temp(1),temp2(1)),cells) = 1;
    else
        switch_on(temp(1),cells) = 1;
    end
elseif (~isempty(temp2))
    switch_on(temp2(1),cells) = 1;
end    
mean_xist(cells) = mean(xr1(round(length(xr1)/2):end)+xr2(round(length(xr1)/2):end));

end
output=[output ;[p, mean(mean_xist), sum(ma,2)'/nr_cells, sum(ba,2)'/nr_cells, ...
    sum(switch_on,2)'/nr_cells, sum(switch_off(1,:),2)/size(find(switch_off(1,:)>0),2), ...
    size(find(switch_off(1,:)>0),2)/nr_cells, sum(switch_off(2,:),2)/size(find(switch_off(2,:)>0),2), ...
    size(find(switch_off(2,:)>0),2)/nr_cells,sum(switch_off(3,:),2)/size(find(switch_off(3,:)>0),2), ...
    size(find(switch_off(3,:)>0),2)/nr_cells, sum(switch_off(4,:),2)/size(find(switch_off(4,:)>0),2), ...
    size(find(switch_off(4,:)>0),2)/nr_cells, sum(switch_off(5,:),2)/size(find(switch_off(5,:)>0),2), ...
    size(find(switch_off(5,:)>0),2)/nr_cells, sum(switch_off(6,:),2)/size(find(switch_off(6,:)>0),2), ...
    size(find(switch_off(6,:)>0),2)/nr_cells, sum(switch_on_stable)/size(find(switch_on_stable>0),2), ...
    size(find(switch_on_stable>0),2)/nr_cells]];

end

dlmwrite(outputFilename,output)
fprintf('scan_par_queue %s -> %s done\n',inputFilename,outputFilename)
end
