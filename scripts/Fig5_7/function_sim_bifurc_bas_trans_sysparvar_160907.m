function [] = function_sim_bifurc_bas_trans_sysparvar_160907(inputFilename,outputFilename)
fprintf('scan_par_queue %s -> %s ',inputFilename,outputFilename)

%inputFilename = 'parameters_wo_Tsix_rep_160818.txt';
par=dlmread(inputFilename);
%par = [100,90,1,0.1733,1.3868,1,1,1,1,1,1,1,1,1,1];
output = [];
output2 = [];
v = 1/1440;
t_start=0;
t_max=500;
output_time_step = 1;
sil_threshold = 10;
second_chr = 0;
k_adv_sil = 1;
overlap_X_T = 229;
%Does Tsix gene overlap with Xist promoter? 1= yes, 0=no
T_over_Xp = 1;



nr_cells = 100;

% Loop over all random parametersets
for loop=1:size(par,1)
%for loop=1:1
p =par(loop,1:15);

if T_over_Xp==0
    p(7)=1; %Avoids rep of Xist promoter by Tsix
end

bs = zeros(t_max/output_time_step+1,nr_cells);
mean_xist_xon = zeros(1,nr_cells);
mean_xist_xoff = zeros(1,nr_cells);
mean_xist_xon_vl = zeros(1,nr_cells);
mean_xist_xoff_vl = zeros(1, nr_cells);
switch_on_Xa = zeros(1,nr_cells);
switch_off_Xi = zeros(1,nr_cells);

for cells = 1:nr_cells

% initial conditions: Tsix ON, Xist OFF
Xist_C1 = zeros(229,1);
Xist_C2 = zeros(229,1);
Tsix_C1 = zeros(357-(229-overlap_X_T),1);
Tsix_C2 = zeros(357-(229-overlap_X_T),1);

%With basal prom trans => Xist prom =0!
p_Xist_C1 = 0;
p_Xist_C2 = 0;

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
%Tsix_silenced=1 => Tsix is silenced in beginning (Xi), Tsix_silenced = 0=> Tsix not silenced in beginning
Tsix_silenced = 0;

const_par = [t_start, t_max, sil_threshold, output_time_step, second_chr, k_adv_sil, Tsix_silenced, overlap_X_T];


[t,xp1,tp1,xr1,tr1,xp2,tp2,xr2,tr2,test] = reaction_bifurc_bas_trans_160907(p, const_par, ...
    Xist_C1, Tsix_C1, Xist_RNA_C1, Tsix_RNA_C1, p_Xist_C1,...
    Xist_C2, Tsix_C2, Xist_RNA_C2, Tsix_RNA_C2, p_Xist_C2);
    
   xoff = xr1;
   toff = tr1;
    
   mean_xist_xoff(cells) = mean(xr1(end-50:end));
   mean_xist_xoff_vl(cells) = mean(xr1(end-100:end-50));
   temp = find(xoff>10);
   if (~isempty(temp))
		switch_on_Xa(cells) = temp(1);
   else
		switch_on_Xa(cells)= t_max;
   end

   % temp = find(xr1>10);
% temp2 = find(xr2>10);
% if (~isempty(temp))
%     if (~isempty(temp2))
%         switch_on(min(temp(1),temp2(1)),cells) = 1;
%     else
%         switch_on(temp(1),cells) = 1;
%     end
% elseif (~isempty(temp2))
%     switch_on(temp2(1),cells) = 1;
% end   
   
   
    
%Initial conditions XIst ON, Tsix OFF   
Xist_C1 = zeros(229,1);
Tsix_C1 = zeros(357-(229-overlap_X_T),1);
p_Xist_C1 = 1;

r = floor(length(Xist_C1)*v*p(1)*0.5);
q = randsample(length(Xist_C1),r);
Xist_C1(q) = 1;

Xist_RNA_C1 = 0.5*floor(p(1)/p(4));
Tsix_RNA_C1 = 0;
%Only silence Tsix on Xi if the repressive model includes silencing of Tsix
%by Xist RNA
if p(13)<1000
    Tsix_silenced = 1;
end
const_par = [t_start, t_max, sil_threshold, output_time_step, second_chr, k_adv_sil, Tsix_silenced, overlap_X_T];

%const_par = [t_start, t_max, sil_threshold, output_time_step, second_chr, k_adv_sil];
%p2 = p;
%p2(12) = t_max;
%PROBLEM!
%p2(13) = 0.0001;
%p2(15) = 0.0001;

[t,xp1,tp1,xr1,tr1,xp2,tp2,xr2,tr2,test] = reaction_bifurc_bas_trans_160907(p, const_par, ...
    Xist_C1, Tsix_C1, Xist_RNA_C1, Tsix_RNA_C1, p_Xist_C1,...
    Xist_C2, Tsix_C2, Xist_RNA_C2, Tsix_RNA_C2, p_Xist_C2);

%switch_off(cells) = switch_time;
xon = xr1;
ton = tr1;
mean_xist_xon(cells) = mean(xr1(end-50:end));
mean_xist_xon_vl(cells) = mean(xr1(end-100:end-50));   
temp = find(xon<10);
   if (~isempty(temp))
		switch_off_Xi(cells) = temp(1);
   else
		switch_off_Xi(cells)= t_max;
   end

bs(:,cells) = ((xon>10)&(xoff<10));



end %End of loop over cells
%output=[output ;[p, mean(mean_xist), sum(ma,2)'/nr_cells, sum(ba,2)'/nr_cells, ...
%    sum(switch_on,2)'/nr_cells]];
bs_sum = sum(bs,2)'/nr_cells;
frac_bs=mean(bs_sum(end-50:end));
output = [output; [p, mean(mean_xist_xoff_vl), mean(mean_xist_xon_vl), mean(mean_xist_xoff), std(mean_xist_xoff), mean(mean_xist_xon), std(mean_xist_xon), frac_bs, mean(switch_on_Xa), mean(switch_off_Xi)]];

end % end of loop over all monoallelic parameter sets

dlmwrite(outputFilename,output)
fprintf('scan_par_queue %s -> %s done\n',inputFilename,outputFilename)
end
