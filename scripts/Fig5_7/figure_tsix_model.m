clear all;
%% load simulation with 2C model

filename='../../simulations/Fig5_7/sim_2C_full_model_totXA2_bs_sets_rand_log_k12_13_14_15_complete_170829.txt';
sim{1}=dlmread(filename);

filename='../../simulations/Fig5_7/sim_2C_model_no_X_rep_totXA2_bs_sets_rand_log_k12_13_14_15_complete_170817.txt';
sim{2}=dlmread(filename);

filename='../../simulations/Fig5_7/sim_2C_model_wo_polcoll_totXA2_bs_rand_log_k12_13_14_15_complete_170829.txt';
sim{3}=dlmread(filename);
%% plotting parameters
graph_size=[2.5 2.5];
graph_size2=[2.5 1.5];
graph_size3=[4 2.5];

lw=1.5;
fs=8;
fst=10;
pos_x=[2:4.6:30];
pos_x2=[2 9.8 12.8 15.8];

pos_y=[4:3:20];
p1=[1800, 100,700,700];
figure(1)
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto','Units','pixels')

%% Bar graph: fraction XiXa parameter sets in reduced models -1C model
clf
filenames={'../../simulations/Fig5_7/sim_bifurc_sysparvar_all_rep_xixi_ro_170502.txt','../../simulations/Fig5_7/sim_sys_bifurc_part_T_sil_161118.txt',...
    '../../simulations/Fig5_7/sim_sys_bifurc_no_X_rep_170508.txt','../../simulations/Fig5_7/sim_sys_bifurc_no_pol_coll_161121.txt',...
    '../../simulations/Fig5_7/sim_sys_bifurc_only_sil_T_170424.txt','../../simulations/Fig5_7/sim_sys_bifurc_only_rep_X_by_T_170424.txt',...
    '../../simulations/Fig5_7/sim_sys_bifurc_only_pol_coll_170424.txt'};
ts=0.99;
for z=1:length(filenames)
    temp=dlmread(filenames{z});
    if z==2
        data{z}=temp(temp(:,11)==0,22);
    else
        data{z}=temp(:,22);
    end
    to_plot2(z)=100*sum(data{z}>ts)/length(data{z});
end

axes
bar(to_plot2,0.6,'k')
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(2) graph_size3],...
    'ylim',[0 100],'xlim',[0 8],'XTick',[1:7],'XTickLabel',[])
%ylabel('Parameter Sets [%]','Fontsize',fst)
%title('XaXi stable','Fontsize',fst)

%% Bar graph: fraction XiXa parameter sets in reduced models - 2C model

for a=1:3
    to_plot(a)=100*sum(mean(sim{a}(:,118:127),2)>ts)/size(sim{a},1);
end
to_plot_all=nan(1,7);
to_plot_all([1 3 4])=to_plot
axes
bar(to_plot_all,0.6,'FaceColor','k')
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(1) graph_size3],...
    'ylim',[0 4],'xlim',[0 8],'XTick',[1:7],'XTickLabel',[])
%ylabel('Parameter Sets [%]','Fontsize',fst)
%title('XaXa -> XaXi','Fontsize',fst)

print('../../plots/Fig5/Fig5D','-depsc','-loose')


%% example simulations with different extent of mono-allelic expression
% Compile model cpp function to mex file
% mex reaction_2C_wo_trans_XAtot2_k1XAtdep_161202.cpp

clf
pos_y3=[2 6.5];
pos_x3=[2:3:20]

sim_1_2=sim{1};
ma=mean(sim_1_2(:,118:127),2)>ts;
par_ma=find(ma);
max_bi=max(sim_1_2(ma,138:238),[],2);
time=1:100;
sw_time=sum(sim_1_2(ma,250:349).*repmat(time,sum(ma),1),2)+100*(1-sum(sim_1_2(ma,250:349),2));

tes=find(max_bi<0.3&max_bi>0.1);
q1=37;

axes
b=bar(0:100,100*[sim_1_2(par_ma(q1),138:238);sim_1_2(par_ma(q1),27:127)]','stacked');
b(1).FaceColor=[0.3 0.3 0.3];
b(1).EdgeColor=[0.3 0.3 0.3];
b(2).FaceColor=[0.7 0.7 0.7];
b(2).EdgeColor=[0.7 0.7 0.7];
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x3(4)+1 pos_y3(1) graph_size],...
    'xlim',[0 100],'XTick',[0:24:96],'XTickLabel',[],'YTick',[0 50 100]);
set(gca,'XTickLabel',[0:4]);

xlabel('Time [days]','Fontsize',fst);
ylabel('Cells [%]','Fontsize',fst);
title('100 cells','Fontsize',fst)

% plot single cells

v = 1/1440;
t_start=0;
t_diff=100;
t_before = 0;
output_time_step = 1;
sil_threshold = 10;
second_chr = 1;
k_adv_sil = 1;
seed=1;
const_par = [t_start, t_before, t_diff, sil_threshold, output_time_step, second_chr, k_adv_sil, seed];
nr_cells = 3;

p=sim_1_2(par_ma(q1),1:15);
ti={'Cell 1','Cell 2','Cell 3'};

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

[t,xp1,tp1,xr1,tr1,xp2,tp2,xr2,tr2,test, swo, sw_on_stable] = reaction_2C_wo_trans_XAtot2_k1XAtdep_161202(p, const_par, ...
    Xist_C1, Tsix_C1, Xist_RNA_C1, Tsix_RNA_C1, p_Xist_C1,...
    Xist_C2, Tsix_C2, Xist_RNA_C2, Tsix_RNA_C2, p_Xist_C2);


% plot simulation
col_x2=[177 206 85]/255;
col_x1=[35 155 56]/255;
axes
plot(t,xr1,'-','LineWidth',lw,'Color',col_x2);
hold on
plot(t,xr2,'-','LineWidth',lw,'Color',col_x1);

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x3(cells) 2 graph_size],...
    'xlim',[0 100],'XTick',[0:24:96],'XTickLabel',[0:4],'ylim',[-20 200],'YTick',[0:100:200],'YTickLabel',[]);

if cells==1
    ylabel('# Xist','Fontsize',fst)
    set(gca,'YTickLabel',[0:100:200]);
end
xlabel('Time [days]','Fontsize',fst)
title(ti{cells},'Fontsize',fst)


end

print('../../plots/Fig5/Fig5B_C','-depsc','-loose')

