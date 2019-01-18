clear
% load data of systematic parameter variation of k1,k2,k8
path='../../simulations/Fig5_7/';
filename=[path,'sim_bifurc_sysparvar_all_rep_xixi_ro_170502.txt'];
sim=dlmread(filename);

%% load data of systematic parameter variation with partial Tsix silencing
path='../../simulations/Fig5_7/';
filename=[path,'sim_sys_bifurc_part_T_sil_161118.txt'];
sim3=dlmread(filename);

%% load data of systematic parameter variation without collisions
path='../../simulations/Fig5_7/';
filename=[path,'sim_sys_bifurc_no_pol_coll_161121.txt'];
sim4=dlmread(filename);

%% load data of systematic parameter variation only with collisions
path='../../simulations/Fig5_7/';
filename=[path,'sim_sys_bifurc_only_pol_coll_170424.txt'];
sim5=dlmread(filename);

%% load data of systematic parameter variation only with repression
path='../../simulations/Fig5_7/';
filename=[path,'sim_sys_bifurc_only_rep_X_by_T_170424.txt'];
sim6=dlmread(filename);

%% load data of systematic parameter variation only with silencing
path='../../simulations/Fig5_7/';
filename=[path,'sim_sys_bifurc_only_sil_T_170424.txt'];
sim7=dlmread(filename);
%% load data of systematic parameter variation only with silencing
filename=[path,'sim_sys_bifurc_no_X_rep_170508.txt'];
sim8=dlmread(filename);

%% activation thresholds
low_thresh=0.01;
high_thresh=0.99;
k{1} = unique(sim(:,1))';
k{2} = unique(sim(:,2))';
k{8} = unique(sim(:,8))';
k{11} = unique(sim3(:,11))';
ti{1}='Full Mod wo basal';
ti{2}='No Silencing';
ti{3}='No collisions';

%% analyze k1: generate matrices to read out thresholds
bif_p=1;
var_p=[2 8];
no_bif_p=var_p(~ismember(var_p,bif_p));
[b,c,d]=unique(sim(:,no_bif_p),'rows');

sim_index=zeros(length(c),20);
c2=unique(d);
for d2=1:length(c2)
    test=find(d==c2(d2));
    [g h]=sort(sim(test,bif_p));
    sim_index(d2,:)=test(h)';
end

bif_p_mat = zeros(size(sim_index));
bif_p_mat(:)=sim(sim_index(:),bif_p);
mono_mat = zeros(size(sim_index));
mono_mat(:)=sim(sim_index(:),22);
xa_mat = zeros(size(sim_index));
xa_mat(:)=sim(sim_index(:),18);
xi_mat = zeros(size(sim_index));
xi_mat(:)=sim(sim_index(:),20);

ana_thresh=find_thresholds(mono_mat, bif_p_mat, xa_mat, high_thresh, low_thresh);
sel_sets=find((ana_thresh(:,10)>0)&~isnan(ana_thresh(:,5))&(ana_thresh(:,8)./ana_thresh(:,5)<2));

%% analyze k1: generate matrices to read out thresholds for sim3
bif_p=1;
var_p=[1 2 8 11];
no_bif_p=var_p(~ismember(var_p,bif_p));
[b,c,d]=unique(sim3(:,no_bif_p),'rows');

sim_index3=zeros(length(c),20);
c2=unique(d);
for d2=1:length(c2)
    test=find(d==c2(d2));
    [g h]=sort(sim3(test,bif_p));
    sim_index3(d2,:)=test(h)';
end

bif_p_mat3 = zeros(size(sim_index3));
bif_p_mat3(:)=sim3(sim_index3(:),bif_p);
mono_mat3 = zeros(size(sim_index3));
mono_mat3(:)=sim3(sim_index3(:),22);
xa_mat3 = zeros(size(sim_index3));
xa_mat3(:)=sim3(sim_index3(:),18);
xi_mat3 = zeros(size(sim_index3));
xi_mat3(:)=sim3(sim_index3(:),20);

ana_thresh3=find_thresholds(mono_mat3, bif_p_mat3, xa_mat3, high_thresh, low_thresh);
sel_sets3=find((ana_thresh3(:,10)>0)&~isnan(ana_thresh3(:,5))&(ana_thresh3(:,8)./ana_thresh3(:,5)<2));

%% analyze k1: generate matrices to read out thresholds for sim4
bif_p=1;
var_p=[1 2 8];
no_bif_p=var_p(~ismember(var_p,bif_p));
[b,c,d]=unique(sim4(:,no_bif_p),'rows');

sim_index4=zeros(length(c),20);
c2=unique(d);
for d2=1:length(c2)
    test=find(d==c2(d2));
    [g h]=sort(sim4(test,bif_p));
    sim_index4(d2,:)=test(h)';
end

bif_p_mat4 = zeros(size(sim_index4));
bif_p_mat4(:)=sim4(sim_index4(:),bif_p);
mono_mat4 = zeros(size(sim_index4));
mono_mat4(:)=sim4(sim_index4(:),22);
xa_mat4 = zeros(size(sim_index4));
xa_mat4(:)=sim4(sim_index4(:),18);
xi_mat4 = zeros(size(sim_index4));
xi_mat4(:)=sim4(sim_index4(:),20);

ana_thresh4=find_thresholds(mono_mat4, bif_p_mat4, xa_mat4, high_thresh, low_thresh);
sel_sets4=find((ana_thresh4(:,10)>0)&~isnan(ana_thresh4(:,5))&(ana_thresh4(:,8)./ana_thresh4(:,5)<2));

%% plotting parameters
graph_size=[2.5 2.5];
lw=1.5;
fs=8;
fst=10;
pos_x=[2:5:12 15];
pos_y=[2 :4.6:50];
col_x2=[177 206 85]/255;
col_x1=[35 155 56]/255;

p1=[1900, 100,700,500];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')

%%
% analyze k1: plot example bifurcation diagram

sel_par= find(sim(sim_index(:,1),8)==k{8}(7));
right_thresh = ana_thresh(sel_par,5);
right_thresh2 = ana_thresh(sel_par,8);
match_k2=sim(sim_index(sel_par,1),2);
sel_th=right_thresh(match_k2==k{2}(11));
sel_th2=right_thresh2(match_k2==k{2}(11));

col_thresh=[153 25 21]/255;
axes
sel_par= find(sim(sim_index(:,1),8)==k{8}(7)&sim(sim_index(:,1),2)==k{2}(11));
cols=[0.3 0.3 0.3;0.7 0.7 0.7];
bif_p=1;

%line([sel_th2 sel_th2],[-100 1000],'color',col_thresh2,'Linewidth',lw,'LineStyle',':')
e1=errorbar((bif_p_mat(sel_par,:)),sim(sim_index(sel_par,:),20),sim(sim_index(sel_par,:),21),'.-','color',col_x1,'Linewidth',lw)
hold on
e2=errorbar((bif_p_mat(sel_par,:)),sim(sim_index(sel_par,:),18),sim(sim_index(sel_par,:),19),'.-','color',col_x2,'Linewidth',lw)
line([sel_th sel_th],[-100 1000],'color',col_thresh,'Linewidth',lw,'LineStyle',':')
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(2) pos_y(1) graph_size],...
    'xlim',[0 300],'XTick',[0:100:300],'XTickLabel',[0:50:150],'ylim',[-100 1000]);

xlabel('k_X [h^{-1}]','Fontsize',fst);
ylabel('# Xist','Fontsize',fst);
[le1 le2]=legend([e1 e2],{'Xi','Xa'},'location','Northwest','Box','off');
le1.Position([1:3])=[le1.Position(1)*0.96 le1.Position(2)*1.05 le1.Position(3)*0.4];
for q=1:2
    set(le2(q),'Fontsize',fs);
end


% plot k1/k2 for threshold
axes
bins=linspace(-5,2,21);
plot_par=1;
a = histogram(log2(0.5*sim(:,1)./sim(:,2)),bins);
a2=a.Values;
a2_width=a.BinWidth;
a = histogram(log2(0.5*ana_thresh(:,5)./sim(sim_index(:,1),2)),bins);

b2=a.Values;
plot(bins(2:end)-a2_width/2,100*a2/sum(a2),'-k','Linewidth',lw)
hold on
plot(bins(2:end)-a2_width/2,100*b2/sum(b2),'-','Linewidth',lw,'Color',col_thresh)

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(3) pos_y(1) graph_size],...
    'xlim',[bins(1) bins(end)],'ylim',[0 30],'XTick',[-5:1:2]);
xlabel('k_X/k_T [log2]','Fontsize',fst)
ylabel('Parameter Sets [%]','Fontsize',fst)
[le le2]=legend({'all sets',['activation' 10 'threshold']},'location','Northeast','Box','off');
le.Position=[le.Position(1)*1.27 le.Position(2)*0.9 le.Position(3)*0.3 le.Position(4)*1.5];
for q=1:2
    set(le2(q),'Fontsize',fs);
end

%%
%%%%%%% suppl figure 4

%% transitions regions Full model
figure(2)
clf

sel_par= find(sim(sim_index(:,1),8)==k{8}(7));
right_thresh = ana_thresh(sel_par,5);
right_thresh2 = ana_thresh(sel_par,8);
match_k2=sim(sim_index(sel_par,1),2);
sel_th=right_thresh(match_k2==k{2}(11));
sel_th2=right_thresh2(match_k2==k{2}(11));
col_thresh1=[153 25 21]/255;
col_thresh2=[0.5 0.5 0.5];
axes
sel_par= find(sim(sim_index(:,1),8)==k{8}(7)&sim(sim_index(:,1),2)==k{2}(11));
hold on
line([sel_th sel_th],[-100 1000],'color',col_thresh1,'Linewidth',1,'LineStyle','-')
line([sel_th2 sel_th2],[-100 1000],'color',col_thresh2,'Linewidth',1,'LineStyle','-')

e1=errorbar((bif_p_mat(sel_par,:)),sim(sim_index(sel_par,:),20),sim(sim_index(sel_par,:),21),'.-','color',col_x1,'Linewidth',lw)
e2=errorbar((bif_p_mat(sel_par,:)),sim(sim_index(sel_par,:),18),sim(sim_index(sel_par,:),19),'.-','color',col_x2,'Linewidth',lw)
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(2) graph_size],...
    'xlim',[0 300],'ylim',[-100 1000],'XTick',[0:100:300],'XTickLabel',[0:50:150]);

xlabel('k_X [h^{-1}]','Fontsize',fst);
ylabel('# Xist','Fontsize',fst);
[le1 le2]=legend([e1 e2],{'Xi','Xa'},'location','Northwest','Box','off');
le1.Position([1:3])=[le1.Position(1)*0.88 le1.Position(2)*1.05 le1.Position(3)*0.4];
for q=1:2
    set(le2(q),'Fontsize',fs);
end
title('Full Model','Fontsize',fst);


%% transitions regions wo TI

sel_par= find(sim4(sim_index4(:,1),8)==k{8}(7));
right_thresh_4 = ana_thresh4(sel_par,5);
right_thresh2_4 = ana_thresh4(sel_par,8);
match_k2_4=sim4(sim_index4(sel_par,1),2);
sel_th=right_thresh_4(match_k2_4==k{2}(11));
sel_th2=right_thresh2_4(match_k2_4==k{2}(11));
axes
sel_par= find(sim4(sim_index4(:,1),8)==k{8}(7)&sim4(sim_index4(:,1),2)==k{2}(11));
hold on
line([sel_th sel_th],[-100 1000],'color',col_thresh1,'Linewidth',1,'LineStyle','-')
line([sel_th2 sel_th2],[-100 1000],'color',col_thresh2,'Linewidth',1,'LineStyle','-')

e1=errorbar((bif_p_mat4(sel_par,:)),sim4(sim_index4(sel_par,:),20),sim4(sim_index4(sel_par,:),21),'.-','color',col_x1,'Linewidth',lw)
e2=errorbar((bif_p_mat4(sel_par,:)),sim4(sim_index4(sel_par,:),18),sim4(sim_index4(sel_par,:),19),'.-','color',col_x2,'Linewidth',lw)
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(2) pos_y(2) graph_size],...
    'xlim',[0 300],'ylim',[-100 1000],'XTick',[0:100:300],'XTickLabel',[0:50:150]);

xlabel('k_X [h^{-1}]','Fontsize',fst);
ylabel('# Xist','Fontsize',fst);
[le1 le2]=legend([e1 e2],{'Xi','Xa'},'location','Northwest','Box','off');
le1.Position([1:3])=[le1.Position(1)*1.03 le1.Position(2)*1.05 le1.Position(3)*0.4];
for q=1:2
    set(le2(q),'Fontsize',fs);
end
title('No Collisions','Fontsize',fst);





%%

axes
bins_a=linspace(1,6,5);
bins_a = [0 2 4 6];
temp=histogram(right_thresh2_4./right_thresh_4,bins_a);
a=temp.Values;
a_width=temp.BinWidth;

%bins_b=linspace(1,6,5);
bins_b = bins_a;
temp=histogram(right_thresh2./right_thresh,bins_b);
b=temp.Values;
b_width=temp.BinWidth;

%p=bar(bins(2:end)-a_width/2,[a./sum(a)*100;b./sum(b)*100]','grouped')
p=bar([bins_a(2:end)-a_width/2;bins_b(2:end)-b_width/2]',[a./sum(a)*100;b./sum(b)*100]','grouped')
p(1).FaceColor='k';
p(2).FaceColor=[0.6 0.6 0.6];
p(1).EdgeColor='k';
p(2).EdgeColor=[0.6 0.6 0.6];
p(1).BarWidth=1;

set(gca,'Fontsize',8,'Units','Centimeters','Position',[pos_x(1) pos_y(1) graph_size],...
    'xlim',[0 6],'ylim',[0 100],'XTick',[0:2:6],'TickLength',[0.02 0],...
    'TickDir','out','Linewidth',1,'YTick',[0:20:100])
xlabel('k_{X_{high}} / k_{X_{low}}','Fontsize',fst)
ylabel('Parameter Sets [%]','Fontsize',fst)
[le le2]=legend({'W/o TI','Full'},'location','Northeast','Box','off');
%le.Position=[le.Position(1)*2.1 le.Position(2) le.Position(3)*0.3 le.Position(4)*1.5];
for q=1:2
    set(le2(q),'Fontsize',fs);
end

print('../../plots/FigS4/FigS4B_C','-depsc')