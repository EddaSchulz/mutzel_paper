%% Generate plot of Embryo data
clear
graph_size=[3.5  2.5];
lw=2;
fs=8;
fst=10;

p1=[1900, 100,700,500];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')

clf
data2=xlsread('../../data/Fig2_3_4/Data of E5.0.xlsx');
data2_norm=data2(:,2:4)./repmat(sum(data2(:,2:4),2),1,3);

cm=[140/255 198/255 62/255; 0.1529 0.667 0.882;0.5 0.5 0.5];
cm=[0.3 0.3 0.3;0.7 0.7 0.7;1 1 1];
b=bar(data2_norm*100,'stacked')
b(1).Parent.Parent.Colormap = cm;
b(1).BarWidth=0.5;
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[2 2 graph_size],...
    'xlim',[0 16],'ylim',[0 140],'xtick',[]);
xlabel('Embryos','Fontsize',fst)
ylabel('Cells [%]','Fontsize',fst)

plot_nr_em=[72, 76,32,52,58,48,36,52,58,76,50,76,53,68,42];
 for t=1:15
     tx=text(t,105,['n=',num2str(plot_nr_em(t))],'Fontsize',6)
     set(tx,'Rotation',90)
 end

print('../../plots/Fig3/Fig3C','-depsc')

%% Generate plots of experimental data in Fig. 4
clear
% plot settings
clear
p1=[1900, 100,600,700];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')
clf
graph_size=[2.5 2.5];
pos_x=[2:2.8:20];
pos_y=[2:4.5:20];
lw=2;
fs=8;
fst=10;
%% load data - bi-allelic induction
fish{1}=dlmread('../../data/Fig2_3_4/biallelic_FISH_exp1.txt');
fish{2}=dlmread('../../data/Fig2_3_4/biallelic_FISH_exp2.txt');
fish{3}=dlmread('../../data/Fig2_3_4/biallelic_FISH_exp9.txt');

bi_data=[fish{1}(1,:)./sum(fish{1});fish{2}(1,:)./sum(fish{2});fish{3}(1,:)./sum(fish{3})];
bi_mean=100*mean(bi_data);
bi_std=100*std(bi_data);

mono_data=[sum(fish{1}([2 4],:))./sum(fish{1});sum(fish{2}([2 4],:))./sum(fish{2});sum(fish{3}([2 4],:))./sum(fish{3})];
mono_mean=100*mean(mono_data);
mono_std=100*std(mono_data);
%% Plot 4F - bar graph of FISH patterns
clf
colormap([0 0 0; 0.8 0.8 0.8]);

axes
b=bar([mono_mean(3:5); mono_mean(6:8)]');
hold on;
xpos=[b(1).XData+b(1).XOffset;b(2).XData+b(2).XOffset];

plot(xpos(1,:),100*mono_data(:,3:5)','.','color',[0.5 0.5 0.5],'Markersize',8)
hold on
plot(xpos(2,:),100*mono_data(:,6:8)','.','color',[0.5 0.5 0.5],'Markersize',8)
errorbar(xpos',[mono_mean(3:5); mono_mean(6:8)]',[mono_std(3:5); mono_std(6:8)]','k.','Markersize',0.1,'Linewidth',0.7,'CapSize',2.5)

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(2) graph_size],...
    'XTickLabel',[6 24 48],'ylim',[0 100],'xlim',[0.5 3.5],'XTick',[1:3],'YTick',[0:25:100])
ylabel('Cells [%]','Fontsize',fst);
%xlabel('Time after Dox [h]','Fontsize',fst)
[tt_mono pval]=ttest2(mono_data(:,3:5),mono_data(:,6:8))

axes
b=bar([bi_mean(3:5); bi_mean(6:8)]');
hold on;
xpos=[b(1).XData+b(1).XOffset;b(2).XData+b(2).XOffset];

plot(xpos(1,:),100*bi_data(:,3:5)','.','color',[0.5 0.5 0.5],'Markersize',8)
hold on
plot(xpos(2,:),100*bi_data(:,6:8)','.','color',[0.5 0.5 0.5],'Markersize',8)
errorbar(xpos',[bi_mean(3:5); bi_mean(6:8)]',[bi_std(3:5); bi_std(6:8)]','k.','Markersize',0.1,'Linewidth',0.7,'CapSize',2.5)

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(2) pos_y(2) graph_size],...
    'XTickLabel',[6 24 48],'ylim',[0 100],'xlim',[0.5 3.5],'XTick',[1:3],'YTick',[])
[tt_bi pval]=ttest2(bi_data(:,3:5),bi_data(:,6:8))

[le1 le2]=legend({'Control','Dox'},'location','Northwest','Box','off');
le1.Position([1:3])=[le1.Position(1)*1 le1.Position(2)*1 le1.Position(3)*0.4];
for q=1:2
    set(le2(q),'Fontsize',fs);
end
 print('../../plots/Fig4/Fig4F','-depsc','-loose')
 
%% Fig. 4G - plot TRex data
clf
data1=readtable('../../data/Fig2_3_4/xist_fold.txt');
data1_b6=table2array(data1(1:2:18,4));
data1_cast=table2array(data1(2:2:18,4));

col_x1=[177 206 85]/255;
col_x2=[35 155 56]/255;

axes
hold on
x_pos1=[1 1 1 2 2 2 3 3 3]-0.1;
x_pos2=[1 1 1 2 2 2 3 3 3]+0.1;
  plot(x_pos1, data1_b6, '.','color',col_x2,'Markersize',10)
  plot(x_pos2, data1_cast, '.','color',col_x1,'Markersize',10)
  line([0 4],[0 0],'LineStyle','--','Color',[0.5 0.5 0.5])
  set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(2) pos_y(2) graph_size],...
    'xlim',[0.5 3.5],'XTick',[1:3],'XTickLabel',[6 24 48],'ylim',[-7 14],'YTick',[],'Box','on');
print('../../plots/Fig4/Fig4G_exp','-depsc','-loose')

%% Fig. 4I -  H3K27 experiment
clf
colormap([0 0 0; 0.8 0.8 0.8]);

K27=dlmread('../../data/Fig2_3_4/results_96h_IF_FISH.txt');
dox=[2 4 6];
ctl=[1 3 5];
axes

K27_mean=[mean(K27(:,ctl),2) mean(K27(:,dox),2)];
K27_std=[std(K27(:,ctl),[],2) std(K27(:,dox),[],2)];
b=bar(K27_mean(1:3,:));
hold on;
xpos=[b(1).XData+b(1).XOffset;b(2).XData+b(2).XOffset];
    
  plot(xpos(1,:), K27(1:3,ctl), '.','color',[0.5 0.5 0.5],'Markersize',8)
  plot(xpos(2,:), K27(1:3,dox), '.','color',[0.5 0.5 0.5],'Markersize',8)
  errorbar(xpos',K27_mean(1:3,:),K27_std(1:3,:),'k.','Markersize',0.1,'Linewidth',0.7,'CapSize',2.5)
  [tt pval]=ttest2(K27(:,dox)',K27(:,ctl)')

  set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(1) graph_size],...
    'XTickLabel',{'XaXi','Xa*Xi','XiXi'},'ylim',[0 100],'xlim',[0.5 3.5],'XTick',[1:3],'YTick',[0:25:100])
ylabel('Cells [%]','Fontsize',fst);
[le1 le2]=legend({'Control','Dox'},'location','Northwest','Box','off');
le1.Position([1:3])=[le1.Position(1)*1.2 le1.Position(2)*1 le1.Position(3)*0.4];
for q=1:2
    set(le2(q),'Fontsize',fs);
end
print('../../plots/Fig4/Fig4I','-depsc','-loose')

%% Fig 4C - TX dox speedup experiment
clf
colormap([0 0 0; 0.8 0.8 0.8]);

fish_all{1}=load('../../data/Fig2_3_4/TX_fish_a.txt');
fish_all{2}=load('../../data/Fig2_3_4/TX_fish_b.txt');
fish_all{3}=load('../../data/Fig2_3_4/TX_fish_c.txt');

%plot data
pat{2}=[fish_all{1}(:,1) fish_all{2}(:,1) fish_all{3}(:,1)];
pat{1}=[fish_all{1}(:,2) fish_all{2}(:,2) fish_all{3}(:,2)];
tot=[sum(fish_all{1}(:,1:3),2) sum(fish_all{2}(:,1:3),2) sum(fish_all{3}(:,1:3),2)];

for n=1:2
    sel.data=100*pat{n}./tot;
    axes
    
    
plot_mean=100*mean(pat{n}./tot,2);
plot_std=100*std(pat{n}./tot,[],2);
b=bar([plot_mean(6:10) plot_mean(1:5)]);
hold on;
xpos=[b(1).XData+b(1).XOffset;b(2).XData+b(2).XOffset];

plot(xpos(2,:), sel.data(1:5,:), '.','color',[0.5 0.5 0.5],'Markersize',6)
  hold on
  plot(xpos(1,:), sel.data(6:10,:), '.','color',[0.5 0.5 0.5],'Markersize',6)
  errorbar(xpos',[plot_mean(6:10) plot_mean(1:5)],[plot_std(6:10) plot_std(1:5)],'k.','Markersize',0.1,'Linewidth',0.7,'CapSize',2)

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(n) pos_y(4) graph_size],...
    'XTickLabel',[0:4],'ylim',[0 110],'xlim',[0.5 5.5],'YTick',[],'XTick',[1:5])
if n==1
    %ylabel('Cells [%]','Fontsize',fst);
    %xlabel('Differentiation [days]','Fontsize',fst);
    %set(gca,'YTick',[0:25:100]);
end
[tt pval]=ttest2(pat{n}(1:5,:)',pat{n}(6:10,:)')
end
[le1 le2]=legend({'Control','Dox'},'location','Northwest','Box','off');
le1.Position([1:3])=[le1.Position(1)*1 le1.Position(2)*1 le1.Position(3)*0.4];
for q=1:2
    set(le2(q),'Fontsize',fs);
end

print('../../plots/Fig4/Fig4C','-depsc','-loose')

%% Fig. S3D - TX EdU experiment
graph_size=[2.5 2.5];
edu_table=readtable('../../data/Fig2_3_4/edu_data_table.txt')
edu_table_add=readtable('../../data/Fig2_3_4/edu_additional_data_table.txt');

edu_array=table2array(edu_table(:,2:end))+table2array(edu_table_add(:,2:end));
nr_cells_dox=edu_array(1:3,2:2:end)+edu_array(4:6,2:2:end);
nr_cells_dox_sum=nr_cells_dox(:,1:3)+nr_cells_dox(:,4:6)+nr_cells_dox(:,7:9);

edu_pos=100*edu_array(1:3,:)./(edu_array(1:3,:)+edu_array(4:6,:));

xist1=reshape(edu_pos(1,2:2:end),3,3)';
xist1_mean=mean(xist1);
xist1_std=std(xist1);
xist2=reshape(edu_pos(2,2:2:end),3,3)';
xist2_mean=mean(xist2);
xist2_std=std(xist2);
xist0=reshape(edu_pos(3,2:2:end),3,3)';
xist0_mean=mean(xist0);
xist0_std=std(xist0);

clf
b=bar([xist1_mean;xist2_mean]');
hold on
xpos=[b(1).XData+b(1).XOffset;b(2).XData+b(2).XOffset];
plot(xpos(2,:), xist2', '.','color',[0.5 0.5 0.5],'Markersize',8)
hold on
plot(xpos(1,:), xist1' , '.','color',[0.5 0.5 0.5],'Markersize',8)

errorbar(xpos,[xist1_mean;xist2_mean],[xist1_std;xist2_std],'k.','Markersize',0.1,'Linewidth',0.7,'CapSize',2)
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(n) pos_y(4) graph_size],...
    'XTickLabel',[0:4],'ylim',[0 130],'xlim',[0.5 3.5],'YTick',[0:50:100],'XTick',[1:3],'XTickLabel',[6 24 48])
ylabel('EdU^+ Cells [%]','Fontsize',fst);
xlabel('Time after Dox [h]','Fontsize',fst);
plot_nr=nr_cells_dox_sum';
for t=1:6
    tx=text(xpos(t),140,['n=',num2str(plot_nr(t))],'Fontsize',8)
    set(tx,'Rotation',90)
end
set(tx,'Rotation',90)
[le1 le2]=legend({'1x Xist','2x Xist'},'location','Northwest','Box','off');
le1.Position([1:3])=[le1.Position(1)*1.45 le1.Position(2)*1 le1.Position(3)*0.4];
for q=1:2
    set(le2(q),'Fontsize',fs);
end
[tt pval]=ttest2(xist1,xist2)

print('../../plots/FigS3/FigS3D','-depsc','-loose')
