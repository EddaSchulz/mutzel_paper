%plot settings
graph_size=[2.5 2.5];
lw=1.5;
fs=8;
fst=10;
pos_x=[2:3:20];
pos_y=[2:4.5:20];
figure(1)
clf
p1=[1800, 100,700,700];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')
%% load simulation data
sim=dlmread('../../simulations/Fig5_7/sim_2C_Xist_Tsix_8kb_indiv_cells_dir_Oct2018.txt');
sim_all=dlmread('../../simulations/Fig5_7/sim_2C_Xist_Tsix_8kb_Oct2018.txt');
[a b1]=ismember(sim(1,1:15),sim_all(:,1:15),'rows');
sim_wt=dlmread('../../simulations/Fig5_7/sim_2C_model_no_X_rep_totXA2_bs_sets_rand_log_k12_13_14_15_complete_170817.txt');
%% plot trajectories for parameter set

clf
col_x1=[177 206 85]/255;
col_x2=[35 155 56]/255;

% plot 3 individual cells
ti={'Cell 1','Cell 2','Cell 3'};

n=56;
q=1;
x1=sim(q+[1:3]+q,26:126);
x2=sim(q+[1:3]+q,137:237);
for z=1:3
    axes
    plot(0:1:100,sim(q+z,26:126),'LineWidth',lw,'Color',col_x1)
    hold on
    plot(0:1:100,sim(q+z,137:237),'LineWidth',lw,'Color',col_x2)
    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(z) 2 graph_size],...
        'xlim',[-3 96],'XTick',[0 48 96],'XTickLabel',[0 2 4],'ytick',[],'ylim',[0 1.1*max([x1(:); x2(:)])]);
    if z==1
        set(gca,'ytick',[0:250:1000])
        ylabel('Xist [# molecules]','Fontsize',fst)
        xlabel('Time [days]','Fontsize',fst)
    end
    title(ti{z},'Fontsize',fst)
end

% plot fraction in 100 cells
axes

b=bar(0:100,100*[sim_all(b1,138:238);sim_all(b1,27:127)]','stacked')
b(1).FaceColor=[0.3 0.3 0.3];
b(1).EdgeColor=[0.3 0.3 0.3];
b(2).FaceColor=[0.7 0.7 0.7];
b(2).EdgeColor=[0.7 0.7 0.7];

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(4)+1 2 graph_size],...
        'xlim',[-3 96],'XTick',[0 48 96],'XTickLabel',[0 2 4],'ytick',[0:50:100],'ylim',[0 100]);
 ylabel('Cells [%]','Fontsize',fst)
xlabel('Time [days]','Fontsize',fst)
[le le2]=legend({'bi-allelic','mono-allelic'},'location','Northeast','Box','off');
le.Position=[le.Position(1)*1.31 le.Position(2)*0.9 le.Position(3)*0.3 le.Position(4)*1.5];
title('100 cells','Fontsize',fst)
for q=1:2
    set(le2(q),'Fontsize',fs);
end

% plot distribution of mono-allelic up-regulation across all parameter sets
% for the human model
axes
temp=histogram(100*mean(sim_all(:,108:127),2),10);
a=temp.Values;
b=temp.BinEdges(2:end)-0.5*temp.BinWidth;
yl=[0 5];
a_plot=100*a/sum(a);
a_plot(1)=a_plot(1)*0.01*yl(2);

bar(b,a_plot,1,'FaceColor',[0.5 0.5 0.5])

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(4)+1 pos_y(2) graph_size],...
    'xlim',[0 100],'ylim',yl,'XTick',[0:20:100],'YTick',[0:5],'YTickLabel',[0:4 100]);

xlabel('XaXi [% cells]','Fontsize',fst)
ylabel('Parameter Sets [%]','Fontsize',fst)

% plot distribution of mono-allelic up-regulation across all parameter sets
% for the mouse model
axes
temp=histogram(100*mean(sim_wt(:,108:127),2),10);
a=temp.Values;
b=temp.BinEdges(2:end)-0.5*temp.BinWidth;
yl=[0 5];
a_plot=100*a/sum(a);
a_plot(1)=a_plot(1)*0.01*yl(2);

bar(b,a_plot,1,'FaceColor',[0.5 0.5 0.5])

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(4)+1 pos_y(3) graph_size],...
    'xlim',[0 100],'ylim',yl,'XTick',[0:20:100],'YTick',[0:5],'YTickLabel',[0:4 100]);

xlabel('XaXi [% cells]','Fontsize',fst)
ylabel('Parameter Sets [%]','Fontsize',fst)



print('../../plots/FigS6/FigS6','-depsc','-loose')
