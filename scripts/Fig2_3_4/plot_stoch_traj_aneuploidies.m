clear
%% Generate plots in Fig. 2
%plot settings
graph_size=[2.5 2.5];
graph_size2=[1.5,2.5];
lw=1.5;
fs=8;
fst=10;
pos_x=[2:3:20];
pos_x2=[2:2:20];
pos_y=[2:4.5:20];
figure(1)
clf
p1=[1800, 100,700,700];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')
%% load simulated data

wt = dlmread('../../simulations/Fig2_3_4/MA_sets_stoch_sim_Oct2018.txt');
ia = 1:100:size(wt,1);
nr_cells = size(wt,1)./length(ia);
xi = zeros(size(wt));
xi(wt>=(wt(:,2).*wt(:,21))./5) = 1;

nr_xi_wt = [xi(:,34:134)+xi(:,275:375)];

ma_wt = zeros(size(nr_xi_wt,1)./nr_cells, size(nr_xi_wt,2));
ba_wt = zeros(size(nr_xi_wt,1)./nr_cells, size(nr_xi_wt,2));
rat_wt = NaN(size(nr_xi_wt,1)./nr_cells, size(nr_xi_wt,2));
ma_sets_wt = [];

for i = 1:length(ia)
    temp=nr_xi_wt(ia(i):ia(i)+99,:);
    ma_wt(i,:) = sum(temp==1);
    ba_wt(i,:) = sum(temp==2);
    par_wt(i,:)=wt(ia(i),1:33);
end
wt_sets=mean(ma_wt(:,81:101),2)>0.8;
max_ba = max(ba_wt,[],2);
%% Fig. 2C+D - plot trajectories
clf
col_x1=[177 206 85]/255;
col_x2=[35 155 56]/255;
    
ti={'Cell 1','Cell 2','Cell 3'};

n=56;
q=5;
x1=wt(ia(n)+[1:3]+q,34:134);
x2=wt(ia(n)+[1:3]+q,275:375);
for z=1:3
    axes
    plot(0:1:100,wt(ia(n)+z+q,34:134),'LineWidth',lw,'Color',col_x1)
    hold on
    plot(0:1:100,wt(ia(n)+z+q,275:375),'LineWidth',lw,'Color',col_x2)
    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x2(z) 2 graph_size2],...
        'xlim',[-3 96],'XTick',[0 48 96],'XTickLabel',[0 2 4],'ytick',[],'ylim',[0 1.1*max([x1(:); x2(:)])]);
    if z==1
        set(gca,'ytick',[0:100:500])
        ylabel('Xist [# molecules]','Fontsize',fst)
        xlabel('Time [days]','Fontsize',fst)
    end
    title(ti{z},'Fontsize',fst)
end
%
axes

b=bar(0:100,[ba_wt(n,:);ma_wt(n,:)]','stacked')
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

print('../../plots/Fig2/Fig2C_D','-depsc','-loose')

%% Fig. 2H+I +S1B - aneuploidies

cm=flipud([255 255 204;194 230 153;120 198 121;35 132 67;0 0 0]./255);
figure(1)
clf
colormap(cm)
figure(2)
clf
colormap(cm)
nr_par_sets=50;
graph_size=[2.5 2.5];
pos_y=[2 6.5 11 15.5 19.5:3.8:30];
ti={'XO/XY (2n)','XX (2n)','XXX (2n)','XXXX (2n)','XXX (3n)','XXXX (4n)','XXX (3n)','XXXX (4n)'};

filenames={'../../simulations/Fig2_3_4/stoch_sim_2n1X_Oct2018.txt','../../simulations/Fig2_3_4/stoch_sim_WT_mouse_in_vivo_fits.txt',...
    '../../simulations/Fig2_3_4/stoch_sim_2n3X_Oct2018.txt','../../simulations/Fig2_3_4/stoch_sim_2n4X_Oct2018.txt',...
    '../../simulations/Fig2_3_4/stoch_sim_3n3X_XA_rep_Oct2018.txt','../../simulations/Fig2_3_4/stoch_sim_4n4X_XA_rep_Oct2018.txt',...
    '../../simulations/Fig2_3_4/stoch_sim_3n3X_XA_dil_Oct2018.txt','../../simulations/Fig2_3_4/stoch_sim_4n4X_XA_dil_Oct2018.txt'}; 

nr_par_sets=50;
for z=1:length(filenames)
    temp=dlmread(filenames{z});
    sim=zeros(size(temp,1),437);
    sim(:,1:size(temp,2))=temp;
    xi = zeros(size(sim));
    xi(sim>=(sim(:,2).*sim(:,21))./5) = 1;

    nr_xi = reshape(sum(xi(:,[134,235,336,437]),2),100,100);
    to_plot=zeros(100,5);
    for w=1:5
        to_plot(:,w)=sum(nr_xi==(w-1))';  
    end
    to_plot=fliplr(to_plot);
    if z<=6
        figure(1)
    else
        figure(2)
    end
    axes
    [~, b]=max(mean(to_plot));
    bar(sortrows(to_plot(1:nr_par_sets,:),b),'stacked','EdgeColor','none')
    if z<=4
        ps=[pos_x2(z) pos_y(2)];
    elseif z<=6
        ps=[pos_x2(z-2) pos_y(1)];
    else
        ps=[pos_x2(z-4) pos_y(1)];
    end
    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters',...
    'Position',[ps graph_size2],...
    'ylim',[0 100],'YTick',[],'xlim',[0 nr_par_sets],'XTick',[]);
    title(ti{z},'Fontsize',fst)
    if ismember(z,[1,5,7])
        set(gca,'YTick',[0 50 100]);
        ylabel('Cells [%]','Fontsize',fst)
        xlabel('parameter sets','Fontsize',fst)
    end
end

figure(1)
print('../../plots/Fig2/Fig2H_I','-depsc','-loose')

figure(2)
print('../../plots/Fig2/Fig2J','-depsc','-loose')

