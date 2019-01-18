%% Generate simulation plots in Fig. 4
%plot settings
clear
graph_size=[2.5 2.5];
lw=1.5;
fs=8;
fst=10;
pos_x=[2:2.8:20];
pos_y=[2:4.5:20];
figure(1)
clf
p1=[1800, 100,700,700];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')
%% load simulated data

filenames={'../../simulations/Fig2_3_4/stoch_sim_WT_mESC_fits.txt','../../simulations/Fig2_3_4/stoch_sim_Dox_ind_Oct2018.txt',...
    '../../simulations/Fig2_3_4/stoch_sim_BA_exp_Oct2018.txt','../../simulations/Fig2_3_4/stoch_sim_BA_exp_Oct2018.txt'}; 

for z=1:length(filenames)
    data = dlmread(filenames{z});
    ia = 1:100:size(data,1);
    par{z}= data(ia,1:31);
    nr_cells = size(data,1)./length(ia);
    xi = zeros(size(data));
    xi(data>=(data(:,2).*data(:,21))./5) = 1;
    if z==1||z==3
        x1=34:134;
        x2=135:235;
    elseif z==2
        x1=58:158;
        x2=183:283;
    elseif z==4
        x1=236:336;
        x2=337:437;
    end
    nr_xi{z} = [xi(:,x1)+xi(:,x2)];
    ma{z} = zeros(size(nr_xi{z},1)./nr_cells, size(nr_xi{z},2));
    ba{z} = zeros(size(nr_xi{z},1)./nr_cells, size(nr_xi{z},2));
    for i = 1:length(ia)
        temp=nr_xi{z}(ia(i):ia(i)+99,:);
        ma{z}(i,:) = sum(temp==1);
        ba{z}(i,:) = sum(temp==2);
    end
end

%% Fig. 4B: dox speedup
clf
n=1;
colormap([0 0 0;0.3 0.3 0.3;0.7 0.7 0.7]);

for p=1:2
    if p==1
        dat=ma;
    else
        dat=ba;
    end
    axes
    bar(0:4,[dat{1}(n,1:24:100+par{1}(n,10));dat{2}(n,1:24:100+par{2}(n,10))]')

    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(p) pos_y(1) graph_size],...
    'XTickLabel',[0:4],'ylim',[0 110],'xlim',[-0.5 4.5],'YTick',[])

    if p==1
        set(gca,'YTick',[0:25:100]);
    end
end
[le1 le2]=legend({'Control','Dox'},'location','Northwest','Box','off');
le1.Position([1:3])=[le1.Position(1)*1 le1.Position(2)*1 le1.Position(3)*0.4];

print('../../plots/Fig4/Fig4B','-depsc','-loose')

%% Fig. S3A - plot boxplots
clf
tp = false(100,101);
for u=1:size(tp,1)
    tp(u,1:24:100+par{1}(u,10)) = 1;
end
%tp = 
%tp = repmat(1:24:100,100,1) + par{1}(:,10);
all_ma=[reshape(ma{1}(tp),100,5);reshape(ma{2}(tp),100,5)];
all_ba=[reshape(ba{1}(tp),100,5);reshape(ba{2}(tp),100,5)];
day=repmat(0:4,size(ma{1},1)+size(ma{2},1),1);
pat=[repmat(1,size(ma{1},1),5); repmat(2,size(ma{2},1),5)];
for p=1:2
    if p==1
        dat=all_ma;
    else
        dat=all_ba;
    end
    axes
    bp=boxplot(dat(:),{day(:), pat(:)},'factorgap',10,'color',[0 0 0;0.7 0.7 0.7],'Symbol','.','OutlierSize',4);
    set(bp(:,:),'LineWidth',1);
    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(p) pos_y(1) graph_size],...
    'ylim',[-10 110],'xtick',1.5:2.9:50,'xticklabel',0:4,'ytick',[],'xlim',[0 15]);
    if p==1
        set(gca,'ytick',0:50:100)
        ylabel('Cells [%]','Fontsize',fst);
        %xlabel('Differentiation [days]','Fontsize',fst);
    end
end

print('../../plots/FigS3/FigS3A','-depsc','-loose')

%% Fig. 4G - plot BA induction simulation
figure(1)
clf
col_x1=[177 206 85]/255;
col_x2=[35 155 56]/255;

n=1;
d1=34:134;
d2=135:235;
c1=236:336;
c2=337:437;
data=dlmread('../../simulations/Fig2_3_4/stoch_sim_BA_exp_Oct2018.txt');
ia = 1:100:size(data,1);
nr_cells = size(data,1)./length(ia);
mean_xist = mean(data(ia(n):ia(n)+100,:));
ti=[56 72 96];
plot(1:3,log2(mean_xist(d1(ti+1))./mean_xist(c1(ti+1))),'o','Markersize',5,'MarkerFaceColor',col_x2,'Color',col_x2)
hold on
plot(1:3,log2(mean_xist(d2(ti+1))./mean_xist(c2(ti+1))),'o','Markersize',5,'MarkerFaceColor',col_x1,'Color',col_x1)

line([0 100],[0 0],'LineStyle','--','Color',[0.5 0.5 0.5])
%
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(1) graph_size],...
    'xlim',[0.5 3.5],'XTick',[1:3],'XTickLabel',[6 24 48],'ylim',[-7 14],'YTick',[-5:5:10],'Box','on');
%xlabel('Time after Dox [h]','Fontsize',fst)
ylabel(['Xist in', 10,'Dox vs. Ctl [log2]'],'Fontsize',fst)

print('../../plots/Fig4/Fig4G_sim','-depsc','-loose')
%% Fig. S3B - boxplot
clf
for c=1:length(ia)
    mean_xist = mean(data(ia(c):ia(c)+99,:));
    to_plot(1:3,c)=mean_xist(d1(ti+1))./mean_xist(c1(ti+1));
    to_plot(4:6,c)=mean_xist(d2(ti+1))./mean_xist(c2(ti+1));
end
    day=repmat([6 24 48]',2,size(to_plot,2));
    pat=[repmat(1,3,size(to_plot,2)); repmat(2,3,size(to_plot,2))];
    
    bp=boxplot(log2(to_plot(:)),{day(:), pat(:)},'factorgap',15,'color',[col_x2;col_x1],'Symbol','.','OutlierSize',4)
    line([0 100],[0 0],'LineStyle','--','Color',[0.5 0.5 0.5])
    set(bp(:,:),'LineWidth',1);
    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(1) graph_size],...
    'ylim',[-9 17],'xtick',1.5:2.9:50,'xticklabel',[6 24 48],'ytick',[-10:10:10],'xlim',[0 9]);
    ylabel(['Xist in', 10,'Dox vs. Ctl [log2]'],'Fontsize',fst);
    xlabel('Time after Dox [h]','Fontsize',fst)
 
print('../../plots/FigS3/FigS3B','-depsc','-loose')
    



