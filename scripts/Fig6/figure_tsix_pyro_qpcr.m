pyro=dlmread('../../data/Fig6/pyro_tsix_matlab.txt');

figure(1)
%% plotting parameters
graph_size=[2.5 2.5];
lw=1.5;
fs=8;
fst=10;
pos_x=[2 5 8 8];
pos_x=[2:5:22 26];
pos_x=[2 5 10 13];

pos_y=[2 :4.6:15.8 19];
p1=[1800, 1400,700,800];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')
%% plot Pyro data
figure(1)
clf
time=[0 0.5 1 2 4 8];
temp=0.01*reshape(pyro(:,1),6,3);
f0=mean(temp(1,:));
f0=repmat(temp(1,:),6,1);
up=temp.*(1-f0)./(f0.*(1-temp));

%up=0.01*temp./mean(temp(1,:));
temp=0.01*reshape(pyro(:,2),6,3);
f0=repmat(temp(1,:),6,1);
down=temp.*(1-f0)./(f0.*(1-temp));

%down=temp./mean(temp(1,:));
col1=[14 98 130]/255;
col2=[161 194 205]/255;
col2=col2*0.8;

errorbar(time',mean(log2(down),2),std(log2(down),[],2),'.-','Linewidth',lw,'Markersize',16,'Color',col2);
hold on
errorbar(time',mean(log2(up),2),std(log2(up),[],2),'.-','Linewidth',lw,'Markersize',16,'Color',col1);

set(gca,'YTick',[-2:1],'XTick',[0:2:8],'XTickLabel',[],'TickLength',[0.02 0],'xlim',[-0.5 8.5],'ylim',[-2 0.5],...
'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(2) pos_y(3) graph_size]);
%xlabel('Time after Dox [h]','Fontsize',fst);
ylabel({'Rel. Tsix [log2]'},'Fontsize',fst);
title('TX XX (Pyro)','Fontsize',fst)
[le le2]=legend({'Tsix 3''','Tsix 5'''},'Location','SouthEast','Fontsize',fst,'Box','off')
le.Position=[1.25*le.Position(1) 1.06*le.Position(2) le.Position(3)*0.5 le.Position(4)];
for q=1:2
    set(le2(q),'Fontsize',fs);
end

[a p]=ttest2(up',down')
    
%% plot  qPCR
axes
data=load('../../data/Fig6/130703_matlab.txt');
ref=[1:3];
data(2:2:end,:)=[];
data_norm=(repmat(geomean(data(:,ref),2),1,size(data,2))-data);
data_norm(17:18,:)=NaN; % samples degraded

genenames={'Arpo','Gapdh','rrm2','spliced Tsix','Tsix 3','Xist','Tsix 5','Rnf12','Jpx','Tsx'};
gene=[4 5 7];
clear mean_data std_data


for g1=1:length(gene)
    temp=(reshape(data_norm(:,gene(g1)),12,3));
    %Only XO
    rel_data{g1}=temp(1:6,:)-repmat(temp(1,:),6,1);
    mean_data(:,g1)=nanmean(temp(1:6,:),2);
    te(:,g1)=ttest2(temp(4,:),temp(5,:))
    std_data(:,g1)=nanstd(temp(1:6,:),[],2);
end

mean_plot=mean_data(1:6,:)-repmat(mean_data(1,:),6,1);
std_plot=std_data(1:6,:);
set(0,'defaultAxesColorOrder',[[145 38 143]/255; col2; col1]);

errorbar(repmat(time',1,3),mean_plot,std_plot,'.-','Linewidth',lw,'Markersize',16)

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(2) pos_y(2) graph_size],...
    'XTickLabel',[0:2:8],'ylim',[-1.6 0.5],'xlim',[-0.5 8.5],'yticklabel',[-1.5:0.5:1.5])
ylabel('Rel. Tsix [log2]','Fontsize',fst);
xlabel('Time after Dox [h]','Fontsize',fst);
title('TX XO (qPCR)','Fontsize',fst)
[le le2]=legend({'spliced Tsix','Tsix 3''','Tsix 5''','spliced Tsix','bg'},'Location','SouthEast','Fontsize',fst,'Box','off')
le.Position=[1.33*le.Position(1) 1.04*le.Position(2) le.Position(3)*0.5 le.Position(4)];
for q=1:2
    set(le2(q),'Fontsize',fs);
end

spliced_3=0
[a1 p1]=ttest2(rel_data{1}',rel_data{2}')
spliced_5=0
[a2 p2]=ttest2(rel_data{1}',rel_data{3}')
d3_5=0
[a3 p3]=ttest2(rel_data{2}',rel_data{3}')

print('../../plots/Fig6/Fig6C','-depsc','-loose')






