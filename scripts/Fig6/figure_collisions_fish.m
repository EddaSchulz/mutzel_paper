clear
clf
% plotting parameters
graph_size=[2.5 2.5];
lw=1.5;
fs=8;
fst=10;
pos_x=[2:3:20];

pos_y=[2  5];
p1=[800, 1400,700,800];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')

%% load data

mult_factor = 5; % cutoff for positive/negative signals = # of stdev above background mean
no_outliers_bg = 3; % delete the first n outliers to avoid contamination from signal
filename{1} = ('../../data/Fig6/quant_TXY_Rep1_all.txt');
filename{2} = ('../../data/Fig6/quant_TXY_Rep2_all.txt');
filename{3} = ('../../data/Fig6/quant_TXY_Rep3_all.txt');
filename{4} = ('../../data/Fig6/quant_TXYdA_Rep1_all.txt');
filename{5} = ('../../data/Fig6/quant_TXYdA_Rep2_all.txt');
filename{6} = ('../../data/Fig6/quant_TXYdA_Rep3_all.txt');
filename{7} = ('../../data/Fig6/quant_TXY_Rep1_bgd.txt');
filename{8} = ('../../data/Fig6/quant_TXY_Rep2_bgd.txt');
filename{9} = ('../../data/Fig6/quant_TXY_Rep3_bgd.txt');
filename{10} = ('../../data/Fig6/quant_TXYdA_Rep1_bgd.txt');
filename{11} = ('../../data/Fig6/quant_TXYdA_Rep2_bgd.txt');
filename{12} = ('../../data/Fig6/quant_TXYdA_Rep3_bgd.txt');


for f=1:length(filename)
    temp = dlmread(filename{f});
    data{f} = reshape(temp(:,2),3,size(temp,1)/3)';
    
end

wt{1} = data{1};
wt{2} = data{2};
wt{3} = data{3};
wt{4}=[data{1};data{2};data{3}];
mut{1} = data{4};
mut{2} = data{5};
mut{3} = data{6};
mut{4}=[data{4};data{5};data{6}];

%% calculate the detection threshold using background measurements
bg_all = [data{7};data{8};data{9};data{10};data{11};data{12}];
bg_wt{1} = data{7}; 
bg_mut{1} = data{10}; 
bg_wt{2} = data{8}; 
bg_mut{2} = data{11}; 
bg_wt{3} = data{9}; 
bg_mut{3} = data{12}; 


% remove top 1% of cells from background distribution
[a b] = sort(mean(bg_all,2),'descend');
bg_all(b(1:round(0.01*size(bg_all,1))),:)=[];
thresh = mean(bg_all)+mult_factor*std(bg_all);

for i = 1:3
    [a b] = sort(mean(bg_wt{i},2),'descend');
    bg_wt{i}(b(1:round(0.01*size(bg_wt{i},1))),:)=[];
    thresh_wt{i} = mean(bg_wt{i})+mult_factor*std(bg_wt{i});
    [a b] = sort(mean(bg_mut{i},2),'descend');
    bg_mut{i}(b(1:round(0.01*size(bg_mut{i},1))),:)=[];
    thresh_mut{i} = mean(bg_mut{i})+mult_factor*std(bg_mut{i});
    temp = [bg_wt{i};bg_mut{i}];
    [a b] = sort(mean(temp,2),'descend');
    temp(b(1:round(0.01*size(temp,1))),:)=[];
    thresh_both{i} = mean(temp)+mult_factor*std(temp);
end

%% scatter plot

yl_xist=[0 0.7]*max([wt{1}(:,1);mut{1}(:,1);wt{2}(:,1);mut{2}(:,1);wt{3}(:,1);mut{3}(:,1)]);
yl_tsix(1,:)=[-0.05 1.1]*max([wt{1}(:,2);mut{1}(:,2);wt{2}(:,2);mut{2}(:,2);wt{3}(:,2);mut{3}(:,2);]);
yl_tsix(2,:)=[-0.05 1.1]*max([wt{1}(:,3);mut{1}(:,3);wt{2}(:,3);mut{2}(:,3);wt{3}(:,3);mut{3}(:,3)]);
    

for i = 1:4
    figure(i);
    clf
    y_lab={'Tsix 5'' [a.u.]','Tsix 3'' [a.u.]'}; 
    for m=[1 2]
        for n=[1 2]
            if m==1
                data_x=wt{i}(:,1);
                data_y=wt{i}(:,n+1);
            else
                data_x=mut{i}(:,1);
                data_y=mut{i}(:,n+1);
            end
            axes
            plot(data_x,data_y,'k.','Markersize',1)
            hold on
            if i==4
                plot(yl_xist,[thresh(n+1) thresh(n+1)],'Linewidth',0.5,'Color',[0.5 0.5 0.5])
                plot([thresh(1) thresh(1)],yl_tsix(n,:),'Linewidth',0.5,'Color',[0.5 0.5 0.5])
            else
                if m==1
                    plot(yl_xist,[thresh_wt{i}(n+1) thresh_wt{i}(n+1)],'Linewidth',0.5,'Color',[0.5 0.5 0.5])
                    plot([thresh_wt{i}(1) thresh_wt{i}(1)],yl_tsix(n,:),'Linewidth',0.5,'Color',[0.5 0.5 0.5])
                else
                    plot(yl_xist,[thresh_mut{i}(n+1) thresh_mut{i}(n+1)],'Linewidth',0.5,'Color',[0.5 0.5 0.5])
                    plot([thresh_mut{i}(1) thresh_mut{i}(1)],yl_tsix(n,:),'Linewidth',0.5,'Color',[0.5 0.5 0.5])
                end
            end
            set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(m) pos_y(abs(n-3)) graph_size],...
                'ylim',yl_tsix(n,:),'xlim',yl_xist,'XTick',[0:5:15]*10000)
            if m==2
                set(gca,'yTick',[])
            else
                ylabel(y_lab{n},'Fontsize',fst);
            end
            if n==1
                set(gca,'XTick',[]);
            else
                xlabel('Xist [a.u.]','Fontsize',fst);
            end      
        end
    end
end
%% plot boxplots


for i = 1:4
    data_all{i}=[wt{i}; mut{i}];
    geno{i} = [ones(size(wt{i},1),1); 2*ones(size(mut{i},1),1)];
    xist{i} = zeros(size(data_all{i},1),1);
    if i==4
        xist{i}(data_all{i}(:,1)>thresh(1)) = 1;
    else
        xist{i}(data_all{i}(:,1)>max([thresh_wt{i}(1), thresh_mut{i}(1)])) = 1;
    end
end
%Number of Xist positive cells in WT Replicate 1
%size(find(geno{1}==1&xist{1}==1))
%Number of Xist positive cells in MUT Replicate 1
%size(find(geno{1}==2&xist{1}==1))
%Quantification of Tsix5' in all Xist+ in WT Rep1
%data_all{1}(geno{1}==1&xist{1}==1,2)

ti={'Tsix 5'' [log10]','Tsix 3'' [log10]'};
xist_color=[35 155 56]/255;
for i = 1:4
    figure(i);
    for n=1:2
       axes
        bp=boxplot(log10(data_all{i}(:,n+1)),{geno{i}, xist{i}},'factorgap',30,'color',[0 0 0;xist_color],'Symbol','.');
        hold on 
        if i==4
            plot([0 6],log10([thresh(n+1) thresh(n+1)]),':','Linewidth',1,'Color',[0.5 0.5 0.5])
        else
            plot([0 6],log10([thresh_both{i}(n+1) thresh_both{i}(n+1)]),':','Linewidth',1,'Color',[0.5 0.5 0.5])
        end
        set(bp(:,:),'LineWidth',1);
        h=findobj(gca,'tag','Upper Whisker');
        set(h,'LineStyle','-')
        h=findobj(gca,'tag','Lower Whisker');
        set(h,'LineStyle','-')
        set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters',...
            'Position',[pos_x(4) pos_y(abs(n-3)) graph_size],'xlim',[0 6],'XTick',[1.5 4.5],'XTickLabel',{'TXY','TXY\DeltaA'});
        if n==1
            set(gca,'XTickLabel',[])
        end
        if n==2
            set(gca,'ylim',[4 6])
        end
        ylabel(ti{n},'Fontsize',fst);
    end

end

figure(1);
print('../../plots/FigS5/FigS5A_D','-depsc','-loose')
figure(2);
print('../../plots/FigS5/FigS5B_E','-depsc','-loose')
figure(3);
print('../../plots/FigS5/FigS5C_F','-depsc','-loose')
figure(4);
print('../../plots/Fig6/Fig6E_F','-depsc','-loose')
