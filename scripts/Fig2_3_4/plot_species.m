% Generate plots in Fig. 3 and S1A
clear
%% Plot settings
graph_size=[2.5 2.5];
lw=1.5;
fs=8;
fst=10;
pos_x=[2:6:16];
pos_y=[2 :3.6:50];

p1=[1900, 100,700,800];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')

%% load simulation data
clf
col_x1=[177 206 85]/255;
col_x2=[35 155 56]/255;
graph_size=[2.5 2.5];
pos_x=[2:3:16];

wt_sum = dlmread('../../simulations/Fig2_3_4/SUMMARY_stoch_sim_Oct2018.txt');
sw = wt_sum(:,998);
par_wt = wt_sum(:,1:33);
ba_wt = wt_sum(:,275:375);
ma_wt = wt_sum(:,34:134);
wt_sets=find(mean(ma_wt(:,81:101),2)>0.8);
sel_wt=wt_sets;
max_ba = 100*max(ba_wt,[],2);
%% Fig. S2A: plot maximal bi-allelic fraction vs switch-on-to-silencing ratio
clf
axes
plot(sw(sel_wt)./par_wt(sel_wt,8),max_ba(sel_wt),'k.','Markersize',1)

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1)+0.5 pos_y(2)+1.5 graph_size],...
    'xlim',[0 30],'ylim',[-5 100],'YTick',[0:25:100],'XTick',[0:10:30])

xlabel('Switch-ON/sil_{tXA}','Fontsize',fst)
ylabel('Bi-allelic [% max]','Fontsize',fst);

print('../../plots/FigS2/FigS2A','-depsc','-loose')

%% Fig. 3d: fit in vivo data
clf
output=[];
%% Fit mouse data
data_mouse = dlmread('../../data/Fig2_3_4/mouse_in_vivo_data_cellnr.txt',' ',1,0);
tp = data_mouse(:,1);
% define range to offset the starting point
offset=[96:120];

% Creat a matrix with all possible combinations of probabilities for BA, MA,
% neg (must add up to 1) in steps of 0.01
p = [0:0.01:1 ; 0:0.01:1];
all_p = [];
for p1=0:0.01:1
    for p2=0:0.01:1
        if (p1+p2)<=1
            all_p = [all_p;p1 p2 1-(p1+p2)];
        end
    end
end

log_lik = zeros(length(wt_sets),length(offset));
nr_sim_cells = 100;

parpool(3);
% C: Neg, MA, BA, cells in data, R: different timepoints 
        data = data_mouse(:,2:end);
parfor f = 1:length(offset)
    f
    tp_new=tp-offset(f);
    tic
    temp2 = zeros(length(wt_sets),1);
    for m = 1:length(wt_sets)
        m
        
        % MA, BA neg cells in simulation
        sim=round(nr_sim_cells*[1-wt_sum(wt_sets(m),tp_new+34)'-wt_sum(wt_sets(m),tp_new+275)' wt_sum(wt_sets(m),tp_new+34)' wt_sum(wt_sets(m),tp_new+275)']);
        
        lik = zeros(size(data,1), size(all_p,1));
        %Sum over all timepoints
        for j = 1:size(data,1)
            % for each probability combination calculate the likelihood that this is
            % the true probability given the data and given the simulation, then
            % multiply
            for n=1:size(all_p,1)
                %This is the likelihood to observe data and simulation
                %given that all_p(n,:) are the true probabilities?
                lik(j,n) = mnpdf(data(j,:),all_p(n,:))*mnpdf(sim(j,:),all_p(n,:));
            end
            % Then sum over all p combinations and take the -log => WHY sum over all p:
            % Integration over all probability distributions???
            % combinations
            temp = -log(sum(lik(j,:)));
            % now sum over all time points (small value is high
            % likelyhood) (log(a)+log(b) = log(a.b) => likelihood to
            % observe data and simulation at all timepoints given the true
            % probability
            temp2(m) = temp2(m)+temp;
            
            %log_lik(m,f) = log_lik(m,f)+temp;
        end
        %clear lik;
    end
    log_lik(:,f) = temp2;
    %clear temp2;
    toc
end

[a b]=min(log_lik(:));
[sel_par d]=ind2sub(size(log_lik),b);
sel_off=offset(d);
%sel_par = 583;
%sel_off = offset(16);
% write best fit to file
out_par=[3:8 11:14 21:23];
output = [sel_off wt_sum(wt_sets(sel_par),out_par)];
% write best 100 fits into file
[min_lik_set ind_off] = min(log_lik, [],2);
[sorted_lik ind_best] = sort(min_lik_set,1, 'ascend');
best_fits = [wt_sum(wt_sets(ind_best(1:100)),:)];
%write offset for each parameter set in column 10
best_fits(:,10) = offset(ind_off(ind_best(1:100)))';
%Write lik of each fit into Column 15
best_fits(:,15) = sorted_lik(1:100); % = min_res_set(ind_best(1:10))
%dlmwrite('../../simulations/Fig2_3_4/par_best_fits_lik_mouse_in_vivo_for_mutant_sim_Oct2018.txt', best_fits(:,1:33));
%dlmwrite('../../simulations/Fig2_3_4/WT_sim_best_fits_lik_mouse_in_vivo_Oct2018.txt', best_fits);

%% Fit rabbit data
clear a; clear b;
data_rabbit = dlmread('../../data/Fig2_3_4/rabbit_in_vivo_data_cellnr.txt',' ',1,0);
tp_rb = data_rabbit(:,1);

% extract the measured time points from the simulations
offset_rb=[43:67];

log_lik_rb = zeros(length(wt_sets),length(offset_rb));
nr_sim_cells = 100;

% C: Neg, MA, BA, cells in data, R: different timepoints 
data_rb = data_rabbit(:,2:end);
parfor f = 1:length(offset_rb)
    tp_new_rb=tp_rb-offset_rb(f);
    tic
    temp2 = zeros(length(wt_sets),1);
    for m = 1:length(wt_sets)
        sim=round(nr_sim_cells*[1-wt_sum(wt_sets(m),tp_new_rb+34)'-wt_sum(wt_sets(m),tp_new_rb+275)' wt_sum(wt_sets(m),tp_new_rb+34)' wt_sum(wt_sets(m),tp_new_rb+275)']);
        lik = zeros(size(data_rb,1), size(all_p,1));
        %Sum over all timepoints
        for j = 1:size(data_rb,1)
            for n=1:size(all_p,1)
                lik(j,n) = mnpdf(data_rb(j,:),all_p(n,:))*mnpdf(sim(j,:),all_p(n,:));
            end
            temp = -log(sum(lik(j,:)));
            temp2(m) = temp2(m)+temp;
        end
    end
    log_lik_rb(:,f) = temp2;
    toc
end

%Plot jth best fit
% [a b]=max(log_lik_rb(:));
% [sel_par2 d]=ind2sub(size(log_lik_rb),b);
% for j = 1:1
%     log_lik_rb(sel_par2,d) = Inf;
%     [a b]=min(log_lik_rb(:));
%     [sel_par2 d]=ind2sub(size(log_lik_rb),b);
% end

[a b]=min(log_lik_rb(:));
[sel_par2 d]=ind2sub(size(log_lik_rb),b);
sel_off2=offset_rb(d);
output = [output; sel_off2 wt_sum(wt_sets(sel_par2),out_par)];

%sel_off2 = 58;
%sel_par2 = 509;
%% plot mouse
axes
b=bar(0:100,100*[ba_wt(wt_sets(sel_par),:);ma_wt(wt_sets(sel_par),:)]','stacked');
b(1).FaceColor=[0.3 0.3 0.3];
b(1).EdgeColor=[0.3 0.3 0.3];
b(2).FaceColor=[0.7 0.7 0.7];
b(2).EdgeColor=[0.7 0.7 0.7];

hold on
plot(tp-sel_off,100*(data_mouse(:,4)./sum(data_mouse(:,2:end),2)),'o','MarkerFaceColor',[0.3 0.3 0.3],'LineWidth',1,'Color',[0.9 0.9 0.9])
plot(tp-sel_off,100*sum(data_mouse(:,3:4),2)./sum(data_mouse(:,2:end),2),'o','MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',1,'Color','k')
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(1) pos_y(2) graph_size],...
        'xlim',[-3 96],'XTick',[],'XTickLabel',[0 2 4],'ytick',[0:50:100],'ylim',[0 100]);
%ylabel('Cells [%]','Fontsize',fst)
%% plot rabbit
axes
b=bar(0:100,100*[ba_wt(wt_sets(sel_par2),:);ma_wt(wt_sets(sel_par2),:)]','stacked')
b(1).FaceColor=[0.3 0.3 0.3];
b(1).EdgeColor=[0.3 0.3 0.3];

b(2).FaceColor=[0.7 0.7 0.7];
b(2).EdgeColor=[0.7 0.7 0.7];

hold on
p1=plot(tp_rb-sel_off2,100*(data_rabbit(:,4)./sum(data_rabbit(:,2:end),2)),'o','MarkerFaceColor',[0.3 0.3 0.3],'LineWidth',1,'Color',[0.9 0.9 0.9]);
p2=plot(tp_rb-sel_off2,100*sum(data_rabbit(:,3:4),2)./sum(data_rabbit(:,2:end),2),'o','MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',1,'Color','k');
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(1) pos_y(1) graph_size],...
        'xlim',[-3 96],'XTick',[0:48:96],'XTickLabel',[0 2 4],'ytick',[0:50:100],'ylim',[0 100]);
xlabel('Time [days]','Fontsize',fst)
print('../../plots/Fig3/Fig3D','-depsc','-loose')

print('../../plots/Fig3/Fig3D','-depsc','-loose')
%% Fit mESC
clf;
clear a; clear b;
%C:1BA,C2:MA,C3:Neg,C4:XO
fish_all{1}=load('../../data/Fig2_3_4/TX_fish_a.txt');
fish_all{2}=load('../../data/Fig2_3_4/TX_fish_b.txt');
fish_all{3}=load('../../data/Fig2_3_4/TX_fish_c.txt');
tp_esc = [24:24:96];
%C1:Neg,C2:MA,C3:BA
data_esc=[fish_all{1}(7:10,[3 2 1])+fish_all{2}(7:10,[3 2 1])+fish_all{3}(7:10,[3 2 1])];

offset_esc=[0:24];
log_lik_esc = zeros(length(wt_sets),length(offset_esc));
parfor f = 1:length(offset_esc)
    tp_new_esc=tp_esc-offset_esc(f);
    tic
    temp2 = zeros(length(wt_sets),1);
    for m = 1:length(wt_sets)
        sim=round(nr_sim_cells*[1-wt_sum(wt_sets(m),tp_new_esc+34)'-wt_sum(wt_sets(m),tp_new_esc+275)' wt_sum(wt_sets(m),tp_new_esc+34)' wt_sum(wt_sets(m),tp_new_esc+275)']);
        lik = zeros(size(data_esc,1), size(all_p,1));
        %Sum over all timepoints
        for j = 1:size(data_esc,1)
            for n=1:size(all_p,1)
                lik(j,n) = mnpdf(data_esc(j,:),all_p(n,:))*mnpdf(sim(j,:),all_p(n,:));
            end
            temp = -log(sum(lik(j,:)));
            temp2(m) = temp2(m)+temp;
        end
    end
    log_lik_esc(:,f) = temp2;
    toc
end

[a b]=min(log_lik_esc(:));
[sel_par_esc d]=ind2sub(size(log_lik_esc),b);
sel_off_esc=offset_esc(d);
output = [output; sel_off_esc wt_sum(wt_sets(sel_par_esc),out_par)];
% write best 100 fits into file
[min_lik_esc ind_off_esc] = min(log_lik_esc, [],2);
[sorted_lik_esc ind_best_esc] = sort(min_lik_esc,1, 'ascend');
best_fits_esc = [wt_sum(wt_sets(ind_best_esc(1:100)),:)];
%write offset for each parameter set in column 10
best_fits_esc(:,10) = offset_esc(ind_off_esc(ind_best_esc(1:100)))';
%dlmwrite('../../simulations/Fig2_3_4/par_best_fits_lik_mESC_for_mutant_sim_Oct2018.txt', best_fits_esc(:,1:33));
%dlmwrite('../../simulations/Fig2_3_4/WT_sim_best_fits_lik_mESC_Oct2018.txt', best_fits_esc);

% plot mESC
axes
b=bar(0:100,100*[wt_sum(wt_sets(sel_par_esc),275:375);wt_sum(wt_sets(sel_par_esc),34:134)]','stacked');
b(1).FaceColor=[0.3 0.3 0.3];
b(1).EdgeColor=[0.3 0.3 0.3];
b(2).FaceColor=[0.7 0.7 0.7];
b(2).EdgeColor=[0.7 0.7 0.7];

hold on
plot(tp_esc-sel_off_esc,100*(data_esc(:,3)./sum(data_esc,2)),'o','MarkerFaceColor',[0.3 0.3 0.3],'LineWidth',1,'Color',[0.9 0.9 0.9])
plot(tp_esc-sel_off_esc,100*sum(data_esc(:,2:3),2)./sum(data_esc,2),'o','MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',1,'Color','k')
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(1) pos_y(3) graph_size],...
        'xlim',[-3 100],'XTick',[0:48:96],'XTickLabel',[0 2 4],'ytick',[0:50:100],'ylim',[0 100]);
ylabel('Cells [%]','Fontsize',fst)
xlabel('Time [days]','Fontsize',fst)

print('../../plots/FigS1/FigS1A','-depsc','-loose')
dlmwrite('../../simulations/Fig2_3_4/best_fit_pars_lik.txt',output,'\t');

%% Load simulation data on human
clear data;
x1=34:234;
x2=235:435;
filenames={'../../simulations/Fig2_3_4/stoch_sim_damp_XR_upreg_Oct2018.txt','../../simulations/Fig2_3_4/stoch_sim_damp_XR_damp_Oct2018.txt'};
for z=1:length(filenames)
    data{z} = dlmread(filenames{z});
    ia = 1:100:size(data{z},1);
    nr_cells = size(data{z},1)./length(ia);
    xi = zeros(size(data{z}));
    xi(data{z}>=(data{z}(:,2).*data{z}(:,21))./5) = 1;
    
    nr_xi{z} = [xi(:,x1)+xi(:,x2)];
    ma{z} = zeros(size(nr_xi{z},1)./nr_cells, size(nr_xi{z},2));
    ba{z} = zeros(size(nr_xi{z},1)./nr_cells, size(nr_xi{z},2));

    for i = 1:length(ia)
        temp=nr_xi{z}(ia(i):ia(i)+99,:);
        ma{z}(i,:) = sum(temp==1);
        ba{z}(i,:) = sum(temp==2);
    end
end

%% Fig. 3e: plot human simulations
figure(1)
clf
figure(2)
clf
col_x1=[177 206 85]/255;
col_x2=[35 155 56]/255;

pos_x=[2:3:16];
pos_y=[2 :3:50];

for z=1:2
figure(1)
n=[2 17];
q=[1 1];
chr1=data{z}(ia(n(z))+q,x1);
chr2=data{z}(ia(n(z))+q,x2);
axes
plot(0:200,chr1,'LineWidth',lw,'Color',col_x1)
hold on
plot(0:200,chr2,'LineWidth',lw,'Color',col_x2)
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(z) pos_y(2) graph_size],...
    'xlim',[-3 200],'XTick',[],'XTickLabel',[0:2:8],'ytick',[],'ylim',[-20 1.1*max([chr1(:); chr2(:)])]);
if z==1
   ylabel('Xist [# molecules]','Fontsize',fst)
   set(gca,'ytick',[0:250:1000]);
else
[le le2]=legend({'Xist1','Xist2'},'location','Northeast','Box','off');
le.Position=[le.Position(1)*1.5 le.Position(2)*0.9 le.Position(3)*0.3 le.Position(4)*1.5];
%title('100 cells','Fontsize',fst)
for q=1:2
    set(le2(q),'Fontsize',fs);
end
end
axes

b=bar(0:200,[ba{z}(n(z),:);ma{z}(n(z),:)]','stacked');
b(1).FaceColor=[0.3 0.3 0.3];
b(1).EdgeColor=[0.3 0.3 0.3];
b(2).FaceColor=[0.7 0.7 0.7];
b(2).EdgeColor=[0.7 0.7 0.7];

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(z) pos_y(1) graph_size],...
        'xlim',[-3 200],'XTick',[0:48:200],'XTickLabel',[0:2:8],'ytick',[],'ylim',[0 110]);
if z==1
    ylabel('Cells [%]','Fontsize',fst)
    set(gca,'ytick',[0:50:100])
else
xlabel('Time [days]','Fontsize',fst)
[le le2]=legend({'bi-allelic','mono-allelic'},'location','Northeast','Box','off');
le.Position=[le.Position(1)*1.8 le.Position(2)*0.9 le.Position(3)*0.3 le.Position(4)*1.5];
%title('100 cells','Fontsize',fst)
for q=1:2
    set(le2(q),'Fontsize',fs);
end
end
figure(2)
axes
all_dat=[ma{z}(:,1:48:200);ba{z}(:,1:48:200)];
day=repmat(0:4,size(ma{z},1)+size(ba{z},1),1);
pat=[repmat(1,size(ma{z},1),5); repmat(2,size(ba{z},1),5)];

bp=boxplot(all_dat(:),{day(:), pat(:)},'factorgap',10,'color',[0.7 0.7 0.7;0 0 0],'Symbol','.','OutlierSize',2);
set(bp(:,:),'LineWidth',1);
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(z) pos_y(2) graph_size],...
    'ylim',[-10 110],'xtick',1.5:2.9:50,'xticklabel',0:2:8,'ytick',[],'xlim',[0 15]);
if z==1
ylabel('Cells [%]','Fontsize',fst);
set(gca,'ytick',[0:50:100])
xlabel('Time [days]','Fontsize',fst);

end
end
figure(1)
print('../../plots/Fig3/Fig3E','-depsc','-loose')
figure(2)
print('../../plots/FigS2/FigS2B','-depsc','-loose')
