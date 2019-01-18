clear
% load results from ODE simulation
sim = dlmread('../../simulations/Fig2_3_4/tXA_cXR_model_scan_for_stoch_sim.txt');
mono=(sim(:,34)./sim(:,35))<0.1|(sim(:,34)./sim(:,35))>10;
no_bi=(sim(:,38)<0.1*max(sim(:,34),sim(:,35)))|(sim(:,39)<0.1*max(sim(:,34),sim(:,35)));
male_off = sim(:,36)<0.1*max(sim(:,34),sim(:,35));
male_on = sim(:,40)<0.1*max(sim(:,34),sim(:,35));

sel_par = find(mono&no_bi&male_off&male_on);
sel_par_ma=find(mono);

stoch_sim = dlmread('../../simulations/Fig2_3_4/SUMMARY_stoch_sim_Oct2018.txt');
ma = stoch_sim(:,34:134);
ba = stoch_sim(:,275:375);
%sw = stoch_sim(:,226:326);
stoch_sel_par = find(mean(ma(:,81:101),2)>0.8);
%% explore overlap betwe mono no_bi etc

length(sel_par_ma)
sum(mono&no_bi)
sum(mono&male_on&no_bi)

% all parameters with male_on are also no_bi
%%

graph_size=[2.5 2.5];
lw=1.5;
fs=8;
fst=10;
pos_x=[2:4:20];
pos_y=[2 :4.6:50];

p1=[1900, 100,700,800];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')

%% Paper suppl.: Distribution of mono-allelic up-regulation in stochastic simulation
figure(1)
clf
temp=histogram(100*mean(ma(:,81:101),2),10);

a=temp.Values;
b=temp.BinEdges(2:end)-0.5*temp.BinWidth;
yl=[0 3];
a_plot=100*a/sum(a);
a_plot(1)=a_plot(1)*0.01*yl(2);

bar(b,a_plot,1,'FaceColor',[0.5 0.5 0.5])

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(3) pos_y(1) graph_size],...
    'xlim',[0 100],'ylim',yl,'XTick',[0:20:100],'YTick',[0:1:3],'YTickLabel',[0 1 3 100]);

xlabel('XaXi [% cells]','Fontsize',fst)
ylabel('Parameter Sets [%]','Fontsize',fst)

print('../../plots/FigM/FigM2','-depsc')

%% distribution of max BA expression
clf
max_ba = max(ba(stoch_sel_par,:),[],2);
temp=histogram(100*max_ba,10);

a=temp.Values;
b=temp.BinEdges(2:end)-0.5*temp.BinWidth;
yl=[0 100];
a_plot=100*a/sum(a);

bar(b,a_plot,1,'FaceColor',[0.5 0.5 0.5])

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(3) pos_y(1) graph_size],...
    'xlim',[0 100],'ylim',yl,'XTick',[0:20:100],'YTick',[0:25:100]);

xlabel('Bi-allelic [% max]','Fontsize',fst)
ylabel('Parameter Sets [%]','Fontsize',fst)

print('../../plots/FigM/FigM6','-depsc')
%% plot parameter distributions
cm=[14 114 186;217 84 26;237 178 32]/255;

xl_dist=[1 5;1 5;1 5;1 5;0.01 10;0.01 10;0.01 10;0.01 10;0 48;0 48;50 500;50 500;50 500];
var_par = [3 5 11 13 4 6 12 14 7 8 21 22 23];

var_log = [4 6 12 14 21:23];
sim_log = sim;
sim_log(:,var_log) = log10(sim(:,var_log));
stoch_sim_log = stoch_sim;
stoch_sim_log(:,var_log) = log10(stoch_sim(:,var_log));

figure(1)
clf
p1=[1900, 100,700,800];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')
for n=1:length(var_par)

     axes
     hold on;
     if n<5
        [a,b] = hist(sim_log(:,var_par(n)),20);    
        [a2,b] =hist(sim_log(sel_par,var_par(n)),b);
         plot((b),100*a./sum(a),'-','LineWidth',lw,'Color',cm(1,:))
         plot(b,100*a2./sum(a2),'-','LineWidth',lw,'Color',cm(2,:))
         set(gca,'Units','Centimeters','Position',[pos_x(n) pos_y(3) graph_size],'XTick',[1 3 5])

     elseif n<9
        [a,b] = hist(sim_log(:,var_par(n)),20);    
        [a2,b] =hist(sim_log(sel_par,var_par(n)),b);
         plot((b),100*a./sum(a),'-','LineWidth',lw,'Color',cm(1,:))
         plot(b,100*a2./sum(a2),'-','LineWidth',lw,'Color',cm(2,:))
         set(gca,'Units','Centimeters','Position',[pos_x(n-4) pos_y(2) graph_size],'XTick',[-2 -1 0 1],'XTickLabel',[0.01 0.1 1 10])
     else
         [a2,b] =hist(stoch_sim_log(:,var_par(n)),20);
          plot((b),100*a./sum(a),'-','LineWidth',lw,'Color',cm(1,:))
          set(gca,'Units','Centimeters','Position',[pos_x(n-8) pos_y(1) graph_size])
          if n==9|n==10
              set(gca,'XTick',[1 5 10 15])
          else
              set(gca,'XTick',[1.7 2.3 2.699],'XTickLabel',round(10.^[1.7 2.3 2.699],0))
          end
     end
    [a3,b] =hist(stoch_sim_log(stoch_sel_par,var_par(n)),b);
    plot(b,100*a3./sum(a3),'-','LineWidth',lw,'Color',cm(3,:))
    
  %  xlim(log10(xl_dist(n,:)))
    ylim([0 120*max([a2/sum(a2) a3/sum(a3)])])
     
    %title(['p',num2str(var_par(n))])
    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,...
    'xlim',[b(1)-abs(b(2)-b(1)) b(end)+abs(b(2)-b(1))]);
    if n<5|n==9|n==10
        xlabel(['p',num2str(var_par(n))],'Fontsize',fst)
    else
        xlabel(['p',num2str(var_par(n)),' [log]'],'Fontsize',fst)
    end
    if n==4
        [le le2]=legend({'all','XaXi stable (ODE)','XaXa->XaXi (Gillespie)'},'location','Northeast','Box','off');
        le.Position=[le.Position(1)*1.35 le.Position(2) le.Position(3)*1 le.Position(4)*1];
        for q=1:2
            set(le2(q),'Fontsize',fs);
        end
    end
    if n==1|n==5|n==9
        ylabel('Parameter sets [%]','Fontsize',fst)
    end
end

print('../../plots/FigM/FigM3','-depsc')

%% supplement: plot phase diagrams

figure(2)
clf
pert=[0.5 1 2];
p1=[1000, 100,700,400];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')

graph_size=[2.5 2.5];

%pos_x=[2:3.5:20];
pos_x=[9 5.5 2];
pos_y=[2:5:30];
xl=[0 1];
yl=[0 0.7];
xt=[];
reds=[232 152 117;217 83 25;140 44 2]/255;
blues=[102 170 215;0 114 189;11 72 112]/255;
greys=[0.7 0 0;0 0 0; 0 0.5 0.7];  
n=22;% 15
p=stoch_sim(stoch_sel_par(n),1:32);
ode_ss=sim(ismember(sim(:,[3:6 11:14]),p([3:6 11:14]),'rows'),:);
p6=p(6);
p14_all=p(14)*pert;
xist_set=0:0.001:1;
cXR_set=0:0.001:1;
dot_col='k';
dot_size=20;
lin_col=greys(2,:);

xa_tot = 0.5*((1-ode_ss(34)^p(3)/(ode_ss(34)^p(3)+p(4)^p(3)))+(1-ode_ss(35)^p(3)/(ode_ss(35)^p(3)+p(4)^p(3))));
f_xa_tot = (xa_tot^p(11)/(xa_tot^p(11)+p(12)^p(11)));

xa_tot_bi =1-ode_ss(35)^p(3)/(ode_ss(35)^p(3)+p(4)^p(3));
f_xa_tot_bi = (xa_tot_bi^p(11)/(xa_tot_bi^p(11)+p(12)^p(11)));

for k=1:3

axes
p14=p14_all(2);
switch k
    case 1
        f_xa=(1^p(11)/(1^p(11)+p(12)^p(11)));
    case 2
        xa_tot = 0.5*((1-ode_ss(34)^p(3)/(ode_ss(34)^p(3)+p(4)^p(3)))+(1-ode_ss(35)^p(3)/(ode_ss(35)^p(3)+p(4)^p(3))));
        f_xa = (xa_tot^p(11)/(xa_tot^p(11)+p(12)^p(11)));
    case 3
        xa_tot_bi =1-ode_ss(35)^p(3)/(ode_ss(35)^p(3)+p(4)^p(3));
        f_xa = (xa_tot_bi^p(11)/(xa_tot_bi^p(11)+p(12)^p(11)));
end
cXR_from_xist=1-xist_set.^p(5)./(xist_set.^p(5)+p6'.^p(5));
xist_from_cXR=(1-cXR_set.^p(13)./(cXR_set.^p(13)+p14'.^p(13))).*f_xa;
plot(cXR_set,xist_from_cXR,'LineWidth',lw,'Color',[141 198 63]./255)
hold on
plot(cXR_from_xist,xist_set,'LineWidth',lw,'Color',[146 39 143]./255)

xist_from_cXR_from_xist=(1-cXR_from_xist.^p(13)./(cXR_from_xist.^p(13)+p14'.^p(13))).*f_xa;
diff=xist_set-xist_from_cXR_from_xist>0;
[ind]=find(diff(:,2:end)-diff(:,1:(end-1)));
plot_ss_x=nan(size(cXR_from_xist));
plot_ss_x(ind)=cXR_from_xist(ind);
plot_ss_y=nan(size(xist_from_cXR_from_xist));
plot_ss_y(ind)=xist_from_cXR_from_xist(ind);
plot(plot_ss_x',plot_ss_y','.','Markersize',dot_size,'Color',dot_col)
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(k) pos_y(1) graph_size],...
    'ylim',yl,'XTick',xt,'YTick',[]);
%xlabel('cXR','Fontsize',fst)
%ylabel('Xist','Fontsize',fst)
end
print('../../plots/FigM/FigM4','-depsc')


% plot phase diagrams for perturbed parameters and parameter distributions
figure(3)
clf
p1=[1900, 100,700,800];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')

ti={'XaXa','XaXi','XiXi'};
% effect of p14
n=4;

for k=1:3

axes
p14=p14_all;
p6=p(6)*[ 1 1 1];
switch k
    case 1
        f_xa=(1^p(11)/(1^p(11)+p(12)^p(11)));
    case 2
        xa_tot = 0.5*((1-ode_ss(34)^p(3)/(ode_ss(34)^p(3)+p(4)^p(3)))+(1-ode_ss(35)^p(3)/(ode_ss(35)^p(3)+p(4)^p(3))));
        f_xa = (xa_tot^p(11)/(xa_tot^p(11)+p(12)^p(11)));
    case 3
        xa_tot_bi =1-ode_ss(35)^p(3)/(ode_ss(35)^p(3)+p(4)^p(3));
        f_xa = (xa_tot_bi^p(11)/(xa_tot_bi^p(11)+p(12)^p(11)));
end
cXR_from_xist=1-xist_set.^p(5)./(xist_set.^p(5)+p6'.^p(5));
xist_from_cXR=(1-cXR_set.^p(13)./(cXR_set.^p(13)+p14'.^p(13))).*f_xa;

plot(cXR_from_xist,xist_set,'Color',lin_col,'LineWidth',lw)
hold on
set(gca,'colororder',greys,'nextplot','add');
plot(cXR_set,xist_from_cXR,'LineWidth',lw)
xist_from_cXR_from_xist=(1-cXR_from_xist.^p(13)./(cXR_from_xist.^p(13)+p14'.^p(13))).*f_xa;
% indicate steady states
diff=xist_set-xist_from_cXR_from_xist>0;
[ind]=find(diff(:,2:end)-diff(:,1:(end-1)));
plot_ss_x=nan(size(cXR_from_xist));
plot_ss_x(ind)=cXR_from_xist(ind);
plot_ss_y=nan(size(xist_from_cXR_from_xist));
plot_ss_y(ind)=xist_from_cXR_from_xist(ind);
plot(plot_ss_x(1,:)',plot_ss_y(1,:)','.','Markersize',dot_size,'Color',greys(1,:))
plot(plot_ss_x(3,:)',plot_ss_y(3,:)','.','Markersize',dot_size,'Color',greys(1,:))
plot(plot_ss_x(2,:)',plot_ss_y(2,:)','.','Markersize',dot_size,'Color',dot_col)

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(k) pos_y(1) graph_size],...
    'ylim',yl,'XTick',xt,'YTick',[]);

title(ti{k},'Fontsize',fst)

end

% effect of p12
n=3;

for k=1:3

axes
p14=p14_all(2);
p6=p(6)*[1 1 1];
p4=p(4);
p12=p(12)*pert([3 2 1]);
switch k
    case 1
        f_xa=(1^p(11)./(1^p(11)+p12.^p(11)));
    case 2
        xa_tot = 0.5*((1-ode_ss(34)^p(3)./(ode_ss(34)^p(3)+p4.^p(3)))+(1-ode_ss(35)^p(3)./(ode_ss(35)^p(3)+p4.^p(3))));
        f_xa = (xa_tot.^p(11)./(xa_tot.^p(11)+p12.^p(11)));
    case 3
        xa_tot_bi =1-ode_ss(35)^p(3)./(ode_ss(35)^p(3)+p4.^p(3));
        f_xa = (xa_tot_bi.^p(11)./(xa_tot_bi.^p(11)+p12.^p(11)));
end
cXR_from_xist=1-xist_set.^p(5)./(xist_set.^p(5)+p6'.^p(5));
xist_from_cXR=(1-cXR_set.^p(13)./(cXR_set.^p(13)+p14'.^p(13))).*f_xa';
plot(cXR_from_xist,xist_set,'color',lin_col,'LineWidth',lw)
hold on
set(gca,'colororder',greys(1:3,:),'nextplot','add');
plot(cXR_set,xist_from_cXR,'LineWidth',lw)
xist_from_cXR_from_xist=(1-cXR_from_xist.^p(13)./(cXR_from_xist.^p(13)+p14'.^p(13))).*f_xa';
%indicate steady states
diff=xist_set-xist_from_cXR_from_xist>0;
[ind]=find(diff(:,2:end)-diff(:,1:(end-1)));
plot_ss_x=nan(size(cXR_from_xist));
plot_ss_x(ind)=cXR_from_xist(ind);
plot_ss_y=nan(size(xist_from_cXR_from_xist));
plot_ss_y(ind)=xist_from_cXR_from_xist(ind);
plot(plot_ss_x(1,:)',plot_ss_y(1,:)','.','Markersize',dot_size,'Color',greys(1,:))
plot(plot_ss_x(2,:)',plot_ss_y(2,:)','.','Markersize',dot_size,'Color',dot_col)
plot(plot_ss_x(3,:)',plot_ss_y(3,:)','.','Markersize',dot_size,'Color',blues(3,:))

set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x(k) pos_y(2) graph_size],...
    'ylim',yl,'XTick',xt,'YTick',[]);
title(ti{k},'Fontsize',fst)
end

print('../../plots/FigM/FigM5','-depsc')

