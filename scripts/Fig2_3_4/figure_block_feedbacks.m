clear 
addpath('../Fig1');
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
stoch_sel_par = find(mean(ma(:,81:101),2)>0.8);

%%
graph_size=[2.5 2.5];
lw=1.5;
fs=8;
fst=10;
pos_x=[2:4:20];
pos_y=[2 :4.6:50];

p1=[1900, 100,700,800];
%set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')

%% plot ODE bifurcation diagram
n = 22; %145
p=stoch_sim(stoch_sel_par(n),1:23);
ode_ss=sim(ismember(sim(:,[4 6 12 14]),p([4 6 12 14]),'rows'),:);

figure(1)
clf
pos_y=[2 :4:50];

xa_tot = ((1-ode_ss(34)^p(3)/(ode_ss(34)^p(3)+p(4)^p(3)))+(1-ode_ss(35)^p(3)/(ode_ss(35)^p(3)+p(4)^p(3))));
xa_tot_bi =2*(1-ode_ss(34)^p(3)/(ode_ss(34)^p(3)+p(4)^p(3)));
xa_tot_no =2*(1-ode_ss(35)^p(3)/(ode_ss(35)^p(3)+p(4)^p(3)));

set(gcf,'OuterPosition',p1,'PaperPositionMode','auto','Units','pixels')
options = odeset('NonNegative',[1:18]);
 
p_ode=p;
p_ode(:,[1,2,7:10,15:33])=0;
p_ode(26:27)=1;
p_ode(33)=1;
p_noFB=p_ode;
p_noFB(14)=1000;
tspan=0:100;
 

%Simulate the steady state for the full model and the models with blocked feedbacks 
%on one allele with different tXA doses and different initial Xist and cXR levels for each tXA dose
y0=zeros(1,18); 
txa=[0:0.1:3];
ini=[0:0.1:1];
clear yend yend_noFB yend_notxa
 for z=1:length(txa)
     for q=1:length(ini)
        p_ode(15)=txa(z);
        p_noFB(15)=txa(z);
        y0([1 7])=[ini(q) 1-ini(q)];
        [T,y] = ode15s(@model_chr,tspan,y0,options,p_ode);
        yend(z,q)=y(end,1);
        %[z q  y(end,1)]
        [T,y] = ode15s(@model_chr,tspan,y0,options,p_noFB);
        yend_noFB(z,q)=y(end,1);
        %[z q  y(end,1)]
        p_ode(15)=xa_tot;
        [T,y] = ode15s(@model_chr,tspan,y0,options,p_ode);
        %[z q  y(end,1)]
        yend_notxa(z,q)=y(end,1);
     end
 end
 
clf
 pos_x3=[2:3:20];
 axes   
 plot(repmat(txa',1,size(yend,2)),yend,'k.','LineWidth',lw)

 set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x3(1) pos_y(2) graph_size],...
        'xlim',[0 2.2],'ylim',[-0.05 0.7],'xticklabel',[0:2],'ytick',[0:0.2:1]);%,'XTick',[xa_tot_no xa_tot xa_tot_bi]);
    
xlabel('tXA dose','Fontsize',fst)
ylabel('Xist','Fontsize',fst)

 axes   
  plot(repmat(txa',1,size(yend_noFB,2)),yend_noFB,'k.','LineWidth',lw)

 set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x3(2) pos_y(2) graph_size],...
        'xlim',[0 2.2],'xticklabel',[0:2],'ylim',[-0.05 0.7],'xticklabel',[0:2],'ytick',[0:0.2:1],'yticklabel',[]);%,'XTick',[xa_tot_no xa_tot xa_tot_bi]);
xlabel('tXA dose','Fontsize',fst)

    
axes   
 plot(repmat(txa',1,size(yend_notxa,2)),yend_notxa,'k.','LineWidth',lw)

 set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x3(3) pos_y(2) graph_size],...
        'xlim',[0 2.2],'ylim',[-0.05 0.7],'xticklabel',[0:2],'ytick',[0:0.2:1],'yticklabel',[]);%,'XTick',[xa_tot_no xa_tot xa_tot_bi]);
xlabel('tXA dose','Fontsize',fst)



% Simulate the global steady states for the full model and the models with
% blocked feedbacks starting from different Xist, cXR and tXA initial
% levels
p_ode(15)=0;
p_notXA=p_ode;
p_notXA(15)=xa_tot;
p_notXA(26)=2;
clear y_end y_end_noFB y_end_notxa y_end_sym y_end_noFB_sym y_end_notxa_sym

%ini=[0:0.1:1];
ini=[0,0.01,0.1:0.1:1];
i=1;
i2=1;
for q1=1:length(ini)
    for q2=1:length(ini)
        y0=zeros(1,18); 
        y0(1:2)=[ini(q1) ini(q2)];
        y0(5:8)=1-[y0(1) y0(2) y0(1) y0(2)];
        y0(5) = [1-y0(1)^p_ode(3)/(y0(1)^p_ode(3)+p_ode(4)^p_ode(3))];
        y0(6) = [1-y0(2)^p_ode(3)/(y0(2)^p_ode(3)+p_ode(4)^p_ode(3))];
        y0(7) = [1-y0(1)^p_ode(5)/(y0(1)^p_ode(5)+p_ode(6)^p_ode(5))];
        y0(8) = [1-y0(2)^p_ode(5)/(y0(2)^p_ode(5)+p_ode(6)^p_ode(5))];
        [T,y1] = ode15s(@model,tspan,y0,options,p_ode);
        [T,y2] = ode15s(@model,tspan,y0,options,p_noFB);
        [T,y3] = ode15s(@model_mod,tspan,y0,options,p_notXA);
        if q1==q2
             y_end_sym(i,:)=y1(end,:);
             y_end_noFB_sym(i,:)=y2(end,:);
             y_end_notxa_sym(i,:)=y3(end,:);
             i=i+1;
        else
            y_end(i2,:)=y1(end,:);
             y_end_noFB(i2,:)=y2(end,:);
             y_end_notxa(i2,:)=y3(end,:);
             i2=i2+1;
        end
    end
end

axes
plot(y_end(:,1), y_end(:,2),'o','Markersize',5,'MarkerFaceColor','k','Color','k')
hold on
plot(y_end_sym(:,1), y_end_sym(:,2),'o','Markersize',5,'Color','k')
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x3(1) pos_y(1) graph_size],...
        'xlim',[-0.05 0.4],'XTick',[0:0.2:4],'ylim',[-0.05 0.4],'ytick',[0:0.2:4]);
xlabel('Xist 1','Fontsize',fst)
ylabel('Xist 2','Fontsize',fst)
    
axes
plot(y_end_noFB(:,1), y_end_noFB(:,2),'o','Markersize',5,'MarkerFaceColor','k','Color','k')
hold on
plot(y_end_noFB_sym(:,1), y_end_noFB_sym(:,2),'o','Markersize',5,'Color','k')
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x3(2) pos_y(1) graph_size],...
        'xlim',[-0.05 0.4],'XTick',[0:0.2:4],'ylim',[-0.05 0.4],'ytick',[0:0.2:4],'yticklabel',[]);
xlabel('Xist 1','Fontsize',fst)

axes
plot(y_end_notxa(:,1), y_end_notxa(:,2),'o','Markersize',5,'MarkerFaceColor','k','Color','k')
hold on
plot(y_end_notxa_sym(:,1), y_end_notxa_sym(:,2),'o','Markersize',5,'Color','k')
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x3(3) pos_y(1) graph_size],...
        'xlim',[-0.05 0.4],'XTick',[0:0.2:4],'ylim',[-0.05 0.4],'ytick',[0:0.2:4],'yticklabel',[]);
xlabel('Xist 1','Fontsize',fst)

print('../../plots/Fig2/Fig2E_F_G','-depsc')
