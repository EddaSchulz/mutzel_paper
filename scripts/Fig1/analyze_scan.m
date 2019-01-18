% read in simulation results for model comparison
sim = dlmread('../../simulations/Fig1/model_scan.txt');
%% classify each parameter set as mono-allelic, not bi-allelic, 
% expression in males starting from and off state and expression in males starting from and Xist on state
mono=(sim(:,34)./sim(:,35))<0.1|(sim(:,34)./sim(:,35))>10;
no_bi=(sim(:,38)<0.1*max(sim(:,34),sim(:,35)))|(sim(:,39)<0.1*max(sim(:,34),sim(:,35)));
male_off = sim(:,36)<0.1*max(sim(:,34),sim(:,35));
male_on = sim(:,40)<0.1*max(sim(:,34),sim(:,35));
%% analyze maintenance of XaXi state in 1reg models
mod_sum=nan(8,3);
for m=1:8
    ind=zeros(1,8);
    ind(m)=1;
    mod_par=find(ismember(sim(:,25:32),ind,'rows'));
    mod_sum(m,1)=length(mod_par);
    mod_sum(m,2)=sum(mono(mod_par));
end
mod_sum(:,3)=100*mod_sum(:,2)./mod_sum(:,1);
% print results from model comparison for 1reg models: 
% col1: #par sets,
% col2: #MA par sets, 
% col3: % MA par sets (shown in Fig. 1B)
mod_sum
%% analyze maintenance of XaXi state in  2reg models
tot=nan(8,8);
ma=nan(8,8);
for m=1:8
    for n=1:8
    ind=zeros(1,8);
    ind([m n])=1;
    mod_par=find(ismember(sim(:,25:32),ind,'rows'));
    tot(m,n)=length(mod_par);
    ma(m,n)=sum(mono(mod_par));
    end
end
ma_frac=100*ma./tot;
order=[1 5 2 6 3 7 4 8];
mods={'cXA','tXA','cXR','tXR','ecXA','etXA','ecXR','etXR'};
%print results for 2reg models: fraction of MA par sets (Table M1)
mods(order)
ma_frac(order,order)

%% analyze no_bi and male simulations for 2reg models with cXR
cXRmod_sum=nan(8,5);
for m=1:8
    ind=zeros(1,8);
    ind(3)=1;
    ind(m)=1;
    mod_par=(ismember(sim(:,25:32),ind,'rows'));
    cXRmod_sum(m,1)=sum(mod_par);
    cXRmod_sum(m,2)=length(find(mod_par&mono));
    cXRmod_sum(m,3)=sum(no_bi(find(mod_par&mono)));
    cXRmod_sum(m,4)=sum(male_off(find(mod_par&mono)));
    cXRmod_sum(m,5)=sum(male_on(find(mod_par&mono)));
    cXRmod_sum(m,6)=sum(male_on(find(mod_par&mono))&no_bi(find(mod_par&mono)));
end
% print results from 2 reg models with cXR (col1:3 shown in Table M2)
% col1: % XiXi unstable within mono-allelic parameter sets
% col2: % Xa stable within mono-allelic parameter sets
% col3: % Xi unstable within mono-allelic parameter sets
% col4: % XiXi and Xi unstable within mono-allelic parameter sets

cXRmod_frac=100*cXRmod_sum(:,3:6)./repmat(cXRmod_sum(:,2),1,4);
mods(order)
cXRmod_frac(order,:)

%% plot trajectories of ma expression for cxR models: figure 1 main text
tspan=[0:100];
options = odeset('NonNegative',[1:14]);

graph_size=[1.5 1.5];
lw=1.5;
fs=8;
fst=10;
col_x1=[177 206 85]/255;
col_x2=[35 155 56]/255;

clf
q=ones(8,8);

plot_order=[1 3 5 7 2 4 6 8];
pos_x2 = ([2 5:2.3:30]);
pos_y=fliplr([2:1.7:18]);
pos_y = pos_y(plot_order);

y0=[0.01, 1, zeros(1,16)];
y0(3:2:17) = 1;
y0(4:2:10) = 0.01;
y0(12:2:18) = 1;

for m1=1:8
    m2=m1;
        ind=zeros(1,8);
        ind([m1 m2])=1;
        sel_par = ismember(sim(:,25:32),ind,'rows');
        temp = find(sel_par&mono);
        if isempty(temp)
             temp = find(sel_par); 
        end          
        p=sim(temp(q(m1,m2)),:);      
                   
        % plot example trajectory
        
        p(33)=1;
        [T,y] = ode15s(@model,tspan,y0,options,p);
            
        axes
        plot(T,y(:,1),'-','LineWidth',lw,'Color',col_x1);
        hold on;
        plot(T,y(:,2),'-','LineWidth',lw,'Color',col_x2);
        set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x2(1) pos_y(m1) graph_size],...
        'xlim',[-3 96],'XTick',[],'XTickLabel',[],'ylim',[-0.3 1.5],'ytick',[]);
        if m1==8
               set(gca,'XTick',[0 48 96])
        end
        set(gca,'XTickLabel',[0 2 4])
        if ma_frac(m1,m2)>0
            col_text='r';
        else
            col_text='k';
        end
        
        text(48,1.2,[num2str(round(ma_frac(m1,m2),2,'significant')),'%'],'Fontsize',8,'Color',col_text)
end

%plot example tranjectories for bi-allelic expression

for m1=1:8
    m2=3;
        ind=zeros(1,8);
        ind([m1 m2])=1;
        sel_par = ismember(sim(:,25:32),ind,'rows');
        temp = find(sel_par&mono&no_bi);
        if isempty(temp)
             temp = find(sel_par&mono); 
        end          
        p=sim(temp(q(m1,m2)),:);      
                   
        % plot example trajectory
        p(33)=1;
        y0=[0.99, 1, zeros(1,16)];
        y0(3:10) = 0.01;
        y0(11:18) = 1;

        [T,y] = ode15s(@model,tspan,y0,options,p);
            
        axes
        plot(T,y(:,1),'-','LineWidth',lw,'Color',col_x1);
        hold on;
        plot(T,y(:,2),'-','LineWidth',lw,'Color',col_x2);
        line([0 100],[p(34) p(34)],'Color',[0.5 0.5 0.5],'Linestyle',':');
        line([0 100],[p(35) p(35)],'Color',[0.5 0.5 0.5],'Linestyle',':');
        set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x2(4) pos_y(m1) graph_size],...
        'xlim',[-3 96],'XTick',[],'XTickLabel',[],'ylim',[-0.3 1.5],'ytick',[]);
        if m1==8
               set(gca,'XTick',[0 48 96],'XTickLabel',[0 2 4])
        end
        if cXRmod_sum(m1,3)>0
            col_text='r';
        else
            col_text='k';
        end
        text(48,1.2,[num2str(round(100*cXRmod_sum(m1,3)./cXRmod_sum(m1,2),2,'significant')),'%'],'Fontsize',8,'Color',col_text)
        
        %plot male on
        temp = find(sel_par&mono&male_on);
        if isempty(temp)
             temp = find(sel_par&mono); 
        end          
        p=sim(temp(q(m1,m2)),:);      
        p(33)=0;
        y0(2:2:18)=0;
        [T,y] = ode15s(@model,tspan,y0,options,p);
            
        axes
        plot(T,y(:,1),'-','LineWidth',lw,'Color',col_x2);
        hold on;
        line([0 100],[p(34) p(34)],'Color',[0.5 0.5 0.5],'Linestyle',':');
        line([0 100],[p(35) p(35)],'Color',[0.5 0.5 0.5],'Linestyle',':');
        set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x2(5) pos_y(m1) graph_size],...
        'xlim',[-3 96],'XTick',[],'XTickLabel',[],'ylim',[-0.3 1.5],'ytick',[]);
        if m1==8
            set(gca,'XTickLabel',[0 2 4],'XTick',[0 48 96])
        end
        if cXRmod_sum(m1,5)>0
            col_text='r';
        else
            col_text='k';
        end
        text(48,1.2,[num2str(round(100*cXRmod_sum(m1,5)./cXRmod_sum(m1,2),2,'significant')),'%'],'Fontsize',8,'Color',col_text)
        
        %plot male_off
        temp = find(sel_par&mono&male_off);
        if isempty(temp)
             temp = find(sel_par&mono); 
        end          
        p=sim(temp(q(m1,m2)),:);      
                   
        % plot example trajectory

        y0=[0.01, 0, zeros(1,16)];
        y0(3:2:18) = 1;
        [T,y] = ode15s(@model,tspan,y0,options,p);
            
        axes
        plot(T,y(:,1),'-','LineWidth',lw,'Color',col_x2);
        hold on;
        line([0 100],[p(34) p(34)],'Color',[0.5 0.5 0.5],'Linestyle',':');
        line([0 100],[p(35) p(35)],'Color',[0.5 0.5 0.5],'Linestyle',':');
        set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',8,'Units','Centimeters','Position',[pos_x2(3) pos_y(m1) graph_size],...
        'xlim',[-3 96],'XTick',[],'XTickLabel',[],'ylim',[-0.3 1.5],'ytick',[]);
        if m1==8
            set(gca,'XTickLabel',[0 2 4],'XTick',[0 48 96])
        end
        if cXRmod_sum(m1,4)>0
            col_text='r';
        else
            col_text='k';
        end
        
        text(48,1.2,[num2str(round(100*cXRmod_sum(m1,4)./cXRmod_sum(m1,2),2,'significant')),'%'],'Fontsize',8,'Color',col_text)
end

print('../../plots/Fig1/Fig1B_C','-depsc')



