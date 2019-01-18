clear

filename='../../simulations/Fig5_7/sim_2C_k1XAtdep_model_wo_X_rep_mut_sets_WT_170915.txt';
sim{1}=dlmread(filename);

filename='../../simulations/Fig5_7/sim_2C_k1XAtdep_model_wo_X_rep_mut_sets_hetTsix_170915.txt';
sim{2}=dlmread(filename);

filename='../../simulations/Fig5_7/sim_2C_k1XAtdep_model_wo_X_rep_mut_sets_homTsix_170915.txt';
sim{3}=dlmread(filename);

filename='../../simulations/Fig5_7/sim_2C_k1XAtdep_model_wo_X_rep_mut_sets_hetXist_170915.txt';
sim{4}=dlmread(filename);

%%

t_half=nan(100,4);
for a=1:4
    [C,ia,ic]=unique(sim{a}(:,1:15),'rows');
    ind1=16:116;
    ind2=117:217;
    for n=1:length(C)
        ma_frac{a}(n,:)=mean(((sim{a}(ic==n,ind1)<10&sim{a}(ic==n,ind2)>10))|((sim{a}(ic==n,ind1)>10&sim{a}(ic==n,ind2)<10)));
        ba_frac{a}(n,:)=mean((sim{a}(ic==n,ind1)>10&sim{a}(ic==n,ind2)>10));
        chr1_frac{a}(n,:)=mean((sim{a}(ic==n,ind1)>10&sim{a}(ic==n,ind2)<10));
        chr2_frac{a}(n,:)=mean((sim{a}(ic==n,ind1)<10&sim{a}(ic==n,ind2)>10));
        xist_on_frac{a}(n,:)=mean((sim{a}(ic==n,ind1)>10|sim{a}(ic==n,ind2)>10));
        sum_dat{a}(n,1:4)=[ma_frac{a}(n,end) ba_frac{a}(n,end) chr1_frac{a}(n,end) chr2_frac{a}(n,end)];
        temp=find(sum((sim{a}(ic==n,ind1)>10)|(sim{a}(ic==n,ind2)>10))>50);
        t_half(n,a)=temp(1)-1;
        chr1_mean{a}(n,:)=mean(sim{a}(ic==n,ind1));
        chr2_mean{a}(n,:)=mean(sim{a}(ic==n,ind2));
    end
end
%% plotting parameters
graph_size=[2.5 2.5];
lw=2;
fs=8;
fst=10;
pos_x=[2.7:4.2:15.5];
pos_y=[2 6.5 11 15.5 19.5:3.8:30];

p1=[1800, 100,700,1000];
set(gcf,'OuterPosition',p1,'PaperPositionMode','auto')

%%
clf
ti={'wt','hetTsix','homTsix','hetXist'};
lab{1}={'Xi: Wt1','Xi: Wt2','Xi: Wt1&2'};
lab{2}={'Xi: Wt','Xi: Mut','Xi: Wt&Mut'};
lab{3}={'Xi: Mut1','Xi: Mut2','Xi: Mut1&2'};
lab{4}={'Xi: Wt','Xi: Mut','Xi: Wt&Mut'};

col{1}=[0 0 0;0 0 0;0 0 0];
col{2}=[0 0 0;1 0 0;0.5 0 0];
col{3}=[1 0 0;1 0 0;1 0 0];
col{4}=[0 0 0;1 0 0;0.5 0 0];



% boxplot %cells that express Xist from each chomosome
for n=1:4
    ax=axes
    bp=boxplot(100*sum_dat{n}(:,[4 3 2]),'Color',col{n},'Symbol','.','OutlierSize',4,'Labels',lab{n})
    set(bp(:,:),'LineWidth',lw);
    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(n) pos_y(4) graph_size],...
    'ylim',[-10 110]);

    %title(ti{n})
    ax.XTickLabelRotation=45;
    if n==1
        ylabel('Cells [%]','Fontsize',fst);
    end
end

% scatter plot: example parameter set, plot xist level from each chromosome
q=39;
lab{1}={'# Xist (Wt1)','# Xist (Wt2)'};
lab{2}={'# Xist (Wt)','# Xist (Mut)'};
lab{3}={'# Xist (Mut1)','# Xist (Mut2)'};
lab{4}={'# Xist (Wt)','# Xist (Mut)'};
col_x={'k','k','r','k'};
col_y={'k','r','r','r'};


for n=1:4
    ax=axes
    plot(sim{n}(ic==q,217),sim{n}(ic==q,116),'k.')
    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(n) pos_y(5) graph_size],...
    'xlim',[-50 600],'ylim',[-50 600]);
    xlabel(lab{n}{1},'color',col_x{n});
    ylabel(lab{n}{2},'color',col_y{n});
end

% plot time course for example cell in example parameter set
time=0:100;
tmp=find(ic==q);
p=6;
col{1}=[0 0 0;0 0 0];
col{2}=[1 0 0;0 0 0];
col{3}=[1 0 0;1 0 0];
col{4}=[1 0 0;0 0 0];

for n=1:4
    ax=axes
    set(gca,'defaultAxesColorOrder',col{n});
    if n==3
        plot(time, sim{n}(tmp(p),16:116)+10,'LineWidth',1,'Color',col{n}(1,:))
        hold on;
        plot(time, sim{n}(tmp(p),117:217),'LineWidth',1,'Color',col{n}(2,:))
    else
        plot(time, sim{n}(tmp(p),16:116),'LineWidth',1,'Color',col{n}(1,:))
        hold on;
        plot(time, sim{n}(tmp(p),117:217),'LineWidth',1,'Color',col{n}(2,:))
    end 
    set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(n) pos_y(6) graph_size],...
    'xlim',[0 100],'ylim',[-50 600],'XTick',[0:24:96],'XTickLabel',[0:4]);
    xlabel('Time [Days]','Fontsize',fst);
    ylabel('# Xist','Fontsize',fst);
end

% plot comparative kinetics for one parameter set
%q=3;
axes
line([-10 110],[50 50],'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on
pl=plot(time,100*xist_on_frac{1}(q,:),'k-',time,100*xist_on_frac{2}(q,:),'r:',time,100*xist_on_frac{4}(q,:),'r-','Linewidth',lw)
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(1) pos_y(3) 4 2.5],...
    'xlim',[-2 40],'ylim',[-5 105],'XTick',[0:12:96],'Box','on');
xlabel('Time [h]','Fontsize',fst);
ylabel('Cells [%]','Fontsize',fst);
le=legend(pl,{'WT','Tsix^{+/-}','Xist^{+/-}'},'Fontsize',fst)
le.Position(1)=0.34;
le.Position(3)=0.008;
plot(t_half(q,[1])-0.5,[50],'ko','LineWidth',2)
plot(t_half(q,[2 4])-0.5,[50 50],'ro','LineWidth',2)


axes
rel_t_half=t_half(:,2:end)-repmat(t_half(:,1),1,3);
bins=linspace(-20,20,41);
a = histogram(rel_t_half(:,1),bins);
a2=a.Values;
a2_width=a.BinWidth;
a = histogram(rel_t_half(:,2),bins);
b2=a.Values;
a = histogram(rel_t_half(:,3),bins);
c2=a.Values;
plot(bins(2:end)-a2_width/2,a2,'r:',bins(2:end)-a2_width/2,c2,'r','Linewidth',lw)
set(gca,'TickLength',[0.02 0],'TickDir','out','Linewidth',1,'Fontsize',fs,'Units','Centimeters','Position',[pos_x(3) pos_y(3) graph_size],...
    'xlim',[-20 20],'ylim',[0 25],'XTick',[-20:10:20]);
xlabel('\DeltaT_{1/2}','Fontsize',fst);
ylabel('# parameter sets','Fontsize',fst);

print('../../plots/Fig7/Fig7','-depsc','-loose')
