addpath('../');

%% Generate random parameter sets
nr_par=400000;

% generate parameter sets using Latin Hypercube sampling (lhsu), the 12
% hill coefficients n are varied between 1 and 5 and the threshold
% parameters K are varied between 0.01 and 10 in a logarithmic fashion
par_n = lhsu(1*ones(12,1),5*ones(12,1),nr_par);
par_k = 10.^lhsu(-2*ones(12,1),1*ones(12,1),nr_par);

p_all = zeros(nr_par,33);
p_all(:,1:2:23) = par_n;
p_all(:,2:2:24) = par_k;

out=zeros(nr_par,41);
tic
% Loop through all parameter sets
for i=1:nr_par
    i
p=p_all(i,:);
% the first 100 000 have one regulator and the rest has two randomly
% selected regulators regulators (parameters 25-32)
if i<100000
    p(24+randperm(8,1))=1;
else
    p(24+randperm(8,2))=1;
end

%% mono-allelic
tspan=[0:100];
options = odeset('NonNegative',[1:18]);
options2 = optimoptions('fsolve','Display','off','FunctionTolerance',1e-18);
lb=zeros(1,18);
up=ones(1,18);

% Simulate female cells from mono-allelic initial conditions
y0=[0.01, 1, zeros(1,16)];
y0(3:2:17) = 1;
y0(4:2:10) = 0.01;
y0(12:2:18) = 1;
p(33)=1;
parfun=@(x)model_f0(x,p);
% integrate equation system for 100h
[T,y] = ode23tb(@model,tspan,y0,options,p);
% solve equations system for dy/dt=0
[y_ss,fval,exitflag]=fsolve(parfun,y(end,:),options2);
if exitflag<1
    y_ss=nan(size(y_ss));
end

% Simulate male cells from an Xist low state
p(33)=0;
parfun_male=@(x)model_f0(x,p);
y0(2:2:18)=0;
[T,y2] = ode23tb(@model,tspan,y0,options,p);
[y2_ss,fval,exitflag]=fsolve(parfun_male,y2(end,:),options2);
if exitflag<1
    y2_ss=nan(size(y2_ss));
end

% Simulate female cells from bi-allelic initial conditions
y0=[0.99, 1, zeros(1,16)];
y0(3:10) = 0.01;
y0(11:18) = 1;
p(33)=1;
[T,y3] = ode23tb(@model,tspan,y0,options,p);
[y3_ss,fval,exitflag]=fsolve(parfun,y3(end,:),options2);
if exitflag<1
    y3_ss=nan(size(y3_ss));
end

% Simulate male cells from an Xist high state
p(33)=0;
y0(2:2:18)=0;
[T,y4] = ode23tb(@model,tspan,y0,options,p);
[y4_ss,fval,exitflag]=fsolve(parfun_male,y4(end,:),options2);
if exitflag<1
    y4_ss=nan(size(y4_ss));
end

out=[p round(y_ss(1:2),6) round(y2_ss(1:2),6) round(y3_ss(1:2),6) round(y4_ss(1:2),6)];

% write parameter sets and simulation results in file 'model_scan.txt'
dlmwrite('../../simulations/Fig1/model_scan.txt',out,'-append')
end
toc
