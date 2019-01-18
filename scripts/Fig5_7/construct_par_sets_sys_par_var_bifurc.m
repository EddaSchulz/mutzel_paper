% Script to construct parameter sets for simulation of 1C Xist/Tsix model

% FULL MODEL
limits = [10 1000; 10 1000; 0.0001 0.1];
i = 20;
k1 = logspace(log10(limits(1,1)), log10(limits(1,2)), i);
k2 = logspace(log10(limits(2,1)), log10(limits(2,2)), i);
%k8*k3 => k3=1000 => effective krev_rep = k8*1000
k8 = logspace(log10(limits(3,1)), log10(limits(3,2)), i);

p = [1,1,1,0.1733,1.3868,1,1,1,1,1,1,1,1,1,1];
%p7 = 0: Xist repressed by Tsix pol 
p(7) = 0;
% p9 = 1: pol coll exist
p(9) = 1;
%p13 =1 => Tsix silenced by Xist
p(13) = 1;

p(3) = 1000;
p(10) = 0;

for z = 1:size(k1,2)
    p(1) = k1(z);
    for y = 1: size(k2,2)
        p(2) = k2(y);
        for a = 1:size(k8,2)
            p(8)=k8(a);
            dlmwrite('../../simulations/Fig5_7/parameters_sys_bifurc_full_model_161117.txt', p, '-append');
		end
	end
end

%% NO POL COLL
limits = [10 1000; 10 1000];
i = 20;
k1 = logspace(log10(limits(1,1)), log10(limits(1,2)), i);
k2 = logspace(log10(limits(2,1)), log10(limits(2,2)), i);
k8 = logspace(log10(0.0001), log10(0.1), i); 


p = [1,1,1,0.1733,1.3868,1,1,1,1,1,0,1,1,1,1];
%p7 = 0: Xist repressed by Tsix pol
p(7) = 0;
%p9 = 0: no pol coll
p(9) = 0;
%p13 =1 => Tsix silenced by Xist
p(13) = 1;

%No basal prom trans => p3 very high, p10 low
p(3) = 1000;
p(10)=0;


temp = [];
for z = 1:size(k1,2)
    p(1) = k1(z);
    for y = 1: size(k2,2)
        p(2) = k2(y);
        for x = 1: size(k8,2)
            p(8) = k8(x);
            dlmwrite('../../simulations/Fig5_7/parameters_sys_bifurc_no_pol_coll_161121.txt', p, '-append');
        end
    end
end

%% NO (OR PARTIAL) TSIX SILENCING
limits = [10 1000; 10 1000];
i = 20;
k1 = logspace(log10(limits(1,1)), log10(limits(1,2)), i);
k2 = logspace(log10(limits(2,1)), log10(limits(2,2)), i);
k8 = logspace(log10(0.0001), log10(0.1), i); %k8*k3 => k3=1000 => effective k8 = k8*1000
k11 = 0:0.1:1;

p = [1,1,1,0.1733,1.3868,1,1,1,1,1,1,1,1,1,1];
%p7 = 0: Xist repressed by Tsix pol (=> k8 has a meaning), p7 = 1 => no repression of Xist by Tsix pol = very high k8
p(7) = 0;
% p9 = 1: pol coll exist, p9 = 0: no pol coll
p(9) = 1;
%p13 very high => Tsix not completely silenced by Xist
p(13) = 10000;
%No basal prom trans => p3 very high, p10 low
p(3) = 1000;
p(10)=0;


temp = [];
for z = 1:size(k1,2)
    p(1) = k1(z);
    for y = 1: size(k2,2)
        p(2) = k2(y);
        for x = 1: size(k8,2)
            p(8) = k8(x);
            for a = 1:size(k11,2)
                p(11)=k11(a);
                dlmwrite('../../simulations/Fig5_7/parameters_sys_bifurc_part_Tsix_sil_161118.txt', p, '-append');
            end
        end
    end
end


%% NO XIST PROMOTER REPRESSION
limits = [10 1000; 10 1000];
i = 20;
k1 = logspace(log10(limits(1,1)), log10(limits(1,2)), i);
k2 = logspace(log10(limits(2,1)), log10(limits(2,2)), i);


p = [1,1,1,0.1733,1.3868,1,1,1,1,1,0,1,1,1,1];
% p7 = 1 => no repression of Xist by Tsix pol = very high k8
p(7) = 1;
% p9 = 1: pol coll exist
p(9) = 1;
%p13 =1 => Tsix silenced by Xist
p(13) = 1;
%No basal prom trans => p3 very high, p10 low
p(3) = 1000;
p(10)=0;


temp = [];
for z = 1:size(k1,2)
    p(1) = k1(z);
    for y = 1: size(k2,2)
        p(2) = k2(y);
        for x = 1: size(k8,2)
            p(8) = k8(x);
            dlmwrite('../../simulations/Fig5_7/parameters_sys_bifurc_no_X_rep_170508.txt', p, '-append');
        end
    end
end



% ONLY POLYMERASE COLLISIONS
limits = [10 1000; 10 1000];
i = 20;
k1 = logspace(log10(limits(1,1)), log10(limits(1,2)), i);
k2 = logspace(log10(limits(2,1)), log10(limits(2,2)), i);
k8 = logspace(log10(0.0001), log10(0.1), i); %k8*k3 => k3=1000 => effective k8 = k8*1000


p = [1,1,1,0.1733,1.3868,1,1,1,1,1,0,1,1,1,1];
%p7 = 0: Xist repressed by Tsix pol (=> k8 has a meaning), p7 = 1 => no repression of Xist by Tsix pol = very high k8
p(7) = 0;
% p9 = 1: pol coll exist, p9 = 0: no pol coll
p(9) = 0;
%p13 very high => Tsix not completely silenced by Xist
p(13) = 1;
%No basal prom trans => p3 very high, p10 low
p(3) = 1000;
p(10)=0;


p(9)=1;
p(7)=1;
p(13)=10000;
p(11)=0;
temp = [];
for z = 1:size(k1,2)
    p(1) = k1(z);
    for y = 1: size(k2,2)
        p(2) = k2(y);
        dlmwrite('../../simulations/Fig5_7/parameters_sys_bifurc_only_pol_coll_170424.txt', p, '-append');
    end
end


% ONLY REPRESSION OF XIST BY TSIX
p(9)=0;
p(7)=0;
p(13)=10000;
p(11)=0;
temp = [];
for z = 1:size(k1,2)
    p(1) = k1(z);
    for y = 1: size(k2,2)
        p(2) = k2(y);
        for x = 1: size(k8,2)
            p(8) = k8(x);
            dlmwrite('../../simulations/Fig5_7/parameters_sys_bifurc_only_rep_X_by_T_170424.txt', p, '-append');
        end
    end
end

% ONLY SILENCING OF TSIX BY XIST
p(9)=0;
p(7)=1;
p(13)=1;
p(11)=0;
temp = [];
for z = 1:size(k1,2)
    p(1) = k1(z);
    for y = 1: size(k2,2)
        p(2) = k2(y);
        dlmwrite('../../simulations/Fig5_7/parameters_sys_bifurc_only_sil_T_170424.txt', p, '-append');
    end
end
