function [out] = calc_ma_ba_swon(sim, nr_par, nr_x)

[C, ia, ib] = unique(sim(:,1:nr_par), 'rows', 'stable');
nr_cells = size(sim,1)./length(ia);
t = (size(sim,2)-nr_par)./nr_x;
sim_bi = zeros(size(sim,1), size(sim,2));
sim_bi(sim>=(sim(:,2).*sim(:,21))./5) = 1;

if nr_x<=2
    nr_xi= [sim_bi(:,nr_par+1:nr_par+t)+sim_bi(:,nr_par+t+1:nr_par+2*t)];
elseif nr_x==3
    nr_xi= [sim_bi(:,nr_par+1:nr_par+t)+sim_bi(:,nr_par+t+1:nr_par+2*t)+sim_bi(:,nr_par+2*t+1:nr_par+3*t)];
elseif nr_x==4
    nr_xi= [sim_bi(:,nr_par+1:nr_par+t)+sim_bi(:,nr_par+t+1:nr_par+2*t)+sim_bi(:,nr_par+2*t+1:nr_par+3*t)+sim_bi(:,nr_par+3*t+1:nr_par+4*t)];
end

frac_ma = zeros(size(nr_xi,1)./nr_cells, size(nr_xi,2));
frac_ba = zeros(size(nr_xi,1)./nr_cells, size(nr_xi,2));
frac_tri = zeros(size(nr_xi,1)./nr_cells, size(nr_xi,2));
frac_tetra = zeros(size(nr_xi,1)./nr_cells, size(nr_xi,2));
rat = NaN(size(nr_xi,1)./nr_cells, size(nr_xi,2));

[r c] = find(sim_bi(:,nr_par+1:nr_par+t));
firstIndex1 = accumarray(r,c,[size(sim_bi(:,nr_par+1:nr_par+t),1),1],@min,t);
[r c] = find(sim_bi(:,nr_par+t+1:nr_par+2*t));
firstIndex2 = accumarray(r,c,[size(sim_bi(:,nr_par+t+1:nr_par+2*t),1),1],@min,t);
temp = min(firstIndex1, firstIndex2);
if nr_x>=3
    [r c] = find(sim_bi(:,nr_par+2*t+1:nr_par+3*t));
    firstIndex3 = accumarray(r,c,[size(sim_bi(:,nr_par+2*t+1:nr_par+3*t),1),1],@min,t);
    temp = min(temp,firstIndex3);
    if nr_x>=4
        [r c] = find(sim_bi(:,nr_par+3*t+1:nr_par+4*t));
        firstIndex4 = accumarray(r,c,[size(sim_bi(:,nr_par+3*t+1:nr_par+4*t),1),1],@min,t);
        temp = min(temp,firstIndex4);
    end
end
sw = temp;
swon = zeros(size(nr_xi,1)./nr_cells,1);

for i = 1:length(ia)
    for j = 1:size(nr_xi,2)
        frac_ma(i,j) = length(find(nr_xi((i-1)*nr_cells+1:nr_cells*i,j)==1))/nr_cells;
        frac_ba(i,j) = length(find(nr_xi((i-1)*nr_cells+1:nr_cells*i,j)==2))/nr_cells;
        frac_tri(i,j) = length(find(nr_xi((i-1)*nr_cells+1:nr_cells*i,j)==3))/nr_cells;
        frac_tetra(i,j) = length(find(nr_xi((i-1)*nr_cells+1:nr_cells*i,j)==4))/nr_cells;
        swon(i) = mean(sw((i-1)*nr_cells+1:nr_cells*i));
    end
end

out = [sim(ia,1:nr_par) frac_ma frac_ba frac_tri frac_tetra swon];
end
