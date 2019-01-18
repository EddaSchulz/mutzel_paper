
filenames = dir('../../simulations/Fig2_3_4/stoch_sim_files/TEMP_stoch_sim_Oct2018_*.txt');

for i = 0:(size(filenames,1)-1)

    clearvars -except i
    i
    current = sprintf('../../simulations/Fig2_3_4/stoch_sim_files/TEMP_stoch_sim_Oct2018_%d.txt',i);
    if exist(current, 'file')
        sim=dlmread(current,',');
    else
        'something wrong here - File does not exist'
        i
    end
    
    
    nr_p = 33;
    nr_x = 2;
    nr_cells = 100;
    
    
    
    t = (size(sim,2)-nr_p)./nr_x;
    out = calc_ma_ba_swon(sim, nr_p, nr_x);
    frac_ma = out(:,nr_p+1:nr_p+t);
    frac_ba = out(:,nr_p+t+1: nr_p+2*t);
    swon = out(:,end);
    
    sel_par = find(mean(frac_ma(:,81:101),2)>0.8);
        sel_sets = [];
    for j = 1:length(sel_par)
        sel_sets = [sel_sets; sim((sel_par(j)-1)*nr_cells+1:sel_par(j)*nr_cells,:)];
    end
    %Uncomment to construct summary files
    %dlmwrite('../../simulations/Fig2_3_4/SUMMARY_stoch_sim_Oct2018.txt', out, '-append');
    %dlmwrite('../../simulations/Fig2_3_4/MA_sets_stoch_sim_Oct2018.txt', sel_sets, '-append');  
end
