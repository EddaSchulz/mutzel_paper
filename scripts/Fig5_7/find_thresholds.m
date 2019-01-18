function out=find_thresholds(mono_mat, bif_p_mat, xa_mat, high_thresh, low_thresh)

nr_cols=size(mono_mat,2);
[r,c] = find(mono_mat>high_thresh);
sel_cols=accumarray(r,c,[],@min);
sel_cols2=sel_cols-1;

sel_cols3=accumarray(r,c,[],@max);
sel_cols3(sel_cols3==nr_cols)=0;
sel_cols4=sel_cols3+1;
sel_cols4(sel_cols3<1)=0;

[r,c] = find(mono_mat>low_thresh);
sel_cols5=zeros(size(mono_mat,1),1);
sel_cols5(1:max(unique(r)))=accumarray(r,c,[],@min);
sel_cols6=sel_cols5-1;
sel_cols6(sel_cols<1)=0;

sel_cols5(sel_cols5>sel_cols2)=0;
sel_cols5(sel_cols==0)=0;

sel_cols7=accumarray(r,c,[],@max);
sel_cols7(sel_cols7==nr_cols)=0;
sel_cols8=sel_cols7+1;
sel_cols8(sel_cols7<1)=0;
sel_cols8(sel_cols3<1)=0;
sel_cols7(sel_cols7<sel_cols4)=0;
sel_cols7(sel_cols3==0)=0;

mo_low=nan(length(sel_cols6),1);
mo_low(sel_cols6>0)=bif_p_mat(sub2ind(size(bif_p_mat),find(sel_cols6>0),sel_cols6(sel_cols6>0)));
bi_low=nan(length(sel_cols),1);
bi_low(sel_cols>1)=bif_p_mat(sub2ind(size(bif_p_mat),find(sel_cols>1),sel_cols(sel_cols>1)));
bi_high=nan(length(sel_cols3),1);
bi_high(sel_cols3>0)=bif_p_mat(sub2ind(size(bif_p_mat),find(sel_cols3),sel_cols3(sel_cols3>0)));
mo_high=nan(length(sel_cols8),1);
mo_high(sel_cols8>0)=bif_p_mat(sub2ind(size(bif_p_mat),find(sel_cols8),sel_cols8(sel_cols8>0)));

trans1=mo_low;
trans2=bi_low;
trans3=bi_high;
trans4=mo_high;

trans1((sel_cols5>1)&((sel_cols-sel_cols6)>0))=bif_p_mat(sub2ind(size(bif_p_mat),...
    find((sel_cols5>1)&((sel_cols-sel_cols6)>0)),sel_cols5((sel_cols5>1)&((sel_cols-sel_cols6)>0))));
trans2(sel_cols2>0&((sel_cols-sel_cols6)>0))=bif_p_mat(sub2ind(size(bif_p_mat),...
    find(sel_cols2>0&((sel_cols-sel_cols6)>0)),sel_cols2(sel_cols2>0&((sel_cols-sel_cols6)>0))));
trans3(sel_cols4>0&((sel_cols8-sel_cols3)>0))=bif_p_mat(sub2ind(size(bif_p_mat),...
    find(sel_cols4>0&((sel_cols8-sel_cols3)>0)),sel_cols4(sel_cols4>0&((sel_cols8-sel_cols3)>0))));
trans4(sel_cols7>0&((sel_cols8-sel_cols3)>0))=bif_p_mat(sub2ind(size(bif_p_mat),...
    find(sel_cols7>0&((sel_cols8-sel_cols3)>0)),sel_cols7(sel_cols7>0&((sel_cols8-sel_cols3)>0))));

xstat_low = nan(length(sel_cols6),1);
xstat_low(sel_cols6>0) = xa_mat(sub2ind(size(bif_p_mat),find(sel_cols6>0),sel_cols6(sel_cols6>0)))>10;
xstat_high = nan(length(sel_cols8),1);
xstat_high(sel_cols8>0) = xa_mat(sub2ind(size(bif_p_mat),find(sel_cols8>0),sel_cols8(sel_cols8>0)))>10;

out=[mo_low trans1 trans2 bi_low bi_high trans3 trans4 mo_high xstat_low xstat_high];

end