function  [sig,pval] = perm_code(data1,data2,iter)

    %if only comparing with one dataset, just create empty vector data2

    size1 = size(data1,1);
    size2 = size(data2,1);
    
    null_delta = zeros(iter,1);
    
    data_combo = vertcat(data1,data2);
    grp_combo(1:size1,1) = 1;
    grp_combo(size1+1:size1+size2,1) = 2;
    orig_delta = nanmean(data_combo(grp_combo==1))-nanmean(data_combo(grp_combo==2));
    
    for x = 1:iter
        rand_vec = rand(size1+size2,1);
        [~,sort_rand] = sort(rand_vec);
        grp_rand = grp_combo(sort_rand); 
        null_delta(x,1) = nanmean(data_combo(grp_rand==1))-nanmean(data_combo(grp_rand==2));
    end
    
    pos_thr = prctile(null_delta,97.5);
    neg_thr = prctile(null_delta,2.5);
    
    pos_sig = double(null_delta> pos_thr);
    neg_sig = double(null_delta < neg_thr);
    both_sig = [pos_sig,neg_sig]; 

    %pval = sum(both_sig)/iter;
    % 
    % pos_pval = sum(null_delta>orig_delta)/iter;
    % neg_pval = sum(null_delta<orig_delta)/iter;
    % 
    % both_pval = [pos_pval,neg_pval];
    % pval = min(both_pval);
    pval = sum(abs(null_delta)>=abs(orig_delta))/iter;
    sig = double(pval<0.05);
end