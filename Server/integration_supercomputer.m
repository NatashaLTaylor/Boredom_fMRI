function [ci,q,part,modz,hc,f] = integration_supercomputer(subnum,gamma,beta)
%INTEGRATION       Creates time-resolved network topological measures
%  [ci,q,p,z,hc,f] = integration_plus1(data,gamma,beta);
%
%  Requirements: https://sites.google.com/site/bctnet/
%
%  Input:      data     time-series organized in 'nodes x nodes x time' matrix
%              gamma(1)    tuning parameter for louvain algorithm (low = large modules & high = small modules)
%              beta(0.75)     similarity measure used to determine clustering assignment (hungarian algorithm) - requires munkres.m & apcluster.m
%
%  Output:     ci       time-resolved community assignment
%              q        time-resolved modularity
%              part        time-resolved participation coefficient
%              modz        time-resolved module-degree z-score
%              hc       cartographic profile
%              f        flexibility

%directories
addpath(genpath('/project/def-jdancker/n24taylo/HCP_neuroimaging/code/graph_code/'));
%addpath '/project/def-jdancker/n24taylo/HCP_neuroimaging/graph_code/' %path to all functions
dir = '/project/def-jdancker/n24taylo/HCP_neuroimaging/MTD/';
cd '/project/def-jdancker/n24taylo/HCP_neuroimaging/MTD/'
    %define variables
    filename = sprintf('%d%s',subnum,'_mtd.mat');
    load(filename);
    data = mtd;

    [nodes,~,time] = size(data);

    ci = zeros(nodes,time); q = zeros(time,1); part = zeros(nodes,time); modz = zeros(nodes,time);

    for t = 1:time 
        
        if t == 1
            [ci(:,t),q(t,1)] = community_louvain(data(:,:,t),gamma,1:1:nodes,'negative_asym');
            part(:,t) = participation_coef_sign(data(:,:,t),ci(:,1));
            modz(:,t) = module_degree_zscore(data(:,:,t),ci(:,1));
        else
            [ci(:,t),q(t,1)] = community_louvain(data(:,:,t),gamma,ci(:,t-1),'negative_asym');
            part(:,t) = participation_coef_sign(data(:,:,t),ci(:,t));
            modz(:,t) = module_degree_zscore(data(:,:,t),ci(:,1));
        end
    end

    %cartographic profile

    xbins = 0:0.01:1; ybins = 5:-.1:-5;
    hc = zeros(size(xbins,2),size(ybins,2),time);
    xNumBins = numel(xbins); yNumBins = numel(ybins);

    for t = 1:time
        Xi = round(interp1(xbins,1:xNumBins,part(:,t),'linear','extrap'));
        Yi = round(interp1(ybins,1:yNumBins,modz(:,t),'linear','extrap'));
        Xi = max(min(Xi,xNumBins),1);
        Yi = max(min(Yi,yNumBins),1);
        hc(:,:,t) = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
    end

    
    nDer = time - 1;
    
%     dyn_mod = zeros(nMod,time);
    number_mod = zeros(time,1);

    for t = 1:time
        temp = tabulate(ci(:,t));
        dyn_mod(1:size(temp,1),t) = temp(:,1);
        number_mod(t,1) = nnz(dyn_mod(:,t)); %number of modules per window
    end    
   
    ignore = double(number_mod>5);
    nMod = max(max(ci(:,ignore==0)));
    dyn_mod(nMod+1:end,:) = [];
    
    %1-of-k encoding
    encode = zeros(nodes,nMod,time);

    for t = 1:time
        if ignore(t,1)==0
            C = ci(:,t);
            R = 1:numel(C);
            A = zeros(numel(C),max(C));
            A(sub2ind(size(A),R',C)) = 1;
            encode(:,1:size(A,2),t) = A;
        else
            encode(:,:,t) = 0;
        end
    end

   
    %dice between consecutive time points
    dice_coef_encode = zeros(nMod,nMod,time);

    for t = 1:nDer
        dice_coef_encode(:,:,t) = bsxfun(@corr,encode(:,:,t),encode(:,:,t+1));
    end


    %threshold & cost
    cost = 1/(double(dice_coef_encode>beta));

    %hungarian algorithm
    assignment = zeros(time,nMod);

    for t = 1:time
        [assignment(t,:),~] = munkres(cost(:,:,t));
    end

    %hungarian un-twisting 
    dyn_mod2 = zeros(nMod,time);  

    for t = 1:time
        for k = 1:nMod
            if number_mod(t,1) < k
                dyn_mod2(k,t) = NaN;
            end
        end
    end

    dyn_mod2(:,1,:) = dyn_mod(:,1,:); %dyn_mod2 starting with dyn_mod's first assignment


    % recoding
    tally = max(dyn_mod2(:,1));            

    for w = 2:time-1
        for k = 1:nMod
            for l = 1:nMod
                if dyn_mod2(k,w-1) == 0
                    dyn_mod2(k,w-1) = tally+1;
                    tally = tally+1;
                end

                if assignment(w-1,k) == l
                    dyn_mod2(l,w) = dyn_mod2(k,w-1);
                end
            end
        end
    end

    ci_temp = zeros(nodes,time);

    for t = 1:time
        for j = 1:nodes
            for k = 1:nMod
                if ci(j,t) == dyn_mod(k,t)
                    ci_temp(j,t) = dyn_mod2(k,t);
                end
            end
        end
    end

    

    %%temporal sorting

    %ci_new_name_1_of_k
    encode2 = zeros(nodes,tally,time);

    for t = 1:time
        for j = 1:nodes
            for h = 1:tally
                if ci_temp(j,t) == h
                    encode2(j,h,t) = 1;
                end
            end
        end
    end


    %ci_signature
    encode2_sum = nanmean(encode2,3);
    encode2_perc = encode2_sum / tally;


    %ci_correlation matrix -- can we use a better similarity metric here?
    corr_sig = corr(encode2_perc);

    
    %use affinity propogation to cluster
    [ap,~,~,~] = apcluster(corr_sig,nanmin(corr_sig(:)));

    
    
    %identity of clusters estimated by affinity propogation
    unique_ap = unique(ap);
    
    %rename ci_new_name according to ap clusters
    dyn_mod3 = dyn_mod2;

    for t = 1:time
        for k = 1:nMod
            for h = 1:tally
                if dyn_mod2(k,t) == h
                    dyn_mod3(k,t) = ap(h,1);
                end
            end
        end
    end

    
    % rename parcels according to affinity clustering
    ci_temp2 = zeros(nodes,time);
    ci_new = zeros(nodes,time);
    

    for t = 1:time
        for j = 1:nodes
            for k = 1:nMod
                if ci_temp(j,t) == dyn_mod2(k,t)
                    ci_temp2(j,t) = dyn_mod3(k,t);
                end
            end
        end
    end
    
    
    %% rename parcels according to order of appearance
    
    for t = 1:time
        for j = 1:nodes
            for x = 1:size(unique_ap,1)
                y = unique_ap(x);
                if ci_temp2(j,t) == y
                    ci_new(j,t) = find(unique_ap==y);
                end
            end
        end
    end
    
    ci = ci_new;
        
    f = flexibility(ci_new');
    
    %% save outputs
    cd '/project/def-jdancker/n24taylo/HCP_neuroimaging/MTD/Graph_Analysis/'
    save_filename1 = sprintf('%d%s',subnum,'_pc.mat');
    save([save_filename1],'part')
    save_filename2 = sprintf('%d%s',subnum,'_modularity.mat');
    save([save_filename2],'q')
    save_filename3 = sprintf('%d%s',subnum,'_modz.mat');
    save([save_filename3],'modz')
    save_filename4 = sprintf('%d%s',subnum,'_flex.mat');
    save([save_filename4],'f')
    save_filename5 = sprintf('%d%s',subnum,'_cart.mat');
    save([save_filename5],'hc')

    quit
end