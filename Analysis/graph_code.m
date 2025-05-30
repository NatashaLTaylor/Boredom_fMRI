function [ci,q,part,modz,hc,f] = graph_code(subnum)
%% Graph Theoretical Analysis
%INTEGRATION       Creates time-resolved network topological measures
% function [ci,q,part,modz,hc,f] = integration_plus5(data,gamma,beta)
%  [ci,q,p,z,hc,f] = integration_plus1(data,gamma,beta);
%
%  This code takes a time-resolved connectivity matrix and 
%  estimates community structure, modularity and the cartographic profile
%  for each region within a region x region x time connectivity matrix.
%  See https://arxiv.org/abs/1511.02976 for more details.
%
%  Requirements: https://sites.google.com/site/bctnet/
%
%  Input:      data     time-series organized in 'nodes x nodes x time' matrix
%              gamma    tuning parameter for louvain algorithm (low = large modules & high = small modules)
%              beta     similarity measure used to determine clustering assignment (hungarian algorithm) - requires munkres.m & apcluster.m
%
%  Output:     ci       time-resolved community assignment
%              q        time-resolved modularity
%              p        time-resolved participation coefficient
%              z        time-resolved module-degree z-score
%              hc       cartographic profile
%              f        flexibility

    %% Load variables
filename = sprintf('%d%s',subnum,'_mtd.mat')
data=load(filename);

% loadfile = sprintf('%s%d%s','load(''',subnum,'_mtd.mat'')');
% eval(loadfile); %loads the mtd file
gamma = 1;
beta = 0.75;

[nodes,~,time] = size(data);


community_louvain;
participation_coef_sign;
module_degree_zscore;
    %% Nested function community_louvain
    function [M,Q]=community_louvain(W,gamma,M0,B)
        W=double(W);                                % convert to double format
        n=length(W);                                % get number of nodes
        s=sum(sum(W));                              % get sum of edges
    
    if ~exist('B','var') || isempty(B)
        type_B = 'modularity';
    elseif ischar(B)
        type_B = B;
    else
        type_B = 0;
        if exist('gamma','var') && ~isempty(gamma)
            warning('Value of gamma is ignored in generalized mode.')
        end
    end
    if ~exist('gamma','var') || isempty(gamma)
        gamma = 1;
    end
    
    if strcmp(type_B,'negative_sym') || strcmp(type_B,'negative_asym')
        W0 = W.*(W>0);                          %positive weights matrix
        s0 = sum(sum(W0));                      %weight of positive links
        B0 = W0-gamma*(sum(W0,2)*sum(W0,1))/s0; %positive modularity
        
        W1 =-W.*(W<0);                          %negative weights matrix
        s1 = sum(sum(W1));                      %weight of negative links
        if s1                                   %negative modularity
            B1 = W1-gamma*(sum(W1,2)*sum(W1,1))/s1;
        else
            B1 = 0;
        end
    elseif min(min(W))<-1e-10
        err_string = [
            'The input connection matrix contains negative weights.\nSpecify ' ...
            '''negative_sym'' or ''negative_asym'' objective-function types.'];
        error(sprintf(err_string))              
    end
    if strcmp(type_B,'potts') && any(any(W ~= logical(W)))
        error('Potts-model Hamiltonian requires a binary W.')
    end
    
    if type_B
        switch type_B
            case 'modularity';      B = (W-gamma*(sum(W,2)*sum(W,1))/s)/s;
            case 'potts';           B =  W-gamma*(~W);
            case 'negative_sym';    B = B0/(s0+s1) - B1/(s0+s1);
            case 'negative_asym';   B = B0/s0      - B1/(s0+s1);
            otherwise;              error('Unknown objective function.');
        end
    else                            % custom objective function matrix as input
        B = double(B);
        if ~isequal(size(W),size(B))
            error('W and B must have the same size.')
        end
    end
    if ~exist('M0','var') || isempty(M0)
        M0=1:n;
    elseif numel(M0)~=n
        error('M0 must contain n elements.')
    end
    
    [~,~,Mb] = unique(M0);
    M = Mb;
    
    B = (B+B.')/2;                                          % symmetrize modularity matrix
    Hnm=zeros(n,n);                                         % node-to-module degree
    for m=1:max(Mb)                                         % loop over modules
        Hnm(:,m)=sum(B(:,Mb==m),2);
    end
    
    Q0 = -inf;
    Q = sum(B(bsxfun(@eq,M0,M0.')));                        % compute modularity
    first_iteration = true;
    while Q-Q0>1e-10
        flag = true;                                        % flag for within-hierarchy search
        while flag
            flag = false;
            for u=randperm(n)                               % loop over all nodes in random order
                ma = Mb(u);                                 % current module of u
                dQ = Hnm(u,:) - Hnm(u,ma) + B(u,u);
                dQ(ma) = 0;                                 % (line above) algorithm condition
                
                [max_dQ,mb] = max(dQ);                      % maximal increase in modularity and corresponding module
                if max_dQ>1e-10                             % if maximal increase is positive
                    flag = true;
                    Mb(u) = mb;                             % reassign module
                    
                    Hnm(:,mb) = Hnm(:,mb)+B(:,u);           % change node-to-module strengths
                    Hnm(:,ma) = Hnm(:,ma)-B(:,u);
                end
            end
        end
        [~,~,Mb] = unique(Mb);                              % new module assignments
        
        M0 = M;
        if first_iteration
            M=Mb;
            first_iteration=false;
        else
            for u=1:n                                       % loop through initial module assignments
                M(M0==u)=Mb(u);                             % assign new modules
            end
        end
        
        n=max(Mb);                                          % new number of modules
        B1=zeros(n);                                        % new weighted matrix
        for u=1:n
            for v=u:n
                bm=sum(sum(B(Mb==u,Mb==v)));                % pool weights of nodes in same module
                B1(u,v)=bm;
                B1(v,u)=bm;
            end
        end
        B=B1;
        
        Mb=1:n;                                             % initial module assignments
        Hnm=B;                                              % node-to-module strength
        
        Q0=Q;
        Q=trace(B);                                         % compute modularity
    end
    end

    function [Ppos,Pneg]=participation_coef_sign(W,Ci)

            n=length(W);                                %number of vertices
    
    Ppos = pcoef( W.*(W>0));
    Pneg = pcoef(-W.*(W<0));
    
        function P=pcoef(W_)
            S   = sum(W_,2);                    %strength
            Gc  = (W_~=0)*diag(Ci);             %neighbor community affiliation
            Sc2 = zeros(n,1);                   %community-specific neighbors
            
            for i = 1:max(Ci)
                Sc2 = Sc2 + (sum(W_.*(Gc==i),2).^2);
            end
            
            P = ones(n,1) - Sc2./(S.^2);
            P(isnan(P)) = 0;
            P(~P) = 0;                            %p_ind=0 if no (out)neighbors
        end
    end

    function Z=module_degree_zscore(W,Ci,flag)
         if ~exist('flag','var')
            flag=0;
        end
    
        switch flag
            case 0  % no action required
            case 1  % no action required
            case 2; W=W.';
            case 3; W=W+W.';
        end
    
        n=length(W);                        %number of vertices
        Z=zeros(n,1);
            for i=1:max(Ci)
                Koi=sum(W(Ci==i,Ci==i),2);
                Z(Ci==i)=(Koi-mean(Koi))./std(Koi);
            end
            Z(isnan(Z))=0;
    end

    
    


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


    

    %% cartographic profile

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

quit
end


