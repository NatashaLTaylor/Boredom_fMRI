function [mtd,mtd_flat]=mtd_analysis_supercomputer(subnum,window,direction,trim)
%"mtd_analysis_supercomputer($subnum,3,0,1)"
%% running MTD on one subject
dir = '/project/def-jdancker/n24taylo/HCP_neuroimaging/trial_MTD/';
%addpath(genpath('/project/def-jdancker/n24taylo/HCP_neuroimaging/trial_MTD/'))
%addpath /project/def-jdancker/n24taylo/HCP_neuroimaging/code/
%subnum = getenv('SLURM_JOB_ID'); %subject name from subject JOB array
filename = sprintf('%d%s',subnum,'_rfMRI_REST1_LR_timeseries.mat')
load(filename); %variable name ts
ts_corr = corr(ts); %calculate the Pearson's correlation between regions
nNodes = size(ts,2)
fig=figure;
imagesc(ts_corr)
saveas(fig,'ts_correlation.png')
%running MTD with middle calc and trimming zeros
%% Coupling inbeded in this function
%params
% window=3;%chosen from Shine et al., 2015 paper
% direction=0;
% trim=1;
   % check inputs and define variables
    
    if nargin==2
        direction=0; trim=0;
    elseif nargin==1
        trim=0;
    end

    [t,nodes] = size(ts);

    %calculate temporal derivative
    td = diff(ts);

    %standardize data
    data_std = std(td);

    for i = 1:nodes
         td(:,i) = td(:,i) / data_std(1,i);
    end


    % [...] = zscore(X,FLAG,DIM) standardizes X by working along the dimension
    %    DIM of X. Pass in FLAG==0 to use the default normalization by N-1, or 1 to use N.


    %functional coupling score
    fc = bsxfun(@times,permute(td,[1,3,2]),permute(td,[1,2,3]));



    %temporal smoothing (credit: T. C. O'Haver, 2008.)
    mtd_temp = zeros(nodes,nodes,t-1);
    
    for j = 1:nodes
        for k = 1:nodes
            mtd_temp(j,k,:) = smooth(squeeze(fc(:,j,k)),window);
        end
    end


    %window type (0 = middle; 1 = forward facing)

    mtd = zeros(nodes,nodes,t);

    if direction == 1
        mtd(:,:,1:t-round(window/2+1)) = mtd_temp(:,:,round(window/2+1):end);
    elseif direction == 0
        mtd(:,:,1:t-1) = mtd_temp;
    end



    %trim ends (0 = no; 1 = yes)?

    if trim == 1 && direction == 0
        mtd(:,:,t-round(window/2):end) = [];
        mtd(:,:,1:round(window/2)) = [];
    elseif trim == 1 && direction == 1
        mtd(:,:,(t-window):end) = [];
    end
size(mtd)


%flatten across mtd 
template = tril(ones(nNodes)-eye(nNodes));
for tt=1:size(mtd,3)
    mtd_temp = mtd(:,:,tt);
    mtd_flat(:,tt) = mtd_temp(template==1);
end

%save output
%addpath /project/def-jdancker/n24taylo/HCP_neuroimaging/MTD/
%cd /project/def-jdancker/n24taylo/HCP_neuroimaging/MTD/
save_filename1 = sprintf('%d%s',subnum,'_mtd.mat')
save([save_filename1],'mtd','-v7.3')
save_filename2 = sprintf('%d%s',subnum,'_flatmtd.mat')
save([save_filename2],'mtd_flat','-v7.3')

quit
end

