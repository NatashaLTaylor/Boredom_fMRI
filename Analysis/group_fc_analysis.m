function [group_fc] = group_fc_analysis
%% Group Functional Connectivity
load('subject_ID.mat'); %loads subject ID with ts extracted
%cd to timeseries direct
cd '/project/def-jdancker/n24taylo/HCP_neuroimaging/timeseries/'

group_fc = zeros(960,502,502); %no. subjects x nROIs

for i=1:960
    filename = sprintf('%d%s',subject_ID(i,2),'_rfMRI_REST1_LR_timeseries.mat');
    load(filename); %variable name ts
    ts_corr = corr(ts); %calculate the Pearson's correlation between regions
    group_fc(i,:,:) = ts_corr; %save corr into one group-level file
end

%save output
save('group_fc.mat','group_fc');

quit
end
