function [group_modz] = group_modz
%% Group Functional Connectivity
load('subject_ID.mat'); %loads subject ID with ts extracted
%cd to timeseries direct
cd '/project/def-jdancker/n24taylo/HCP_neuroimaging/MTD/Graph_Analysis/'

group_modz = zeros(960,502,1195); %no. subjects x nROIS xtime

for i=1:960
    filename = sprintf('%d%s',subject_ID(i,2),'_modz.mat');
    load(filename); 
    group_modz(i,:,:) = modz; %save corr into one group-level file
end

%save output
save('group_modz.mat','group_modz');

quit
end