function [group_mod] = group_modularity
%% Group Functional Connectivity
load('subject_ID.mat'); %loads subject ID with ts extracted
%cd to timeseries direct
cd '/project/def-jdancker/n24taylo/HCP_neuroimaging/MTD/Graph_Analysis/'

group_mod = zeros(960,1195); %no. subjects x time

for i=1:960
    filename = sprintf('%d%s',subject_ID(i,2),'_modularity.mat');
    load(filename); 
    group_mod(i,:) = q; %save corr into one group-level file
end

%save output
save('group_modularity.mat','group_mod');

quit
end