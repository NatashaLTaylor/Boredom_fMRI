function [group_flex] = group_flexibility
%% Group Functional Connectivity
load('subject_ID.mat'); %loads subject ID with ts extracted
%cd to timeseries direct
cd '/project/def-jdancker/n24taylo/HCP_neuroimaging/MTD/Graph_Analysis/'

group_flex = zeros(960,502); %no. subjects x nROIS

for i=1:960
    filename = sprintf('%d%s',subject_ID(i,2),'_flex.mat');
    load(filename); 
    group_flex(i,:) = f; %save corr into one group-level file
end

%save output
save('group_flexibility.mat','group_flex');

quit
end