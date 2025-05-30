function [group_pc] = group_pc
%% Group Functional Connectivity
load('subject_ID.mat'); %loads subject ID with ts extracted
%cd to timeseries direct
cd '    '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/HCP_boredom_backup/Analysis''

group_pc = zeros(947,502,1195); %no. subjects x nROIs
%101107_pc.mat; 119833_pc; 140420;
for i = 1:947
    try
        % Construct filename
        filename = sprintf('%d%s', member_bored(i, 1), '_pc.mat');
        % Load the file
        load(filename); 
        
        % Check if 'part' has the correct size
        if isequal(size(part), [502, 1195]) % Adjust size to match [nROIs, timepoints] or your specific dimensions
            group_pc(i, :, :) = part; % Assign to group matrix
        else
            group_pc(i,:,:)=zeros;
            warning('Size mismatch in file: %s. Skipping...', filename);
        end
    catch ME
        % Handle file loading or other errors
        warning('Error processing file: %s. Error: %s', filename, ME.message);
    end
end

%save output
save('group_pc.mat','group_pc');

quit
end
