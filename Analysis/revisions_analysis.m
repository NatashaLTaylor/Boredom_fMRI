%% Boredom Analysis - revision analysis 23/05/2026

DATA_DIR="/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/HCP_boredom_backup";
cd(DATA_DIR)
% load in boredom score and subject_ids
% load bored score
load(fullfile(DATA_DIR,'Behaviour_var/all_boredscore.mat')) %1st col = IDs,  2nd = sex bin (1 M, 0 F), 3rd = age, 4th = ASRVIII.ASR_083 ('I'm easily bored')

% Find IDs that do not match boredom IDs
load('subject_ID.mat')
a = [551, 560, 562, 587, 600, 601]; %subject ID to remove (did not have their data)
subject_data_ID = subject_ID;
subject_data_ID(a,:)=[];
not_subject = ismember(Boredscore(:,1),subject_data_ID(:,2));
member_data = Boredscore(:,1).*not_subject;
locs_remove = find(member_data==0);

b = ismember(subject_data_ID(:,2),Boredscore(:,1));
other_locs = find(b==0); %these need to be removed from the data avg/std, as we don't have their behav data?

member_bored = Boredscore;
member_bored(locs_remove,:)=[]; %boredom score of subjects (for the data we have so far)
subject_id_bored = subject_data_ID;
subject_id_bored(other_locs,:) = []; %remove subject with missing behav data from imaging data




% create subject_id + pc file

% load in pc 
load("group_pc.mat");


% create subject_id + flexibility file


% create subject_id + fc (over all 400 ROIs)





% this is the original - grp 0 and grp 2 separation
% (FC_analysis_grp0_grp2.mat)
load("FC_grp_analysis.mat"); %loads them all in fc data
% join with the IDs
subject_id_bored = member_bored(:,1);
group_fc_data_flat = [subject_id_bored,flat_fc];
% save under: 
writematrix(group_fc_data_flat,'orig_grp0_2_fc_data_flat.csv');




%------------------------------------%
%% FOLLOW-UP: ALL GRP FC COMPARISON
%------------------------------------%
load("/Analysis/analysis_workspace.mat")
% using spearman ranks:

%grp 0='never bored'
%grp 1 = 'sometimes bored'
%grp 2 = 'often bored'

score = member_bored(:,4);
score = score(:);
nEdges = size(flat_fc, 2); %flatten across all ROI edges
flat_fc(idx,:)=[];

%% Preallocate
rho_vals = nan(nEdges,1);
p_vals   = nan(nEdges,1);

%% Spearman correlation for each edge
for e = 1:nEdges

    edge_values = flat_fc(:,e);

    % Remove NaNs if present
    valid_idx = ~isnan(edge_values) & ~isnan(score);

    if sum(valid_idx) > 5

        [rho,p] = corr( ...
            edge_values(valid_idx), ...
            score(valid_idx), ...
            'Type','Spearman');

        rho_vals(e) = rho;
        p_vals(e)   = p;
    end
end
%% Non-multiple comparison correction:
sig_idx_nofdr = p_vals<0.05;
fprintf('Number significant edges: %d\n', sum(sig_idx_nofdr));

sig_nofdr_rho = double(sig_idx_nofdr).*rho_vals;
%% Multiple comparison correction (FDR)
p_fdr = mafdr(p_vals, 'BHFDR', true);

%% Significant edges
sig_idx = p_fdr < 0.05;

fprintf('Number significant edges: %d\n', sum(sig_idx));

%% Positive relationships only
positive_sig = sig_idx & rho_vals > 0;

fprintf('Positive significant edges: %d\n', sum(positive_sig));

%% Example outputs
sig_rho = rho_vals(positive_sig);
sig_p   = p_fdr(positive_sig);


%% Avg FC for each subject

idx = find(isnan(member_bored(:,4)));

fc_z_flat = atanh(flat_fc);

avg_fc = nanmean(fc_z_flat,2);
avg_fc(idx)=[];
score(idx)=[];

[rho,p] = corr(avg_fc,score,'Type','Spearman');
fprintf('Spearman rho = %.3f\n', rho);
fprintf('p = %.5f\n', p);


fc_all_flat = [subject_id_bored,flat_fc];
writematrix(fc_all_flat,'fc_subject_all.csv');


%% Average each ROI edge for FC - all subjects

nROIs = 502;

% Replicate your exact indexing
template = find(tril(ones(nROIs)) - eye(nROIs));
[row_idx, col_idx] = ind2sub([nROIs, nROIs], template);

% a single subject's edge_vector ---
% roi_mean = zeros(nROIs, 1);  % 502 x 1 output
% for r = 1:nROIs
%     edge_mask = (row_idx == r) | (col_idx == r);
%     roi_mean(r) = mean(edge_vector(edge_mask));
% end

% across all subject:
template = find(tril(ones(nROIs)) - eye(nROIs));
[row_idx, col_idx] = ind2sub([nROIs, nROIs], template);

roi_membership = zeros(length(template), nROIs);
for r = 1:nROIs
    roi_membership(:, r) = (row_idx == r) | (col_idx == r);
end

% Result: [nSubjects x 502]
roi_means = flat_fc * roi_membership ./ sum(roi_membership, 1);

fc_means_roi_all = [subject_id_bored,roi_means];
writematrix(fc_means_roi_all,"fc_subjects_all_mean_roi.csv");



%% Variance plots for FC for each group (boxplot)

nROIs = 502;
template = find(tril(ones(nROIs)) - eye(nROIs));
[row_idx, col_idx] = ind2sub([nROIs, nROIs], template);

% edge_data is [nSubjects x nEdges]
nSubjects = size(flat_fc, 1);
roi_std_fc_all = zeros(nSubjects, nROIs);  % [nSubjects x 502]

for s = 1:nSubjects
    for r = 1:nROIs
        edge_mask = (row_idx == r) | (col_idx == r);
        roi_std_fc_all(s, r) = std(flat_fc(s, edge_mask));
    end
end




%% Box plots for differences in FC avg for the networks - across all 3 groups :

flat_fc %all subjects data
member_bored(:,4) %contains the 0, 1 ,2 variables for the grouping


% with the individual data points plotted
col_grp0 = [193, 170, 211] / 255;  % light purple (Never)
col_grp1 = [168 213 186] /255; % light green (Sometimes)
col_grp2 = [20, 90, 50] / 255;     % dark green (Often)

fig2_nonsig_fc =figure('Position', [100 100 1400 600], 'Theme', 'light');
set(gcf, 'color', 'w')

for n = 1:nNets
    subplot(2, 6, n); hold on;
    
    data_grp0 = grp0_avg(:, n);
    data_grp1 = grp1_avg(:,n);
    data_grp2 = grp2_avg(:, n);
    
    % Skip if no edges exist in this network
    if all(isnan(data_grp0))
        title(net_labels{n}, 'FontSize', 9);
        text(0.5, 0.5, 'No edges', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        continue;
    end
    
    % --- Calculate Statistics ---
    [~, p_val, ~, stats] = ttest2(data_grp0, data_grp2);
    
    if p_val < 0.001
        p_str = 'p < 0.001***';
    elseif p_val < 0.01
        p_str = sprintf('p = %.3f**', p_val);
    elseif p_val < 0.05
        p_str = sprintf('p = %.3f*', p_val);
    else
        p_str = sprintf('p = %.3f', p_val);
    end
    
    % --- Boxplot ---
    all_data = [data_grp0; data_grp2];
    group_label = [ones(length(data_grp0), 1); 2 * ones(length(data_grp2), 1)];
    
    bp = boxplot(all_data, group_label, 'Labels', {'Never', 'Often'}, 'Widths', 0.6);
    set(findobj(gca, 'Tag', 'Whisker'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Lower Adjacent Value'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Upper Adjacent Value'), 'Visible', 'off');
    
    % Colour the boxes
    h = findobj(gca, 'Tag', 'Box');
    if length(h) >= 2
        patch(get(h(2), 'XData'), get(h(2), 'YData'), col_grp0, 'FaceAlpha', 0.5);
        patch(get(h(1), 'XData'), get(h(1), 'YData'), col_grp2, 'FaceAlpha', 0.5);
    end
    
    % --- Individual Subject Data Points (jittered) ---
    jitter_width = 0.15;
    
    jitter_grp0 = 1 + jitter_width * (rand(length(data_grp0), 1) - 0.5);
    scatter(jitter_grp0, data_grp0, 20, col_grp0, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.2, 'MarkerFaceAlpha', 0.4)
    
    jitter_grp2 = 2 + jitter_width * (rand(length(data_grp2), 1) - 0.5);
    scatter(jitter_grp2, data_grp2, 20, col_grp2, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.2, 'MarkerFaceAlpha', 0.4)
    
    % --- SD error bars ---
    mean_grp0 = nanmean(data_grp0);
    sd_grp0   = nanstd(data_grp0);
    mean_grp2 = nanmean(data_grp2);
    sd_grp2   = nanstd(data_grp2);
    
    errorbar(1, mean_grp0, sd_grp0, 'k', 'LineWidth', 2, 'Marker', 'o', ...
        'MarkerSize', 3.5, 'MarkerFaceColor', col_grp0, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    errorbar(2, mean_grp2, sd_grp2, 'k', 'LineWidth', 2, 'Marker', 'o', ...
        'MarkerSize', 3.5, 'MarkerFaceColor', col_grp2, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    
    % Add title with network name and p-value
    title({net_labels{n}, p_str}, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Mean FC (all edges)');
    set(gca, 'FontSize', 8);
end

exportgraphics(fig3_nonsig_fc_allgrp,'/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/Manuscript_CABN/Figures/Avg_FC_NonSig_AllGrp.svg','ContentType','vector');


% variability for all 3 groups
% --- Compute variability for all 3 groups ---
grp0_sd_avg = NaN(1, nNets);
grp1_sd_avg = NaN(1, nNets);
grp2_sd_avg = NaN(1, nNets);

grp0_subj_var = [];
grp1_subj_var = [];
grp2_subj_var = [];

for n = 1:nNets
    roi_mask = (net_idx == n);
    
    block_mask = roi_mask * roi_mask';
    block_mask = tril(block_mask, -1);
    sig_block = sig_fc_mat & block_mask;
    sig_block_vec = sig_block(template);
    
    nSigEdges = sum(sig_block_vec);
    
    if nSigEdges > 0
        fc_grp0 = grp0_fc(:, sig_block_vec);
        fc_grp1 = grp1_fc(:, sig_block_vec);
        fc_grp2 = grp2_fc(:, sig_block_vec);
        
        grp0_sd_avg(n) = mean(nanstd(fc_grp0, 0, 1));
        grp1_sd_avg(n) = mean(nanstd(fc_grp1, 0, 1));
        grp2_sd_avg(n) = mean(nanstd(fc_grp2, 0, 1));
        
        grp0_mean_fc = nanmean(fc_grp0, 1);
        grp1_mean_fc = nanmean(fc_grp1, 1);
        grp2_mean_fc = nanmean(fc_grp2, 1);
        
        grp0_subj_var(:, n) = mean(abs(fc_grp0 - grp0_mean_fc), 2);
        grp1_subj_var(:, n) = mean(abs(fc_grp1 - grp1_mean_fc), 2);
        grp2_subj_var(:, n) = mean(abs(fc_grp2 - grp2_mean_fc), 2);
    else
        grp0_subj_var(:, n) = NaN(size(grp0_fc, 1), 1);
        grp1_subj_var(:, n) = NaN(size(grp1_fc, 1), 1);
        grp2_subj_var(:, n) = NaN(size(grp2_fc, 1), 1);
    end
end

% --- Plotting ---
col_grp0 = [193, 170, 211] / 255;  % light purple (Never)
col_grp1 = [168, 213, 186] / 255;  % light green (Sometimes)
col_grp2 = [20, 90, 50] / 255;     % dark green (Often)

col_grp0_light = 1 - 0.5 * (1 - col_grp0);
col_grp1_light = 1 - 0.5 * (1 - col_grp1);
col_grp2_light = 1 - 0.5 * (1 - col_grp2);

fig_fc_sd_allgrp = figure('Position', [100 100 1400 600], 'Theme', 'light');
set(gcf, 'color', 'w')

for n = 1:nNets
    subplot(2, 6, n); hold on;
    
    data_grp0 = grp0_subj_var(:, n);
    data_grp1 = grp1_subj_var(:, n);
    data_grp2 = grp2_subj_var(:, n);
    
    ylabel('FC Variability');
    set(gca, 'FontSize', 8);
    
    if all(isnan(data_grp0))
        title(net_labels{n}, 'FontSize', 9);
        text(0.5, 0.5, 'No edges', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        axis normal;
        ylim([0 0.5]);
        continue;
    end
    
    % --- Boxplot (3 groups) ---
    all_data = [data_grp0; data_grp1; data_grp2];
    group_label = [ones(length(data_grp0), 1); 2 * ones(length(data_grp1), 1); 3 * ones(length(data_grp2), 1)];
    
    bp = boxplot(all_data, group_label, 'Labels', {'Never', 'Sometimes', 'Often'}, 'Widths', 0.6);
    set(findobj(gca, 'Tag', 'Whisker'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Lower Adjacent Value'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Upper Adjacent Value'), 'Visible', 'off');
    
    % Colour the boxes
    h = findobj(gca, 'Tag', 'Box');
    if length(h) >= 3
        patch(get(h(3), 'XData'), get(h(3), 'YData'), col_grp0, 'FaceAlpha', 0.5);
        patch(get(h(2), 'XData'), get(h(2), 'YData'), col_grp1, 'FaceAlpha', 0.5);
        patch(get(h(1), 'XData'), get(h(1), 'YData'), col_grp2, 'FaceAlpha', 0.5);
    end
    
    % --- Individual subject points (jittered) ---
    jitter_width = 0.15;
    
    jitter_grp0 = 1 + jitter_width * (rand(length(data_grp0), 1) - 0.5);
    scatter(jitter_grp0, data_grp0, 20, col_grp0_light, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
    
    jitter_grp1 = 2 + jitter_width * (rand(length(data_grp1), 1) - 0.5);
    scatter(jitter_grp1, data_grp1, 20, col_grp1_light, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
    
    jitter_grp2 = 3 + jitter_width * (rand(length(data_grp2), 1) - 0.5);
    scatter(jitter_grp2, data_grp2, 20, col_grp2_light, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
    
    % --- SD error bars ---
    mean_grp0 = nanmean(data_grp0);
    sd_grp0   = nanstd(data_grp0);
    mean_grp1 = nanmean(data_grp1);
    sd_grp1   = nanstd(data_grp1);
    mean_grp2 = nanmean(data_grp2);
    sd_grp2   = nanstd(data_grp2);
    
    errorbar(1, mean_grp0, sd_grp0, 'k', 'LineWidth', 1, 'Marker', 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', col_grp0, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    errorbar(2, mean_grp1, sd_grp1, 'k', 'LineWidth', 1, 'Marker', 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', col_grp1, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    errorbar(3, mean_grp2, sd_grp2, 'k', 'LineWidth', 1, 'Marker', 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', col_grp2, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    
    % Title with network name only (no significance)
    title(net_labels{n}, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('FC Variability');
    set(gca, 'FontSize', 8);
    
    axis normal;
    ylim([0 0.5]);
end

exportgraphics(fig_fc_sd_allgrp, '/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/Manuscript_CABN/Figures/VariabilityFC_Sig_AllGrp.svg', 'ContentType', 'vector');




%------------------------------------%
%% FOLLOW-UP: ALL GRP PC COMPARISON
%------------------------------------%

group_pc_time = horzcat(subject_ID,group_pc);

group_pc_avg = squeeze(mean(group_pc,3));
pc_all_subjects = [subject_ID,group_pc_avg];
writematrix(pc_all_subjects,'pc_avg_subject_all.csv');

load("workspace_updated.mat")

% average PC across the 3 groups
part_grp_all = [part_grp_0,part_grp_1,part_grp_2]';
% generate score figure
score=zeros(941,1);
score(1:453,1)=0;
score(454:454+397,1)=1;
score(851:end,1)=2;

nEdges=502;
%% Preallocate
rho_vals = nan(nEdges,1);
p_vals   = nan(nEdges,1);

%% Spearman correlation for each edge
for e = 1:nEdges

    edge_values = part_grp_all(:,e);

    % Remove NaNs if present
    valid_idx = ~isnan(edge_values) & ~isnan(score);

    if sum(valid_idx) > 5

        [rho,p] = corr( ...
            edge_values(valid_idx), ...
            score(valid_idx), ...
            'Type','Spearman');

        rho_vals(e) = rho;
        p_vals(e)   = p;
    end
end
%% Non-multiple comparison correction:
sig_idx_nofdr = p_vals<0.05;
fprintf('Number significant edges: %d\n', sum(sig_idx_nofdr));

sig_nofdr_rho = double(sig_idx_nofdr).*rho_vals;
%% Multiple comparison correction (FDR)
p_fdr = mafdr(p_vals, 'BHFDR', true);

%% Significant edges
sig_idx = p_fdr < 0.05;

fprintf('Number significant edges: %d\n', sum(sig_idx));

%% Positive relationships only
positive_sig = sig_idx & rho_vals > 0;

fprintf('Positive significant edges: %d\n', sum(positive_sig));

%------------------------------------%
%% FOLLOW-UP: ALL GRP FLEXIBILITY COMPARISON
%------------------------------------%
load('group_flexibility.mat');
load('flex_table.mat');

flex_subject_id=dataTable.Subject;

flex_data_subject = [flex_subject_id,group_flex];
% save under: 
writematrix(flex_data_subject,'flexibility_subject_all.csv');


score = str2double(string(dataTable.Group));
nEdges=502;
%% Preallocate
rho_vals = nan(nEdges,1);
p_vals   = nan(nEdges,1);

%% Spearman correlation for each edge
for e = 1:nEdges

    edge_values = group_flex(:,e);

    % Remove NaNs if present
    valid_idx = ~isnan(edge_values) & ~isnan(score);

    if sum(valid_idx) > 5

        [rho,p] = corr( ...
            edge_values(valid_idx), ...
            score(valid_idx), ...
            'Type','Spearman');

        rho_vals(e) = rho;
        p_vals(e)   = p;
    end
end
%% Non-multiple comparison correction:
sig_idx_nofdr = p_vals<0.05;
fprintf('Number significant edges: %d\n', sum(sig_idx_nofdr));

sig_nofdr_rho = double(sig_idx_nofdr).*rho_vals;
%% Multiple comparison correction (FDR)
p_fdr = mafdr(p_vals, 'BHFDR', true);

%% Significant edges
sig_idx = p_fdr < 0.05;

fprintf('Number significant edges: %d\n', sum(sig_idx));

%% Positive relationships only
positive_sig = sig_idx & rho_vals > 0;

fprintf('Positive significant edges: %d\n', sum(positive_sig));






%% LOAD IN LOGISTIC REGRESSION VALUES & PLOT THEM ON BRAIN

readmatrix("/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney\(Staff\)/PhD/7.\ Boredom-\ Danckert/HCP_boredom_backup/model/ordinal_boredom/ordinal_boredom_polr_functional_connectivity_never_vs_often_502roi_fdr.csv ");



%% BRAIN PLOTS REDO DIFFERENT SCALES

data = diff_mat_sig_fc_grp2_grp0(1:400,1)