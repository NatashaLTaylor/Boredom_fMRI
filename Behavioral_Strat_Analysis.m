%% BEHAVIORAL SEPARATION ANALYSIS

% load in the clinical separation for ASR_Attn_T & DSM_ADH_T_Score
% Read the data into a table for analysis
clinicalData = readtable('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/HCP_boredom_backup/model/ordinal_boredom/clinical_behavioral_binary_df_never_vs_often_functional_connectivity.csv');


%% ----- FC COMPARISON ----------%

load('FC_grp_analysis.mat');

%grp based on the clinical high and low
behavior_var='DSM_Adh_T_clinical_binary';
behavior_var='ASR_Attn_T_clinical_binary';
% Extract the functional connectivity data for each group
grp_ids = clinicalData.ASR_Attn_T_clinical_binary; % Assuming 'grp0_fc' is a column in clinicalData

%flat_fc2 %fc we want - only the extreme boredom subjects

grp_clinic_fc = flat_fc2(grp_ids==1,:);
grp_nonclinic_fc = flat_fc2(grp_ids==0,:);

% 1. permuted significant difference
% flatten and loop through
nROIs = 502;
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end
nEdges = size(grp_clinic_fc,2);

% run permutation across fc edges - this is the correct way! - 9/09/24
for i=1:nEdges
    [sig_fc_edges(i,:),pval_fc_edges(i,:)]=perm_code(grp_clinic_fc(:,i),grp_nonclinic_fc(:,i),1000);
end

%significant diff. from permutation
sig_fc_grp_clinic = sig_fc_edges.*mean(grp_clinic_fc,1)';
sig_fc_grp_nonclinic = sig_fc_edges.*mean(grp_nonclinic_fc,1)';

mat_sig_fc_grp_clinic = matify(sig_fc_grp_clinic,502);
mat_sig_fc_grp_nonclinic = matify(sig_fc_grp_nonclinic,502);

%clinical minus non-clinical
diff_mat_sig_fc_grp_clinic_nonclinic = mat_sig_fc_grp_clinic - mat_sig_fc_grp_nonclinic;
%reorder into network grps
reorder_diff_sig_fc_grp_clinic_nonclinic = diff_mat_sig_fc_grp_clinic_nonclinic(voltron_order(:,3),voltron_order(:,3));

filename_save = sprintf("sig_diff_fc_grp_clinic_nonclinic_%s.csv",behavior_var);
csvwrite(filename_save,reorder_diff_sig_fc_grp_clinic_nonclinic);


%plot box plot for each network
load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/schaef_order8_networks.mat')
net_labels = {'Visual','Somato-Motor','Dorsal Attention','Salience Attention','Limbic','Cognitive Control','Default','Tempo-parietal','Cerebellum','Subcortex','Ascending Arousal'};
nNets = 11;

% Define which ROIs belong to which network
net_idx = order_8_nets(:,1);  % network assignment for each of 502 ROIs

% Choose which network pairs you want to examine
% Example: within-network averages (diagonal blocks)
% You can also do between-network pairs

% --- Compute per-subject average FC for significant connections within each network block ---
sig_fc_mat = matify(sig_fc_edges,502);
grpclinic_avg = [];  % subjects x networks
grpnonclinic_avg = [];

for n = 1:nNets
    roi_mask = (net_idx == n);
    
    % Create a mask for this network block (within-network connections)
    block_mask = roi_mask * roi_mask';          % 502x502 logical
    block_mask = tril(block_mask, -1);          % lower triangle only, no diagonal
    
    % Intersect with significant edges
    sig_block = sig_fc_mat & block_mask;        % only significant connections in this block
    
    % Flatten back to edge vector using the same template
    sig_block_vec = sig_block(template);        % logical vector of length nEdges
    
    nSigEdges = sum(sig_block_vec);
    
    if nSigEdges > 0
        % Average across only the significant edges for each subject
        grpclinic_avg(:, n) = mean(grp_clinic_fc(:, sig_block_vec), 2);
        grpnonclinic_avg(:, n) = mean(grp_nonclinic_fc(:, sig_block_vec), 2);
    else
        grpclinic_avg(:, n) = NaN(size(grp_clinic_fc, 1), 1);
        grpnonclinic_avg(:, n) = NaN(size(grp_nonclinic_fc, 1), 1);
    end
end

%% Generate box plots with SD error bars for each network

% Colours
col_grp0 = [193, 170, 211] / 255;  % light purple (Never)
col_grp2 = [20, 90, 50] / 255;     % dark green (Often)

fig_sig_fc=figure('Position', [100 100 1400 500],'Theme','light');
set(gcf,'color','w')
for n = 1:nNets
    subplot(2, 6, n); hold on;
    
    data_grp2 = grpclinic_avg(:, n); %grp2 like 'never' is the 'clinical group'
    data_grp0 = grpnonclinic_avg(:, n);
    
    % Skip if no significant connections in this network
    if all(isnan(data_grp0))
        title(net_labels{n}, 'FontSize', 9);
        text(0.5, 0.5, 'No sig. edges', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        continue;
    end
    
    % Combine for grouped boxplot
    all_data = [data_grp0; data_grp2];
    group_label = [ones(length(data_grp0), 1); 2 * ones(length(data_grp2), 1)];
       % Added 'Symbol', '' to hide default outliers so they don't clutter your SD bars
    bp = boxplot(all_data, group_label, 'Labels', {'Non-Clinical', 'Clinical'}, 'Widths', 0.6, 'Symbol', '');
    
    % Hide default boxplot whiskers and caps 
    set(findobj(gca, 'Tag', 'Whisker'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Lower Adjacent Value'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Upper Adjacent Value'), 'Visible', 'off');

    % Colour the boxes
    h = findobj(gca, 'Tag', 'Box');
    if length(h) >= 2
        patch(get(h(2), 'XData'), get(h(2), 'YData'), col_grp0, 'FaceAlpha', 0.5);
        patch(get(h(1), 'XData'), get(h(1), 'YData'), col_grp2, 'FaceAlpha', 0.5);
    end
    % % --- Individual subject data points (jittered) ---
    %     jitter_width = 0.15;  % controls horizontal spread of points
    % 
    %     % Group 0 (Non-Clinical) at x = 1
    %     jitter_grp0 = 1 + jitter_width * (rand(length(data_grp0), 1) - 0.5);
    %     scatter(jitter_grp0, data_grp0, 15, col_grp0, 'filled', ...
    %         'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    % 
    %     % Group 2 (Clinical) at x = 2
    %     jitter_grp2 = 2 + jitter_width * (rand(length(data_grp2), 1) - 0.5);
    %     scatter(jitter_grp2, data_grp2, 15, col_grp2, 'filled', ...
    %         'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    % --- SD error bars ---
    % Compute mean and SD across subjects for each group
    mean_grp0 = nanmean(data_grp0);
    sd_grp0   = nanstd(data_grp0);
    mean_grp2 = nanmean(data_grp2);
    sd_grp2   = nanstd(data_grp2);
    
    % x-positions: boxplot places groups at x = 1 and x = 2
    errorbar(1, mean_grp0, sd_grp0, 'k', 'LineWidth', 1.5, 'Marker', 'o', ...
        'MarkerSize', 5, 'MarkerFaceColor', col_grp0, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    errorbar(2, mean_grp2, sd_grp2, 'k', 'LineWidth', 1.5, 'Marker', 'o', ...
        'MarkerSize', 5, 'MarkerFaceColor', col_grp2, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    
    title(net_labels{n}, 'FontSize', 9);
    ylabel('Mean FC (sig. edges)');
    set(gca, 'FontSize', 8);
end

filename=sprintf('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/Manuscript_CABN/Figures/Avg_FC_Net_Sig_Grp_Clinical_%s.svg',behavior_var);

exportgraphics(fig_sig_fc,filename,'ContentType','vector');


%% --- Compute FC variability (SD across subjects) for significant edges within each network ---

grp0_sd_avg = NaN(1, nNets);  % 1 x networks (one value per group per network)
grp2_sd_avg = NaN(1, nNets);


% Also store per-subject deviation from group mean (for individual points on the plot)
grp0_subj_var = [];  % subjects x networks
grp2_subj_var = [];

for n = 1:nNets
    roi_mask = (net_idx == n);
    
    block_mask = roi_mask * roi_mask';
    block_mask = tril(block_mask, -1);
    sig_block = sig_fc_mat & block_mask;
    sig_block_vec = sig_block(template);
    
    nSigEdges = sum(sig_block_vec);
    
    if nSigEdges > 0
        % Extract FC values for significant edges: subjects x edges
        fc_grp0 = grp_nonclinic_fc(:, sig_block_vec);
        fc_grp2 = grp_clinic_fc(:, sig_block_vec);
        
        % --- Option A: SD across subjects for each edge, then average ---
        % This gives one summary "variability" value per group
        grp0_sd_avg(n) = mean(nanstd(fc_grp0, 0, 1));  % std across subjects (dim 1), avg across edges
        grp2_sd_avg(n) = mean(nanstd(fc_grp2, 0, 1));
        
        % --- Per-subject "variability" metric ---
        % Absolute deviation of each subject from the group mean at each edge
        % then averaged across edges => one value per subject
        grp0_mean_fc = nanmean(fc_grp0, 1);  % 1 x edges (group mean per edge)
        grp2_mean_fc = nanmean(fc_grp2, 1);
        
        grp0_subj_var(:, n) = mean(abs(fc_grp0 - grp0_mean_fc), 2);  % subjects x 1
        grp2_subj_var(:, n) = mean(abs(fc_grp2 - grp2_mean_fc), 2);
    else
        grp0_subj_var(:, n) = NaN(size(fc_grp0, 1), 1);
        grp2_subj_var(:, n) = NaN(size(fc_grp2, 1), 1);
    end
end

% --- Plotting ---
col_grp0 = [193, 170, 211] / 255;
col_grp2 = [20, 90, 50] / 255;
col_grp0_light = 1 - 0.5 * (1 - col_grp0);
col_grp2_light = 1 - 0.5 * (1 - col_grp2);

fig_fc_sd=figure('Position', [100 100 1400 600], 'Theme', 'light');
set(gcf, 'color', 'w')

for n = 1:nNets
    subplot(2, 6, n); hold on;
    
    data_grp0 = grp0_subj_var(:, n);
    data_grp2 = grp2_subj_var(:, n);
    
    if all(isnan(data_grp0))
        title(net_labels{n}, 'FontSize', 9);
        text(0.5, 0.5, 'No edges', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        continue;
    end
    
    % --- Wilcoxon rank-sum test (non-parametric, good for variance measures) ---
    [p_val, ~, stats] = ranksum(data_grp0, data_grp2);
    
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
    
    bp = boxplot(all_data, group_label, 'Labels', {'Non-Clinical', 'Clinical'}, 'Widths', 0.6);
    set(findobj(gca, 'Tag', 'Whisker'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Lower Adjacent Value'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Upper Adjacent Value'), 'Visible', 'off');
    
    h = findobj(gca, 'Tag', 'Box');
    if length(h) >= 2
        patch(get(h(2), 'XData'), get(h(2), 'YData'), col_grp0, 'FaceAlpha', 0.5);
        patch(get(h(1), 'XData'), get(h(1), 'YData'), col_grp2, 'FaceAlpha', 0.5);
    end
    
    % --- Individual subject points (jittered) ---
    jitter_width = 0.15;
    
    jitter_grp0 = 1 + jitter_width * (rand(length(data_grp0), 1) - 0.5);
    scatter(jitter_grp0, data_grp0, 20, col_grp0_light, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
    
    jitter_grp2 = 2 + jitter_width * (rand(length(data_grp2), 1) - 0.5);
    scatter(jitter_grp2, data_grp2, 20, col_grp2_light, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
    
    % --- SD error bars ---
    mean_grp0 = nanmean(data_grp0);
    sd_grp0   = nanstd(data_grp0);
    mean_grp2 = nanmean(data_grp2);
    sd_grp2   = nanstd(data_grp2);
    
    errorbar(1, mean_grp0, sd_grp0, 'k', 'LineWidth', 1, 'Marker', 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', col_grp0, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    errorbar(2, mean_grp2, sd_grp2, 'k', 'LineWidth', 1, 'Marker', 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', col_grp2, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    
    title({net_labels{n}, p_str}, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('FC Variability (MAD)');
    set(gca, 'FontSize', 8);
end

filename=sprintf('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/Manuscript_CABN/Figures/VariabilityFC_Sig_Clinical_%s.svg',behavior_var);
exportgraphics(fig_fc_sd,filename,'ContentType','vector');




% PLOT FC MATRIX:
%% Figures
load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/HCP_boredom_backup/Code/Github_code/Figures/OrangeBlueColorMapNew.mat');
%C = OrangeBlueColorMap;
C = NewCustomColormap;
f1=figure('Theme','Light')
set(gcf,'color','w')
imagesc(diff_mat_sig_fc_grp_clinic_nonclinic(voltron_order(:,3),voltron_order(:,3)))
title(sprintf('Sig Permuted FC Clinic - Non-Clinic %s', behavior_var))
hold on
line([1,502],[47,47],'Color','black','LineWidth',1) %vis Cent + Peri
line([47,47],[1,502],'Color','black','LineWidth',1)
line([1,502],[117,117],'Color','black','LineWidth',1) %somat mot a/b
line([117,117],[1,502],'Color','black','LineWidth',1) %somat mot a/b
line([1,502],[169,169],'Color','black','LineWidth',1) %dorsatten
line([169,169],[1,502],'Color','black','LineWidth',1) %dorsatten
line([1,502],[220,220],'Color','black','LineWidth',1) %sal ventatten
line([220,220],[1,502],'Color','black','LineWidth',1) %sal ventatten
line([1,502],[244,244],'Color','black','LineWidth',1) % limbic
line([244,244],[1,502],'Color','black','LineWidth',1) % limbic
line([1,502],[305,305],'Color','black','LineWidth',1) %all Cont A-C
line([305,305],[1,502],'Color','black','LineWidth',1) %all Cont A-C
line([1,502],[384,384],'Color','black','LineWidth',1) %all default nets
line([384,384],[1,502],'Color','black','LineWidth',1) %all default nets
%subcortical
line([1,502],[400,400],'Color','black','LineWidth',1) %temp par net
line([400,400],[1,502],'Color','black','LineWidth',1) %temp par net
line([1,502],[414,414],'Color','black','LineWidth',1) %hippocampus +amyg
line([414,414],[1,502],'Color','black','LineWidth',1) %hippocampus +amyg
line([1,502],[430,430],'Color','black','LineWidth',1) %thal
line([430,430],[1,502],'Color','black','LineWidth',1) %thal
line([1,502],[454,454],'Color','black','LineWidth',1) %basal gang.
line([454,454],[1,502],'Color','black','LineWidth',1) %basal gang.
line([1,502],[482,482],'Color','black','LineWidth',1) %cerebellum
line([482,482],[1,502],'Color','black','LineWidth',1) %cerebellum
colormap(C)
xticks([])
yticks([])
colorbar

filename=sprintf('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/HCP_boredom_backup/Analysis_2/Sig_permute_FC_clinical_matrix_%s.svg',behavior_var);
exportgraphics(f1,filename,'ContentType','vector');

%% Plot on the Brain for networks





%% ----- PC COMPARISON ----------%

%grp based on the clinical high and low
behavior_var='DSM_Adh_T_clinical_binary';
% Extract the functional connectivity data for each group
grp_ids = clinicalData.ASR_Attn_T_clinical_binary; % Assuming 'grp0_fc' is a column in clinicalData

%determine matching members:
[lia, locb] = ismember(subject_id_bored(:,2), clinicalData.subject_id);
% lia  = logical array (true where a match exists)
% locb = index into clinicalData for each match (0 if no match)

%get PC for those variables
group_pc_grps_match = group_pc(lia==1,:,:);

%flat_fc2 %fc we want - only the extreme boredom subjects
grp_clinic_pc = group_pc_grps_match(grp_ids==1,:,:);
grp_nonclinic_pc = group_pc_grps_match(grp_ids==0,:,:);

% comparison analysis between two

%avg pc
pc_avg_0 = mean(grp_nonclinic_pc,3);
pc_avg_2 = mean(grp_clinic_pc,3);

% Perform the Kolmogorov-Smirnov test
[h, p] = kstest((pc_avg_0 - mean(pc_avg_0)) / std(pc_avg_0));

if h == 0
    fprintf('The data is normally distributed (p = %.3f).\n', p);
else
    fprintf('The data is not normally distributed (p = %.3f).\n', p);
end

% perform permutation test - as data is not normally distributed
%% Graph Theory Analysis between grp0 grp2
%load('workspace_updated.mat')

part_grp_0 = part_grp_0';
%part_grp_1 = part_grp_1';
part_grp_2 = part_grp_2';
for i = 1:502
    [sig_pc_edges(i,:),pval_pc_edges(i,:)]=perm_code(pc_avg_0(:,i),pc_avg_2(:,i),1000);
end

mean_pc_avg_0 = mean(pc_avg_0);
mean_pc_avg_2 = mean(pc_avg_2);

diff_sig_pc_grp2_grp0 = sig_pc_edges'.*(mean_pc_avg_2- mean_pc_avg_0);

%normall distributed p-value
p_rank_mean = ranksum(mean_pc_avg_0, mean_pc_avg_2);

data=diff_sig_pc_grp2_grp0(:,1:400);

RB_surf_schaef_boredom(data','ASR_Attn_PC_Clinical_grp')


%% ----- FLEXIBILITY COMPARISON ----------%

load('group_flexibility.mat');
load('flex_table.mat');


%grp based on the clinical high and low
behavior_var='DSM_Adh_T_clinical_binary';
% Extract the functional connectivity data for each group
grp_ids = clinicalData.DSM_Adh_T_clinical_binary; % Assuming 'grp0_fc' is a column in clinicalData
subject_id = dataTable.Subject;
%determine matching members:
[lia, locb] = ismember(subject_id, clinicalData.subject_id);
% lia  = logical array (true where a match exists)
% locb = index into clinicalData for each match (0 if no match)

%get PC for those variables
group_flex_grps_match = group_flex(lia==1,:,:);

%flat_fc2 %fc we want - only the extreme boredom subjects
grp_clinic_flex = group_flex_grps_match(grp_ids==1,:,:);
grp_nonclinic_flex = group_flex_grps_match(grp_ids==0,:,:);

% permuted significance
for i = 1:502
    [sig_flex_edges(i,:),pval_flex_edges(i,:)]=perm_code(grp_nonclinic_flex(:,i),grp_clinic_flex(:,i),1000);
end

mean_flex_clinic = mean(grp_clinic_flex,1)';
mean_flex_nonclinic = mean(grp_nonclinic_flex,1)';

diff_sig_flex_grp2_grp0 = sig_flex_edges.*(mean_flex_clinic- mean_flex_nonclinic);

data = diff_sig_flex_grp2_grp0(1:400,1);
RB_surf_schaef_boredom(data','DSM_ADHD_Flex_Clinical_grp')
% running anova:

[h,p,ci,stats]=ttest2(grp_clinic_flex,grp_nonclinic_flex);

