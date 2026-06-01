%% FC Analysis 2 grps:

load('FC_grp_analysis.mat')

grp0_fc = flat_fc2(member_bored2(:,4)==0,:);
grp2_fc = flat_fc2(member_bored2(:,4)==2,:);
% 1. permuted significant difference

% flatten and loop through
nROIs = 502;
for nn=1:nROIs
    template = find(tril(ones(nROIs))-eye(nROIs)); %try taking lower triangle instead tril
end
nEdges = size(grp0_fc,2);

% run permutation across fc edges - this is the correct way! - 9/09/24
for i=1:nEdges
    [sig_fc_edges(i,:),pval_fc_edges(i,:)]=perm_code(grp0_fc(:,i),grp2_fc(:,i),1000);
end

%significant diff. from permutation
sig_fc_grp0 = sig_fc_edges.*mean(grp0_fc,1)';
sig_fc_grp2 = sig_fc_edges.*mean(grp2_fc,1)';

mat_sig_fc_grp0 = matify(sig_fc_grp0,502);
mat_sig_fc_grp2 = matify(sig_fc_grp2,502);

diff_mat_sig_fc_grp2_grp0 = mat_sig_fc_grp2 - mat_sig_fc_grp0;
%reorder into network grps
reorder_diff_sig_fc_grp2_grp0 = diff_mat_sig_fc_grp2_grp0(voltron_order(:,3),voltron_order(:,3));
csvwrite("sig_diff_fc_grp2_grp0.csv",reorder_diff_sig_fc_grp2_grp0);


%% Box plots for network level analysis - revised 29/05/2026

%load in 8 order networks schaef
%first column is them unordered 11 networks, 8 cortical, 9=cerebellum, 10 = subcortex,
%11 = AAS, second column is to reorder them, 3rd column is them reordered

% indices:
%order_8_nets(:,1)==1; % vis Cent + Peri
% order_8_nets(:,1)==2; %somat mot a/b
% order_8_nets(:,1)==3; %dorsatten
%order_8_nets(:,1)==4; %sal ventatten
%order_8_nets(:,1)==5; % limbic
%order_8_nets(:,1)==6; %all Cont A-C
%order_8_nets(:,1)==7; %all default nets
%order_8_nets(:,1)==8; %temp par net
%order_8_nets(:,1)==9; %cerebellum
%order_8_nets(:,1)==10; % subcortex
%order_8_nets(:,1)==11; %AAS nuclei


load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/schaef_order8_networks.mat')

% Reshape sig edges into matrix form to identify which connections are significant
sig_fc_mat = matify(sig_fc_edges, 502);  % binary 502x502 matrix of significant edges

% Define network labels
net_labels = {'Visual','Somato-Motor','Dorsal Attention','Salience Attention','Limbic','Cognitive Control','Default','Tempo-parietal','Cerebellum','Subcortex','Ascending Arousal'};
nNets = 11;

% Define which ROIs belong to which network
net_idx = order_8_nets(:,1);  % network assignment for each of 502 ROIs

% Choose which network pairs you want to examine
% Example: within-network averages (diagonal blocks)
% You can also do between-network pairs

% --- Compute per-subject average FC for significant connections within each network block ---

grp0_avg = [];  % subjects x networks
grp2_avg = [];

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
        grp0_avg(:, n) = mean(grp0_fc(:, sig_block_vec), 2);
        grp2_avg(:, n) = mean(grp2_fc(:, sig_block_vec), 2);
    else
        grp0_avg(:, n) = NaN(size(grp0_fc, 1), 1);
        grp2_avg(:, n) = NaN(size(grp2_fc, 1), 1);
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
    
    data_grp0 = grp0_avg(:, n);
    data_grp2 = grp2_avg(:, n);
    
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
    bp = boxplot(all_data, group_label, 'Labels', {'Never', 'Often'}, 'Widths', 0.6, 'Symbol', '');
    
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
    
    % --- SD error bars ---
    % Compute mean and SD across subjects for each group
    mean_grp0 = nanmean(data_grp0);
    sd_grp0   = nanstd(data_grp0);
    mean_grp2 = nanmean(data_grp2);
    sd_grp2   = nanstd(data_grp2);
    
    % x-positions: boxplot places groups at x = 1 and x = 2
    errorbar(1, mean_grp0, sd_grp0, 'k', 'LineWidth', 1, 'Marker', 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', col_grp0, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    errorbar(2, mean_grp2, sd_grp2, 'k', 'LineWidth', 1, 'Marker', 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', col_grp2, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    
    title(net_labels{n}, 'FontSize', 9);
    ylabel('Mean FC (sig. edges)');
    set(gca, 'FontSize', 8);
end

exportgraphics(fig_sig_fc,'/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/Manuscript_CABN/Figures/Avg_FC_Net_Sig_Grp0_Grp2.svg','ContentType','vector');

%% Generate network plots - ie., Visual network to rest of the cortex etc.. for sig difference edges

% Make the difference matrix symmetric so we can index freely
diff_mat = diff_mat_sig_fc_grp2_grp0;
diff_mat = diff_mat + diff_mat';  % symmetrise (diagonal is already 0)

net_idx = order_8_nets(:,1);
nNets = 11;
nROIs = 502;

net_to_region_avg = zeros(nNets, nROIs);

for n = 1:nNets
    % Get the ROIs belonging to this network
    net_rois = find(net_idx == n);
    
    % For each of the 502 regions, average across the network's ROIs
    % diff_mat(net_rois, :) gives a submatrix of size [nNetROIs x 502]
    % Taking the mean across rows gives a 1 x 502 vector
    net_to_region_avg(n, :) = mean(diff_mat(net_rois, :), 1);
end


data = net_to_region_avg(6,1:400);
RB_surf_schaef_boredom(data,'Cog_Avg_Diff')






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
        fc_grp0 = grp0_fc(:, sig_block_vec);
        fc_grp2 = grp2_fc(:, sig_block_vec);
        
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
        grp0_subj_var(:, n) = NaN(size(grp0_fc, 1), 1);
        grp2_subj_var(:, n) = NaN(size(grp2_fc, 1), 1);
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
    
    bp = boxplot(all_data, group_label, 'Labels', {'Never', 'Often'}, 'Widths', 0.6);
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

exportgraphics(fig_fc_sd,'/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/Manuscript_CABN/Figures/VariabilityFC_Sig_Grp0_Grp2.svg','ContentType','vector');


%% P-Values for within network plots - ttest2 - independent samples

% Recreate the template for the lower triangle (to match your flattened FC arrays)
template = find(tril(ones(nROIs), -1)); 

% --- Compute per-subject average FC for ALL connections within each network ---
grp0_avg = [];  
grp2_avg = [];

for n = 1:nNets
    roi_mask = (net_idx == n);
    
    % Create a mask for this network block and FORCE it to be logical
    block_mask = logical(roi_mask * roi_mask'); % <--- FIXED LINE
    block_mask = tril(block_mask, -1);          % lower triangle only, no diagonal
    
    % Flatten to edge vector using the template
    block_vec = block_mask(template);           % logical vector of length nEdges
    
    nEdgesInNetwork = sum(block_vec);
    
    if nEdgesInNetwork > 0
        % Average across ALL edges in this network for each subject
        grp0_avg(:, n) = mean(grp0_fc(:, block_vec), 2);
        grp2_avg(:, n) = mean(grp2_fc(:, block_vec), 2);
    else
        grp0_avg(:, n) = NaN(size(grp0_fc, 1), 1);
        grp2_avg(:, n) = NaN(size(grp2_fc, 1), 1);
    end
end

%% Generate box plots, SD error bars, and p-values

% with the individual data points plotted
col_grp0 = [193, 170, 211] / 255;  % light purple (Never)
col_grp2 = [20, 90, 50] / 255;     % dark green (Often)

fig2_nonsig_fc = figure('Position', [100 100 1400 600], 'Theme', 'light');
set(gcf, 'color', 'w')

% Initialise results table
results = table();

for n = 1:nNets
    subplot(2, 6, n); hold on;
    
    data_grp0 = grp0_avg(:, n);
    data_grp2 = grp2_avg(:, n);
    
    % Skip if no edges exist in this network
    if all(isnan(data_grp0))
        title(net_labels{n}, 'FontSize', 9);
        text(0.5, 0.5, 'No edges', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        results.Network{n}  = net_labels{n};
        results.p_value(n)   = NaN;
        results.t_stat(n)    = NaN;
        results.df(n)        = NaN;
        results.cohens_d(n)  = NaN;
        results.mean_grp0(n) = NaN;
        results.mean_grp2(n) = NaN;
        results.sd_grp0(n)   = NaN;
        results.sd_grp2(n)   = NaN;
        continue;
    end
    
    % --- Calculate Statistics ---
    [~, p_val, ~, stats] = ttest2(data_grp0, data_grp2);
    
    % --- Effect Size (Cohen's d, pooled SD) ---
    n0 = sum(~isnan(data_grp0));
    n2 = sum(~isnan(data_grp2));
    mean_grp0 = nanmean(data_grp0);
    sd_grp0   = nanstd(data_grp0);
    mean_grp2 = nanmean(data_grp2);
    sd_grp2   = nanstd(data_grp2);
    pooled_sd = sqrt(((n0-1)*sd_grp0^2 + (n2-1)*sd_grp2^2) / (n0 + n2 - 2));
    cohens_d  = (mean_grp0 - mean_grp2) / pooled_sd;
    
    % --- Save to results table ---
    results.Network{n}  = net_labels{n};
    results.p_value(n)   = p_val;
    results.t_stat(n)    = stats.tstat;
    results.df(n)        = stats.df;
    results.cohens_d(n)  = cohens_d;
    results.mean_grp0(n) = mean_grp0;
    results.mean_grp2(n) = mean_grp2;
    results.sd_grp0(n)   = sd_grp0;
    results.sd_grp2(n)   = sd_grp2;
    
    % --- Build display strings ---
    if p_val < 0.001
        p_str = 'p < 0.001***';
    elseif p_val < 0.01
        p_str = sprintf('p = %.3f**', p_val);
    elseif p_val < 0.05
        p_str = sprintf('p = %.3f*', p_val);
    else
        p_str = sprintf('p = %.3f', p_val);
    end
    d_str = sprintf('d = %.2f', cohens_d);
    
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
    errorbar(1, mean_grp0, sd_grp0, 'k', 'LineWidth', 2, 'Marker', 'o', ...
        'MarkerSize', 3.5, 'MarkerFaceColor', col_grp0, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    errorbar(2, mean_grp2, sd_grp2, 'k', 'LineWidth', 2, 'Marker', 'o', ...
        'MarkerSize', 3.5, 'MarkerFaceColor', col_grp2, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    
    % Add title with network name, p-value, and effect size
    title({net_labels{n}, p_str, d_str}, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Mean FC (all edges)');
    set(gca, 'FontSize', 8);
end

% --- Display and save results ---
disp(results);
writetable(results, 'ttest_results_with_effectsize.csv');

exportgraphics(fig2_nonsig_fc,'/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/Manuscript_CABN/Figures/Avg_FC_NonSig_withinnet_sig_Grp0_Grp2.svg','ContentType','vector');




%% --- Compute FC all groups (SD across subjects) for significant edges within each network ---

grp1_id =find(member_bored(:,4)==1);

grp1_fc = flat_fc(grp1_id,:);

% --- 1. Setup and Data Extraction ---
% Assuming grp0_fc, grp1_fc, and grp2_fc are already loaded (Subjects x Edges)
grp0_subj_avg = [];  
grp1_subj_avg = [];
grp2_subj_avg = [];

net_idx = net_idx(:);
all_mask = true(502, 1);

for n = 1:nNets
    roi_mask = (net_idx == n);
    
    % Mask for connections from network n to ALL regions (no sig filter)
    block_mask = (roi_mask * all_mask') | (all_mask * roi_mask');
    block_mask = tril(block_mask, -1);
    block_vec = block_mask(template);
    
    if sum(block_vec) > 0
        grp0_subj_avg(:, n) = nanmean(grp0_fc(:, block_vec), 2);
        grp1_subj_avg(:, n) = nanmean(grp1_fc(:, block_vec), 2);
        grp2_subj_avg(:, n) = nanmean(grp2_fc(:, block_vec), 2);
    else
        grp0_subj_avg(:, n) = NaN(size(grp0_fc, 1), 1);
        grp1_subj_avg(:, n) = NaN(size(grp1_fc, 1), 1);
        grp2_subj_avg(:, n) = NaN(size(grp2_fc, 1), 1);
    end
end

% --- 2. Plotting & Stats Setup ---
col_grp0 = [193, 170, 211] / 255;
col_grp1 = sscanf('A8D5BA', '%2x%2x%2x')' / 255; % Hex #A8D5BA to RGB
col_grp2 = [20, 90, 50] / 255;

fig3_nonsig_fc = figure('Position', [100 100 1500 800], 'Theme', 'light');
set(gcf, 'color', 'w')

results = table();

% --- 3. Main Loop: Stats, Table, and Plotting ---
for n = 1:nNets
    subplot(2, 6, n); hold on;
    
    data_grp0 = grp0_subj_avg(:, n);
    data_grp1 = grp1_subj_avg(:, n);
    data_grp2 = grp2_subj_avg(:, n);
    
    % Skip if no edges exist
    if all(isnan(data_grp0))
        title(net_labels{n}, 'FontSize', 9);
        text(0.5, 0.5, 'No edges', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        results.Network{n} = net_labels{n};
        continue; 
    end
    
    % --- Calculate Statistics (Pairwise t-tests) ---
    [~, p_0v1, ~, stats_0v1] = ttest2(data_grp0, data_grp1);
    [~, p_1v2, ~, stats_1v2] = ttest2(data_grp1, data_grp2);
    [~, p_0v2, ~, stats_0v2] = ttest2(data_grp0, data_grp2);
    
    % --- Means, SDs, and Cohen's d ---
    n0 = sum(~isnan(data_grp0));
    n1 = sum(~isnan(data_grp1));
    n2 = sum(~isnan(data_grp2));
    
    m0 = nanmean(data_grp0); s0 = nanstd(data_grp0);
    m1 = nanmean(data_grp1); s1 = nanstd(data_grp1);
    m2 = nanmean(data_grp2); s2 = nanstd(data_grp2);
    
    % Pooled SDs and Cohen's d
    pool_0v1 = sqrt(((n0-1)*s0^2 + (n1-1)*s1^2) / (n0 + n1 - 2));
    d_0v1 = (m0 - m1) / pool_0v1;
    
    pool_1v2 = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2));
    d_1v2 = (m1 - m2) / pool_1v2;
    
    pool_0v2 = sqrt(((n0-1)*s0^2 + (n2-1)*s2^2) / (n0 + n2 - 2));
    d_0v2 = (m0 - m2) / pool_0v2;
    
    % --- Save to Results Table ---
    results.Network{n} = net_labels{n};
    results.mean_grp0(n) = m0; results.sd_grp0(n) = s0;
    results.mean_grp1(n) = m1; results.sd_grp1(n) = s1;
    results.mean_grp2(n) = m2; results.sd_grp2(n) = s2;
    
    results.p_0v1(n) = p_0v1; results.t_0v1(n) = stats_0v1.tstat; results.d_0v1(n) = d_0v1;
    results.p_1v2(n) = p_1v2; results.t_1v2(n) = stats_1v2.tstat; results.d_1v2(n) = d_1v2;
    results.p_0v2(n) = p_0v2; results.t_0v2(n) = stats_0v2.tstat; results.d_0v2(n) = d_0v2;
    
    % --- Boxplot ---
    all_data = [data_grp0; data_grp1; data_grp2];
    group_label = [ones(length(data_grp0), 1); 2 * ones(length(data_grp1), 1); 3 * ones(length(data_grp2), 1)];
    
    bp = boxplot(all_data, group_label, 'Labels', {'Never', 'Sometimes', 'Often'}, 'Widths', 0.6);
    set(findobj(gca, 'Tag', 'Whisker'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Lower Adjacent Value'), 'Visible', 'off');
    set(findobj(gca, 'Tag', 'Upper Adjacent Value'), 'Visible', 'off');
    
    % Colour the boxes (MATLAB box handles are in reverse order)
    h = findobj(gca, 'Tag', 'Box');
    if length(h) >= 3
        patch(get(h(3), 'XData'), get(h(3), 'YData'), col_grp0, 'FaceAlpha', 0.5);
        patch(get(h(2), 'XData'), get(h(2), 'YData'), col_grp1, 'FaceAlpha', 0.5);
        patch(get(h(1), 'XData'), get(h(1), 'YData'), col_grp2, 'FaceAlpha', 0.5);
    end
    
    % --- Individual Subject Data Points (jittered) ---
    jitter_width = 0.15;
    
    scatter(1 + jitter_width*(rand(length(data_grp0),1)-0.5), data_grp0, 20, col_grp0, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.2, 'MarkerFaceAlpha', 0.4);
    scatter(2 + jitter_width*(rand(length(data_grp1),1)-0.5), data_grp1, 20, col_grp1, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.2, 'MarkerFaceAlpha', 0.4);
    scatter(3 + jitter_width*(rand(length(data_grp2),1)-0.5), data_grp2, 20, col_grp2, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.2, 'MarkerFaceAlpha', 0.4);
    
    % --- SD error bars ---
    errorbar(1, m0, s0, 'k', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', col_grp0, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    errorbar(2, m1, s1, 'k', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', col_grp1, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    errorbar(3, m2, s2, 'k', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', col_grp2, 'MarkerEdgeColor', 'k', 'CapSize', 10);
    
    % --- Format Title Strings ---
    % Helper function to add stars for significance
    get_stars = @(p) repmat('*', 1, (p<0.05) + (p<0.01) + (p<0.001));
    
    %str_0v1 = sprintf('0v1: p=%.3f%s, d=%.2f', p_0v1, get_stars(p_0v1), d_0v1);
    %str_1v2 = sprintf('1v2: p=%.3f%s, d=%.2f', p_1v2, get_stars(p_1v2), d_1v2);
    %str_0v2 = sprintf('0v2: p=%.3f%s, d=%.2f', p_0v2, get_stars(p_0v2), d_0v2);
    
    title(net_labels{n}, 'FontSize', 8, 'FontWeight', 'bold');
    ylabel('Mean FC (all edges)');
    set(gca, 'FontSize', 8);
end

% --- 4. Display and save results ---
disp(results);
writetable(results, 'ttest_results_3groups_effectsize.csv');


exportgraphics(fig3_nonsig_fc,'/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/Manuscript_CABN/Figures/Avg_FC_NonSig_withinnet_sig_Grp0_Grp1_Grp2.svg','ContentType','vector');




%% Graph Theory Measures - between two groups

%load in
dir("*_pc.mat")
%load in grp pc
load('group_pc.mat')
% 0 in spots where pc didnt or the time window was smaller
% only 947 ID's with pc
load("member_bored.mat")
% subject with "NaN" removed - 
locs = isnan(member_bored(:,4));
member_bored_update = member_bored;
member_bored_update(locs==1,:)=[];

% separate into grps for PC
grp0_pc = group_pc(member_bored(:,4)==0,:,:);
grp2_pc = group_pc(member_bored(:,4)==2,:,:);

% Check for rows that are all zeros across all slices
all_zero_rows = all(all(grp2_pc == 0, 2), 3);

% Find the indices of rows that are all zeros
row_indices = find(all_zero_rows);

% Display the result
if ~isempty(row_indices)
    fprintf('The following rows are all zeros across all slices:\n');
    disp(row_indices');
else
    fprintf('No rows are all zeros across all slices.\n');
end

% comparison analysis between two

%avg pc
pc_avg_0 = mean(grp0_pc,3);
pc_avg_2 = mean(grp2_pc,3);

% Perform the Kolmogorov-Smirnov test
[h, p] = kstest((pc_avg_0 - mean(pc_avg_0)) / std(pc_avg_0));

if h == 0
    fprintf('The data is normally distributed (p = %.3f).\n', p);
else
    fprintf('The data is not normally distributed (p = %.3f).\n', p);
end

% perform permutation test - as data is not normally distributed
%% Graph Theory Analysis between grp0 grp2
load('workspace_updated.mat')

part_grp_0 = part_grp_0';
part_grp_1 = part_grp_1';
part_grp_2 = part_grp_2';
for i = 1:502
    [sig_pc_edges(i,:),pval_pc_edges(i,:)]=perm_code(pc_avg_0(:,i),pc_avg_2(:,i),1000);
end

mean_pc_avg_0 = mean(pc_avg_0);
mean_pc_avg_2 = mean(pc_avg_2);

diff_sig_pc_grp2_grp0 = sig_pc_edges'.*(mean_pc_avg_2- mean_pc_avg_0);

%normall distributed p-value
p_rank_mean = ranksum(mean_pc_avg_0, mean_pc_avg_2);









%% Figures
load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/HCP_boredom_backup/Code/Github_code/Figures/OrangeBlueColorMapNew.mat');
%C = OrangeBlueColorMap;
C = NewCustomColormap;

f1=figure('Theme','light')
set(gcf,'color','w')
imagesc(diff_mat_sig_fc_grp2_grp0(voltron_order(:,3),voltron_order(:,3)))
title('Sig Permuted FC Bored - Healthy')
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

exportgraphics(f1,'/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/7. Boredom- Danckert/Manuscript_CABN/Figures/Avg_FC_Sig_Grp0_Grp2_matrix.svg','ContentType','vector');


figure
set(gcf,'color','w')
subplot(1,2,1)
imagesc(mat_sig_fc_grp0(voltron_order(:,3),voltron_order(:,3)))
title('Sig Permuted FC  Healthy')
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
colormap(OrangeBlueColorMap)
xticks([])
yticks([])
colorbar
subplot(1,2,2)
imagesc(mat_sig_fc_grp2(voltron_order(:,3),voltron_order(:,3)))
title('Sig Permuted FC  Bored')
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
colormap(OrangeBlueColorMap)
xticks([])
yticks([])
colorbar

%try get an average - 
mean_diff_sig_fc_grp2_grp0 = mean(diff_mat_sig_fc_grp2_grp0,)


%% figure - participation coeff differences
figure
set(gcf,'color','w')
imagesc(diff_sig_pc_grp2_grp0(voltron_order(1:400,3)))
title('Sig Permuted Avg. PC Bored - Healthy')
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
colormap(OrangeBlueColorMap)
xticks([])
yticks([])
colorbar

% brain plot - 
sig_avg_pc_0 = sig_pc_edges.*mean_pc_avg_0.';

sig_avg_pc_2 = sig_pc_edges.*mean_pc_avg_2.';

diff_sig_pc_grp2_grp0=diff_sig_pc_grp2_grp0.';

RB_surf_schaef(diff_sig_pc_grp2_grp0(1:400,:),'diff_sig_pc_grp2_grp0')

RB_surf_schaef(sig_avg_pc_0(1:400,:),'sig_pc_grp0')

RB_surf_schaef(sig_avg_pc_2(1:400,:),'sig_pc_grp2')