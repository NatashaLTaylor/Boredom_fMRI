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
figure
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
colormap(OrangeBlueColorMap)
xticks([])
yticks([])
colorbar

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