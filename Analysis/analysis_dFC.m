%% Comparison Analysis - Boredom

% FC Analysis 
%run anova between 3 groups of boredom response
% 3 groups being the boredom score, 0,1,2
load("group_fc.mat")
%remove those without boredom score
group_fc(all_locs,:,:) = [];
%all locs to remove [551, 560, 562, 587, 600, 601]
flat_fc = zeros(size(group_fc,1),length(template_id));
for i=1:size(group_fc,1)
    a = squeeze(group_fc(i,:,:));
    flat_fc(i,:)=a(template_id);
    clear a
end

% group-level difference in boredom score vs the average fc across all
% areas
mean_fc_per_sub = mean(flat_fc,2);
[p,tbl,stats]=anova1(mean_fc_per_sub',member_bored(:,4));

[p,tbl,stats]=anova1((group_fc'),member_bored(:,4));

%group-level difference in boredom score vs edge-wise fc for network hubs
avg_network_fc =zeros(947,7,502);
for i=1:7
    for k=1:947
    avg_network_fc(k,i,:)=mean(group_fc(k,voltron_order(:,2)==i,:),2);
    end
end
%plotted out for each network grp
for i=1:7%across each network are there differences
    [p_network(i,:),tbl_network(i,:),stats_network(i,:)]=anova1(squeeze(avg_network_fc(:,i,:))',member_bored(:,4));
end

figure
fc_avg_grp_0 = mean(group_fc(member_bored(:,4)==0,:,:));
imagesc(squeeze(fc_avg_grp_0))

aov = anova(squeeze(avg_network_fc(:,2,:))',member_bored(:,4));

% dFC analysis

%create file of all std of MTD (variability of dFC)
load('subject_ID.mat') %subject IDs that we have data for

for i=602:size(subject_ID,1) %551, 560, 562, 587, 600, 601
    filename1=sprintf('%d%s',subject_ID(i,2),'_std_mtd.mat');
    filename2=sprintf('%d%s',subject_ID(i,2),'_avg_mtd.mat');
    load(filename1);
    load(filename2);
    all_std_dfc(i,:,:)=std_mtd;
    all_avg_dfc(i,:,:)=avg_mtd;
end
%these IDs don't have values [551, 560, 562, 587, 600, 601]
%remove these IDs
a = [551, 560, 562, 587, 600, 601]; %subject ID to remove (did not have their data)
b = all_avg_dfc;
b(a,:,:)=[];
remove_avg_dfc = b;
c = all_std_dfc;
c(a,:,:)=[];
remove_std_dfc = c;


%flatten out across edges
for nn=1:nROI
template = find(tril(ones(nROI))-eye(nROI)); %try taking lower triangle instead tril
end


%% dFC Analysis 
% load bored score
load('all_boredscore.mat') %1st col = IDs,  2nd = sex bin (1 M, 0 F), 3rd = age, 4th = ASRVIII.ASR_083 ('I'm easily bored')

% Find IDs that do not match boredom IDs
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

save('boredom_score_subject.mat','subject_id_bored','member_bored');

%other_locs %subject to be removed from data
remove_avg_dfc(other_locs,:,:) = [];
remove_std_dfc(other_locs,:,:) = [];

% run pearson correlation against boredom score
nROI = size(all_avg_dfc,3);
template = tril(ones(502)-eye(502)); %creates template to flatten across edges
template_id = find(template);
new_vect = a(template_id);

flat_avg_dfc = zeros(size(remove_avg_dfc,1),length(template_id));
for i=1:size(remove_avg_dfc,1)
    a = squeeze(remove_avg_dfc(i,:,:));
    flat_avg_dfc(i,:)=a(template_id);
    clear a
end

%fig avg. dfc per group
mean_dfc_per_grp = mean(remove_avg_dfc(member_bored(:,4)==2,:,:));
figure
set(gcf,'Color','w'); 
imagesc(squeeze(mean_dfc_per_grp(:,voltron_order(:,3),voltron_order(:,3))))
title('Avg. dFC for subjects ASRV 83Q =2')
ylabel('ROI')
xlabel('ROI')

%Anova dFC
p_dfc = zeros(length(flat_avg_dfc),1)
for i=1:length(flat_avg_dfc)
    [p_dfc(i,:),~,~]=anova1(flat_avg_dfc(:,i),member_bored(:,4),'off');
end
%find significance
loc_edges_dfc_sig = find(p_dfc<0.05);
sig_loc_edges = zeros(length(flat_avg_dfc,1));
sig_loc_edges(loc_edges_dfc_sig,1)==1;
matif_sig_edges = matify(sig_loc_edges,502);

avg_dfc_0 = mean(remove_avg_dfc(member_bored(:,4)==0,:,:));
avg_dfc_1 = mean(remove_avg_dfc(member_bored(:,4)==1,:,:));
avg_dfc_2 = mean(remove_avg_dfc(member_bored(:,4)==2,:,:));


%average variability for within each network (17 networks)
for i=1:17
    for k=1:947
    avg_network_dfc(k,i,:)=mean(remove_avg_dfc(k,voltron_order(:,2)==i,:),2);
    end
end




%% Variability in dFC -
% all_std_dfc % 960 x 502 x 502 
% remove subs
remove_std_dfc;

%group networks first
std_network_fc =zeros(947,7,502);

% run pearson correlation against boredom score
nROI = size(remove_std_dfc,3);
template = tril(ones(502)-eye(502)); %creates template to flatten across edges
template_id = find(template);
%new_vect = a(template_id);

flat_std_dfc = zeros(size(remove_std_dfc,1),length(template_id));
for i=1:size(remove_std_dfc,1)
    a = squeeze(remove_std_dfc(i,:,:));
    flat_std_dfc(i,:)=a(template_id);
    clear a
end

%vector id for the networks that are flattened for each ROI
for i=1:size(remove_std_dfc,1)
    a = squeeze(remove_std_dfc(i,:,:));
    flat_std_dfc(i,:)=a(template_id);
    clear a
end
flat_network =

%average variability for within each network (17 networks)
for i=1:17
    for k=1:947
    std_network_fc(k,i,:)=mean(subset_std_dfc(k,voltron_order(:,2)==i,:),2);
    end
end

for i=1:17
    for k=1:947
    a(k,i,:)=mean(subset_std_dfc(k,voltron_order(:,2)==i,voltron_order(:,2)==i),2);
    end
end
%anova of variability of dFC for the 3 groups
[p,tbl,stats]=anova1(subset_std_dfc(:,1,1),member_bored(:,4)); %this would do it for each


for i=1:17
    [p,tbl,stats]=anova1(std_network_fc(:,i,:),member_bored(:,4)); %this would do it for each
end





