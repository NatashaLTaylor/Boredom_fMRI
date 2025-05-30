function [left, right] = RB_surf_schaef(data,name)
%SURF_SCHAEF       Plots results onto Schaefer Atlas
%
% Note: data must be 400x1

struc_L = gifti('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/Code/Conte69.L.midthickness.32k_fs_LR.surf.gii');
struc_R = gifti('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/Code/Conte69.R.midthickness.32k_fs_LR.surf.gii');
left = gifti('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/Code/schaef_left.func.gii');
right = gifti('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/Code/schaef_right.func.gii');

% identity matrix
idx = zeros(401,2);
idx(2:401,1) = 1:1:400;
idx(2:401,2) = data;

%load('C:/Users/natas/Documents/PhD/Code/RB_midWhite_ColorMap.mat')
%C = RB_midWhite;
load('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/PhD/Code/OrangeBlueColorMap.mat')
C = OrangeBlueColorMap;

l = max(abs(data));
limits = [-l, l];


% left
original = left.cdata;
[~,index_net] = ismember(original,idx(:,1));
map_net = idx(:,2);
left.cdata = map_net(index_net);
figure; plot(struc_L,left);
set(gcf,'color','w');
caxis(limits);
view(-90,10);
colormap(C)


% Medial view
figure; plot(struc_L,left);
caxis(limits);
view(90,10);
colormap(C)
set(gcf,'color','w');
lgt = findobj(gcf,'Type','Light');
lgt = lgt(1);
[az,el] = lightangle(lgt);
lightangle(lgt,20,el);

if nargin == 2
    abc = sprintf('%s%s%s','save(left,''left_',name,'.func.gii'',''Base64Binary'');');
    eval(abc)
end


% right
original = right.cdata;
[~,index_net] = ismember(original,idx(:,1));
map_net = idx(:,2);
right.cdata = map_net(index_net);
figure; plot(struc_R,right);
set(gcf,'color','w');
caxis(limits);
view(90,10);
colormap(C)

lgt = findobj(gcf,'Type','Light');
lgt = lgt(1);
[az,el] = lightangle(lgt);
lightangle(lgt,60,el);
% set(p3,'AmbientStrength',0.2);
% set(p3,'SpecularStrength',0.4);


% Medial view
figure; plot(struc_R,right);
set(gcf,'color','w');
caxis(limits);
view(-90,10);
colormap(C)
lgt = findobj(gcf,'Type','Light');
lgt = lgt(1);
[az,el] = lightangle(lgt);
lightangle(lgt,240,el);

if nargin == 2
    abc = sprintf('%s%s%s','save(right,''right_',name,'.func.gii'',''Base64Binary'');');
    eval(abc)
end


end

