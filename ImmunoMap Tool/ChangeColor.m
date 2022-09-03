% This code will change the top 3 colors in ColorsSegments for weighted
% dendrogram
load('ColorsSegments.mat')
%% Define colors
color1 = [248 45 22]/256; % Color for blood
color2 = [242 242 89]/256; % Color for tumor
color3 = [116 33 247]/256; % Color for cross reactive TCRs
%% Substitute color
colorspec1(1).spec = color1;
colorspec1(2).spec = color2;
colorspec1(3).spec = color3;
colorspec2(1).spec = [color1 0.5];
colorspec2(2).spec = [color2 0.5];
colorspec2(3).spec = [color3 0.5];
%% Save file
save('ColorsSegments.mat','colorspec1','colorspec2')