% Gerald Page
% April 11, 2016

%% Leaf Area measurements on single-leaf TIR experiment scans

clear
clc
close all

% load images
drs='.\data\leaf_scans'; % in current directory
images = dir([drs '\*.png']); images = {images.name}; % file names to cell

for k = 2:length(images) % for each image
    files{k} = imread([drs '\' images{k}]);
end

for i = 1:length(files)
    bw{i} = im2bw(files{i});
end

% clear('files');

% remove scale box from corner
for i = 1:length(files)
    bw{i}(2630:3500, 1:750) = 1; % remove scale box
end

% invert bw
for i = 1:length(files)
    bw_inv{i} = 1 - bw{i};
end

% convert to area using pixel size
for i = 1:length(files)
    area(i,2) = sum(sum(bw_inv{i}));
    area(i,3) = (area(i,2)./(11.8462*11.8462))./1000000;
end

csvwrite('single_leaf_areas.csv', area)




