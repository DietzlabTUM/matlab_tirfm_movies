% This and other related scripts were written by Dr. Matthias Schickinger, Dr. Philip Ketterer, Dr. Jonas Funke,
% Anna-Katharina Pumm, Dr. Wouter Engelen and Eva Bertosin

%% Read .tif files

folder_name = uigetdir;
cd(folder_name);

fnames = dir('*.tif');
fnamesfits = dir('*.fits');
numfids = length(fnames);
vals = cell(1,numfids);
vals_norm = cell(1,numfids);

for K = 1:numfids
    vals{K} = imread(fnames(K).name);
    vals_norm{K} = uint16(255*mat2gray(vals{K}));
    K;
end

disp('Done')

%% Classify images and st_dev with gui

i=1;

go_on = 1;

while go_on
    
    disp(['To go : ' num2str(length(fnames)-i)]);
    
    plot_spot_sum_and_stddev_nocb01(i,vals,dev_img, fnames, numfids)
    
    button_stationary = uicontrol('Style', 'pushbutton', 'String', 'Stationary', 'Position', [150 20 100 35], 'Callback', 'class(i) = 1;');
    button_blobb = uicontrol('Style', 'pushbutton', 'String', 'Diffusive', 'Position', [300 20 100 35], 'Callback', 'class(i) = 2;');
    button_moving = uicontrol('Style', 'pushbutton', 'String', 'Moving', 'Position', [450 20 100 35], 'Callback', 'class(i) = 3;');
    button_unclear = uicontrol('Style', 'pushbutton', 'String', 'Unclear', 'Position', [600 20 100 35], 'Callback', 'class(i) = 4;');
    button_back = uicontrol ('Style', 'pushbutton', 'String', 'Back', 'Position', [1080 20 100 35], 'Callback', 'class(i) = 5');
    
    while waitforbuttonpress
    end
    
    pause(0.3)
    
    if isempty(class(i))
        go_on = 0;
    elseif class(i)==0 % Repeat
        disp('Repeat')
        i = i-1;
    elseif class(i)==5
        i = i-2;
    else
        disp([num2str(i) ' classified'])
    end
    
    if i==length(fnames)
        go_on = 0;
        disp('end')
    else
        i=i+1;
    end
    
    close all;
end

close all;


%% Sort classified spots, move to folders etc.

mkdir('stationary');
mkdir('diffusive');
mkdir('moving');
mkdir('unclear');

% enter the file copying thing
for N = 1:numfids
    
    switch class(N)
        case 1
            movefile(fnames(N).name,['stationary/' fnames(N).name])
            movefile(fnamesfits(N).name,['stationary/' fnamesfits(N).name])
        case 2
            movefile(fnames(N).name,['diffusive/' fnames(N).name])
            movefile(fnamesfits(N).name,['diffusive/' fnamesfits(N).name])
        case 3
            movefile(fnames(N).name,['moving/' fnames(N).name])
            movefile(fnamesfits(N).name,['moving/' fnamesfits(N).name])
        case 4
            movefile(fnames(N).name,['unclear/' fnames(N).name])
            movefile(fnamesfits(N).name,['unclear/' fnamesfits(N).name])
    end
    
end


disp('Spots moved')

%%
counter = zeros(4);

for M = 1:numfids
    if class(M) == 1
        counter(1) = counter(1) + 1;
    elseif class(M) == 2
        counter(2) = counter(2) +1;
    elseif class(M) == 3
        counter(3) = counter(3) + 1;
    elseif class(M) == 4
        counter(4) = counter(4) + 1;
    end
end

bar(counter);
title('spot categories: stationary, diffusive, moving, unclear')
disp(counter)