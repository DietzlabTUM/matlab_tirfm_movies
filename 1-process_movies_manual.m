% This and other related scripts were written by Dr. Matthias Schickinger, Dr. Philip Ketterer, Dr. Jonas Funke,
% Anna-Katharina Pumm, Dr. Wouter Engelen and Eva Bertosin

%% startup
clc, clear all, close all
path0 = cd;

%change to correct path
run('/Users/username/Documents/MOVIE_TIRFM/functions/startup.m')


%% LOAD STACK OF MOVIES green channel

channel = cell(1,1);
channel{1} = 'green';

pname=uigetdir(data_dir,'Choose the folder with all .fits files.');
files_ch1 = pickFirstFitsFiles(pname, channel{1});

N_movie = length(files_ch1);

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
mkdir(path_out)

%% LOAD STACK OF MOVIES red channel
channel{2} = 'red';
files_ch2 = pickFirstFitsFiles(pname, channel{2});

%% set parameters and generate movie classes
ch1 = cell(N_movie,1);
ch2 = cell(N_movie,1);

first = 2;
last = -1;
time_per_frame = 50;
r_find = 8;
r_integrate = 8;
N_frames = 100;

sequence_ch1 = 1;

for i=1:N_movie
    ch1{i} = movie(pname, files_ch1{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
    if length(channel) == 2
        ch2{i} = movie(pname, files_ch2{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
    end
end

%% compute average and standard deviation images of full FOV of all Green movies
avg_img_green = cell(N_movie, 2);
stdev_img_green = cell(N_movie);

for i=1:N_movie
    avg_img_green{i, 1} = ch1{i}.average_image(N_frames);
    avg_img_green{i, 3} = ch1{i}.average_image_last(N_frames); % for fitting threshold assignment, drift correction
    mov_temp = fitsread([pname filesep ch1{i}.fname{1}], 'PixelRegion',{[1 length(avg_img_green{1})], [1 length(avg_img_green{1})], [2 500] });
    stdev_img_green{i} = std(mov_temp, 0, 3);
    i
end

clear mov_temp

%% compute average and standard deviation images of full FOV of all red movies
avg_img_red = cell(N_movie, 2);
stdev_img_red = cell(N_movie);

for i=1:N_movie
    avg_img_red{i, 1} = ch2{i}.average_image(N_frames);
    avg_img_red{i, 3} = ch2{i}.average_image_last(N_frames); % for fitting threshold assignment, drift correction
    mov_temp = fitsread([pname filesep ch2{i}.fname{1}], 'PixelRegion',{[1 length(avg_img_red{1})], [1 length(avg_img_red{1})], [2 500] });
    stdev_img_red{i} = std(mov_temp, 0, 3);
    i
end

clear mov_temp
%% manually select all rotors (green channel)
mag = 4; %magnification for easy picking

all_positions = cell(N_movie, 1);

peaks_raw = zeros(0,5);

for i=1:N_movie;
    
    [pos_ch1, N_peaks] = manual_particle_picking(i, stdev_img_green, mag);
    
    all_positions{i,1} = pos_ch1;
    peaks_raw = [peaks_raw; all_positions{i,1}];
    i
end

display(['You have ' num2str(N_peaks) ' peaks'])
cd(path_out)

avg_img = avg_img_green;

save 'data.mat' 'ch1' 'all_positions' 'N_peaks' 'avg_img' 'path_out' 'N_movie'


%% Select all peak coordinates for Red channel based on cross-correlation of channels
all_positions_red = all_positions;
corr_matrix = [];
shift = zeros(N_movie,2);
for i = 1:N_movie
    corr_matrix = normxcorr2(avg_img_green{i}, avg_img_red{i});
    [shift(i,2) shift(i,1)] = find(corr_matrix == max(corr_matrix(:)));
    shift(i,:) = shift(i,:) - length(avg_img_red{i});
    all_positions_red{i} = all_positions_red{i} + [shift(i,:) 0 0 0];
end

save 'data.mat' 'ch1' 'ch2' 'all_positions' 'all_positions_red' 'N_peaks' 'avg_img' 'avg_img_red' 'shift' 'path_out' 'N_movie'

%% export all Green peaks as .tif and .fits
dim_mov_out = 41; %41x41 pixels
dim_mov_out = (dim_mov_out-1)/2;
ROI_size = size(avg_img{1},1);

tic
mkdir('spots_tmp')
cd('spots_tmp')



% go through all movies and export from from each movie all spots
for i = 1:N_movie
    i
    frames_start = 1:ch1{i}.N_frame_per_fits:length(ch1{i}.frames);
    frames_stop = frames_start+ch1{i}.N_frame_per_fits-1;
    frames_stop(end) = length(ch1{i}.frames)+1;
    
    
    for j = 1:length(ch1{i}.fname)
        tmp_fullvideo = fitsread([pname filesep ch1{i}.fname{j}]);
        for k = 1:size(all_positions{i,1},1)
            ymin = all_positions{i,1}(k,1)-dim_mov_out;
            ymax = all_positions{i,1}(k,1)+dim_mov_out;
            xmin = all_positions{i,1}(k,2)-dim_mov_out;
            xmax = all_positions{i,1}(k,2)+dim_mov_out;
            if xmin < 1
                xmin = 1;
            end
            if xmax > ROI_size;
                xmax = ROI_size;
            end
            if ymin < 1
                ymin = 1;
            end
            if ymax > ROI_size;
                ymax = ROI_size;
            end
            
            mov_out = tmp_fullvideo(xmin:xmax, ymin:ymax, :);
            fitswrite(mov_out, ['temp_Mov_' num2str(i) 'Index' num2str(j) 'spot_' num2str(k) 'ch' num2str(1) '.fits'])
            clear mov_out
        end
        clear tmp_fullvideo
        
    end
    
    for l = 1:size(all_positions{i,1},1)
        tmp_video = zeros(41, 41, length(ch1{i}.frames), 'int16');
        [sizeX, sizeY, dim] = size(fitsread(['temp_Mov_' num2str(i) 'Index' num2str(j) 'spot_' num2str(l) 'ch' num2str(1) '.fits']));
        
        for j = 1:length(ch1{i}.fname)
            tmp_video(1:sizeX,1:sizeY,frames_start(j):frames_stop(j)) = fitsread(['temp_Mov_' num2str(i) 'Index' num2str(j) 'spot_' num2str(l) 'ch' num2str(1) '.fits']);
            filename = ['temp_Mov_' num2str(i) 'Index' num2str(j) 'spot_' num2str(l) 'ch' num2str(1) '.fits'];
            
            delete(filename)
        end
        tmp_video(:,:,1) = [];
        sum_mov_out = sum(tmp_video,3);
        sum_mov_out = sum_mov_out-min(sum_mov_out(:));
        max_value = max(max(sum_mov_out));
        factor = 65536/max_value;
        sum_mov_corr = uint16(sum_mov_out.*factor);
        
        fitswrite(tmp_video, ['movie_' num2str(i) 'spot_' num2str(l) 'ch' num2str(1) '.fits'])
        imwrite(sum_mov_corr, ['movie_' num2str(i) 'spot_' num2str(l) 'ch' num2str(1) '.tif'])
        l;
        clear tmp_video
    end
    movefile('movie*', [path_out '/spots'])
end
toc

%% export all Red peaks as .tif and .fits
% tic
% mkdir('spots_tmp')
% cd('spots_tmp')
%
% frames_start = 1:ch2{1}.N_frame_per_fits:length(ch2{1}.frames);
% frames_stop = frames_start+ch2{1}.N_frame_per_fits-1;
% frames_stop(end) = length(ch2{1}.frames)+1;
%
% % go through all movies and export from from each movie all spots
% for i = 1:N_movie
%     for j = 1:length(ch2{i}.fname)
%         tmp_fullvideo = fitsread([pname filesep ch2{i}.fname{j}]);
%         for k = 1:size(all_positions_red{i,1},1)
%             ymin = all_positions_red{i,1}(k,1)-dim_mov_out;
%             ymax = all_positions_red{i,1}(k,1)+dim_mov_out;
%             xmin = all_positions_red{i,1}(k,2)-dim_mov_out;
%             xmax = all_positions_red{i,1}(k,2)+dim_mov_out;
%             if xmin < 1
%                 xmin = 1;
%             end
%             if xmax > ROI_size;
%                 xmax = ROI_size;
%             end
%             if ymin < 1
%                 ymin = 1;
%             end
%             if ymax > ROI_size;
%                 ymax = ROI_size;
%             end
%
%             mov_out = tmp_fullvideo(xmin:xmax, ymin:ymax, :);
%             fitswrite(mov_out, ['temp_Mov_' num2str(i) 'Index' num2str(j) 'spot_' num2str(k) 'ch' num2str(2) '.fits'])
%             clear mov_out
%             j
%         end
%         clear tmp_fullvideo
%     end
%
%     for l = 1:size(all_positions{i,1},1)
%         tmp_video = zeros(41, 41, length(ch1{i}.frames), 'int16');
%         [sizeX, sizeY, dim] = size(fitsread(['temp_Mov_' num2str(i) 'Index' num2str(j) 'spot_' num2str(l) 'ch' num2str(2) '.fits']));
%
%         for j = 1:length(ch1{i}.fname)
%             tmp_video(1:sizeX,1:sizeY,frames_start(j):frames_stop(j)) = fitsread(['temp_Mov_' num2str(i) 'Index' num2str(j) 'spot_' num2str(l) 'ch' num2str(2) '.fits']);
%             filename = ['temp_Mov_' num2str(i) 'Index' num2str(j) 'spot_' num2str(l) 'ch' num2str(2) '.fits'];
%
%             delete(filename)
%         end
%         tmp_video(:,:,1) = [];
%         sum_mov_out = sum(tmp_video,3);
%         sum_mov_out = sum_mov_out-min(sum_mov_out(:));
%         max_value = max(max(sum_mov_out));
%         factor = 65536/max_value;
%         sum_mov_corr = uint16(sum_mov_out.*factor);
%
%         fitswrite(tmp_video, ['movie_' num2str(i) 'spot_' num2str(l) 'ch' num2str(2) '.fits'])
%         imwrite(sum_mov_corr, ['movie_' num2str(i) 'spot_' num2str(l) 'ch' num2str(2) '.tif'])
%         l
%         clear tmp_video
%     end
%     movefile('movie*', [path_out '/spots_red'])
% end
% toc

