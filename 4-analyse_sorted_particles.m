% This and other related scripts were written by Dr. Matthias Schickinger, Dr. Philip Ketterer, Dr. Jonas Funke,
% Anna-Katharina Pumm, Dr. Wouter Engelen and Eva Bertosin

%% load spots

%change to correct path

load('/Users/username/20xx-xx-xx_analysis/data.mat')
load('/Users/username/selected_for_analysis/dev_img.mat')
cd '/Users/username/selected_for_analysis'

%% LOAD STACK OF MOVIES Green Channel

channel = cell(1,1);
channel{1} = 'movie';

pname=pwd;
files_ch1 = pickFirstFitsFiles(pname, channel{1});

%% LOAD STACK OF MOVIES Green Channel
channel{2} = 'green';
for i=1:K
    files_ch2{i} = [files_ch1{i}(1:end-6) '2.fits'];
end

%% set parameters and generate movie classes
K=length(files_ch1);
ch1 = cell(K,1);

first = 2;
last = -1;
time_per_frame = 50;
r_find = 8;
r_integrate = 8;
N_frames = 100;
sequence_ch1 = 1;

for i=1:K
    ch1{i} = movie(pname, files_ch1{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
    if length(channel) == 2
        ch2{i} = movie([path_out '/spots_green'], files_ch2{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
    end
    i
end

disp('Done')

%% find center of rotor -- click in center of rotor

peaks_raw = zeros(0,5);
all_positions = cell(K, 1);

for i=1:K
    [pos_ch1] = ch1{i}.get_peaks_center(r_find, ch1{i}.mov_length-1);
    all_positions{i,1} = pos_ch1;
    peaks_raw = [peaks_raw; all_positions{i,1}];
    disp(i)
end

N_peaks = size(peaks_raw,1);
display(['You have ' num2str(N_peaks) ' peaks'])

%% adjust fits size -- adjust dimensions of fits to fit around rotor
dim = 41; %x/y dimensions of corrected movies
dim = (dim-1)/2;
mkdir ('corrected')

%properties of mask to remove neighboring particles in fits
rad = 15;
x = -dim:dim;
y = -dim:dim;
[xx yy] = meshgrid(x,y);
mask = zeros(size(xx));
mask((xx.^2+yy.^2)<rad^2)=1;   % radius 15, center at the origin

tic
for i = 1:K
    
    Nfr = ch1{i}.mov_length;
    mov_corr_tmp = zeros(100, 100, Nfr);
    tmp = fitsread(ch1{i}.info{1}.Filename);
    
    min_val = min(tmp(:));
    %add values of each frame to mov_corr_tmp, padded with
    %min_val at all sides with thickness dim.
    mov_corr_tmp(1:ch1{i}.sizeX+(2*dim), 1:ch1{i}.sizeX+(2*dim), :) = min_val;
    mov_corr_tmp(dim:ch1{i}.sizeX+dim-1, dim:ch1{i}.sizeY+dim-1, :) = tmp;
    %select only the part with thickness dim around center
    mov_corr = mov_corr_tmp((peaks_raw(i,2)-1):(peaks_raw(i,2)+2*dim-1), (peaks_raw(i,1)-1):(peaks_raw(i,1)+2*dim-1),:);
    
    %normalize each frame in mov_corr and multiply with circular mask
    for j = 1:length(mov_corr)
        mov_corr(:,:,j) = mov_corr(:,:,j) - min(setdiff(mov_corr(:,:,j),min(mov_corr(:,:,j)))); %(min(min(mov_corr(:,:,j)));
        %         mov_corr(:,:,j) = mov_corr(:,:,j) .* mask; %optionally apply circular mask
        mov_corr(:,:,j) = mov_corr(:,:,j) ./ max(max(mov_corr(:,:,j)));
    end
    
    
    fitswrite(int16(round(mov_corr*65535)), [ 'corrected' filesep ch1{i}.fname{1}(1:end-5) '_CORR.fits'])
    disp(i)
end

%set all peak centers to middle of corrected video (dim+1)
peaks_raw(:,1:2) = dim+1;
display ('Done')
toc
%% load corrected fits and generate movie classes of those
clear ch1
cd ('corrected')

channel = cell(1,1);
channel{1} = 'movie';

pname=pwd;
files_ch1 = pickFirstFitsFiles(pname, channel{1});

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
mkdir(path_out)

ch1 = cell(K,1);

sequence_ch1 = 1;

for i=1:K
    disp(i)
    ch1{i} = movie(pname, files_ch1{i}, first, last, sequence_ch1); % pname, fname, first, last, sequence
end

cd(path_out)

save 'data_spec.mat' 'ch1' 'all_positions' 'N_peaks' 'avg_img' 'path_out' 'N_movie' 'peaks_raw'

%% compute average images
avg_img_sub = cell(K, 2);

for i=1:K
    i
    avg_img_sub{i, 1} = ch1{i}.average_image(N_frames);
    avg_img_sub{i, 2} = ch1{i}.average_image_lastLong(N_frames); % for fitting threshold assignment, drift correction
end

%% Fit psf to spots
s_x = 2.5;
s_y = 2.5;
w_fit = 30;

ch1_fit_raw = zeros(N_peaks, 7);
ch1_fit_err_raw = zeros(N_peaks, 7);

h = waitbar(0,'Fitting spots.. please wait');

for i=1:N_peaks
    
    x1 = round(peaks_raw(i,1));
    y1 = round(peaks_raw(i,2));
    [c, c_err, ci, area] = fit_gauss2d_mainaxis_bg(x1, y1, s_x, w_fit, avg_img_sub{i,1});
    ch1_fit_raw(i,:) = c;
    ch1_fit_err_raw(i,:) = c_err;
    
    
    waitbar( i/N_peaks , h, ['Fitting spot... ' num2str(i) ' of ' num2str(N_peaks) ' done']) % update waitbar
end

close(h)

%% DON't SORT OUT

%remove not-accepted spots
ch1_fit = ch1_fit_raw(:, :);
ch1_fit_err = ch1_fit_err_raw(:, :);

peaks = peaks_raw(:,:);
% peaks = [ch1_fit(:,1:2) peaks(:,5)];

%% Get intensity traces 'itraces' plus median filtered itraces
display('Getting intensity traces... please wait')
tic
merged_itraces = cell(K,5);
merged_itraces_ch2 = cell(K,1);
iEndval = cell(K,2);
iEndval_sorted = cell(K,2);
avg_iEndval = zeros(K,2);
avg_iEndval_tenth = zeros(K,2);
for i=1:K
    %get fluorescence intensity traces from position
    ch1_itraces_full = ch1{i}.int_spots_in_frames(1:length(ch1{i}.frames), peaks(i,1:2), r_integrate);
    %     ch2_itraces_full = ch2{i}.int_spots_in_frames(1:length(ch2{i}.frames), peaks(i,1:2), r_integrate);
    display(['Tracing ' channel{1} ' channel in movie #' num2str(i) ' done'])
    toc %
    knumber = cell(size(ch1_itraces_full));
    knumber(:) = {i};
    merged_itraces{i,1} = ch1_itraces_full;
    %     merged_itraces_ch2{i} = ch2_itraces_full;
    merged_itraces{i,3} = knumber;
    
    %add median filtered itraces and average intensity values (over 100 frames) at the end
    tmp = cell(length(ch1_itraces_full),2);
    iEndval{i,1} = zeros(size(ch1_itraces_full,1),1);
    for j=1:length(ch1_itraces_full)
        tmp{j,1} = medfilt1(ch1_itraces_full{j}(:,4),20);
        iEndval{i,1}(j) = mean(ch1_itraces_full{j}(end-100:end,4));
    end
    merged_itraces{i,4} = tmp;
    for ch = 1
        avg_iEndval(i,ch) = mean(iEndval{i,ch});
        [tmp_val, tmp_spot] = sort(iEndval{i,ch});
        iEndval_sorted{i,ch} = [tmp_spot, tmp_val];
        avg_iEndval_tenth(i,ch) = mean(iEndval_sorted{i,ch}(1:ceil(length(iEndval_sorted{i,ch})/10),2));
    end
end
display('itraces complete')
toc %


%%
%pos_in_frame: cell of arrays that for each frame gives starting fit
%coordinates for all spots in respective channel. If both parameters
%return zero, spot is not fitted in that frame.

pos_in_frame = cell(K,2);
for m = 1:K
    % channel 1
    pos_in_frame{m,1} = cell(size(ch1{m}.frames,2),1);
    for j = 1:size(pos_in_frame{m,1},1)
        pos_in_frame{m,1}{j} = zeros(size(merged_itraces{m,1},1),2);
        for s=1:size(pos_in_frame{m,1}{j},1)
            pos_in_frame{m,1}{j}(s,1:2) = (merged_itraces{m,4}{s,1}(j)>=0)*(merged_itraces{m,1}{s}(j,2:3)); % x_0, y_0 remain zero if intensity is below threshold
        end
    end
    m
end
%% VWCM part

vwcm_output = cell(size(ch1));
for m=1:size(ch1,1)
    vwcm_output{m} = cell(size(pos_in_frame{m,1}{1},1),2);
    for j = 1:size(vwcm_output{m},1)
        vwcm_output{m}{j,1} = zeros(length(ch1{m}.frames),8);
        vwcm_output{m}{j,2} = zeros(length(ch1{m}.frames),8);
    end
end

parfor m=1:size(ch1,1)
    % channel 1
    tmp_vwcm_output = cell(length(ch1{m}.frames),1);
    for n = 1:size(tmp_vwcm_output,1)
        tmp_vwcm_output{n} = zeros(size(pos_in_frame{m,1}{n},1),8);
    end
    for p = 1:size(tmp_vwcm_output,1)
        [pos, delta, N, pos_max, v_max, stDev] = ch1{m}.vwcm_in_frame(ch1{m}.frames(p), pos_in_frame{m,1}{p}, 5, 0.01, 1000);
        tmp_vwcm_output{p} = [pos delta N pos_max v_max stDev];
    end
    %write VWCM loop data to vwcm_output cell
    for j=1:size(pos_in_frame{m,1}{1},1)
        for i=1:size(tmp_vwcm_output,1)
            vwcm_output{m}{j,1}(i,:) = tmp_vwcm_output{i}(j,:);
        end
    end
    
    display(['VWCM estimation in movie #' num2str(m) ' done'])
end
cd(path_out)
save -v7.3 'vwcmoutputs.mat' 'vwcm_output'

%% Save all relevant data
data = cell(K,1);
cd(path_out)
for m=1:K %loop through movies
    data{m} = cell(size(merged_itraces{m,1},1),2);
    for s=1:size(data{m},1)
        for ch = 1
            data{m}{s,ch}.pos0 = merged_itraces{m,ch}{s}(:,2:3);
            data{m}{s,ch}.itrace = merged_itraces{m,ch}{s}(:,4);
            data{m}{s,ch}.med_itrace = merged_itraces{m,4}{s,ch}(:);
        end
    end
end
% file that the position data from gF and vwcm estimators will be added to
% save -v7.3 'data_spot_pairs.mat' 'data' 'path_out'

% data needed for processing (batch jobs and later)
save -v7.3 'data_proc.mat' 'pos_in_frame' 'time_per_frame'


% movie objects
save -v7.3 'movie_objects.mat' 'ch1'

% stuff that might be useful for plotting figures
save 'data_plot.mat' 'channel'

% for archiving purposes
save -v7.3 'data_archive.mat' 'avg_img' 'N_frames' 'r_find' 'r_integrate' 'peaks_raw' 'peaks'

%% Output structure part
for m = 1:size(data,1)
    for s = 1:size(data{m},1)
        for c = 1
            % get data from vwcm estimate
            data{m}{s,c}.vwcm.pos = vwcm_output{m}{s,c}(:,1:2);
            data{m}{s,c}.vwcm.delta = vwcm_output{m}{s,c}(:,3);
            data{m}{s,c}.vwcm.N = vwcm_output{m}{s,c}(:,4);
            data{m}{s,c}.vwcm.pos_max = vwcm_output{m}{s,c}(:,5:6);
            data{m}{s,c}.vwcm.v_max = vwcm_output{m}{s,c}(:,7);
            data{m}{s,c}.vwcm.stDev = vwcm_output{m}{s,c}(:,8);
        end
    end
    display(['Structures in movie #' num2str(m) ' done'])
end

% Save data
cd(path_out)
save -v7.3 'data_spot_pairs.mat' 'data' 'path_out' 'merged_itraces_ch2'