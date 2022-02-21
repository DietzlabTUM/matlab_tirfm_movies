%% startup
clc, clear all, close all
path0 = cd;

%change to correct path
run('/Users/username/Documents/MOVIE_TIRFM/functions/startup.m')

%% load fitted data (as hdf5, eg from Picasso)

hdf2data;

load('data.mat')

K = spots;

save 'data_info.mat' path_out spots frames K

%% Choose one of the following 4 options to specify applied AC Field
% option 1: (switched between on and off)

% 'Wait before' : Time before AC field is switched on [s]
% 'Duration'    : Time between electric field is switched [ms] (i.e. 1sec/frequency/2)
% 'Cycles'      : Number of AC field cycles
% 'Wait between': Time for AC field off blocks of AC field on [s]
% 'Time per frame' [ms]; 'Voltage' [V]
% '#Blocks'     : Number of blocks with AC field on

input = {'Wait before:', 'Duration:', 'Cycles', 'Wait between:', 'Time per Frame:', 'Voltage', '#Blocks'};
input_default = {'4', '100', '75', '4', '4', '20', '3'};
tmp = inputdlg(input, 'Parameters', 1, input_default);

waitp = str2double(tmp(1)); 
dur = str2double(tmp(2)); 
iter = str2double(tmp(3)); 
waitb = str2double(tmp(4)); 
tpf = str2double(tmp(5)); 
vol = str2double(tmp(6));
mod = str2double(tmp(7));

eField = [waitp dur iter waitb tpf vol mod];

%% Choose one of the following 4 options to specify applied AC Field
% option 2: (frequency sweep)

% 'Wait before'     : Time before AC field is switched on [s]
% 'Duration(1-4)'   : Time between electric field is switched for each frequency [ms] (i.e. 1sec/frequency/2)
% '#Frames'         : Number frames for each frequency
% 'Time per frame' [ms]; 'Voltage' [V]

input = {'Wait before:', 'Duration1:', 'Duration2:', 'Duration3:', 'Duration4:' '#Frames:', 'Time per Frame:', 'Voltage'};
input_default = {'1', '5', '50', '100', '500', '2000', '4', '20'};
tmp = inputdlg(input, 'Parameters', 1, input_default);

wait = str2double(tmp(1)); 
dur1 = str2double(tmp(2)); 
dur2 = str2double(tmp(3)); 
dur3 = str2double(tmp(4)); 
dur4 = str2double(tmp(5)); 
numframes = str2double(tmp(6));
tpf = str2double(tmp(7));
vol = str2double(tmp(8));

eField = [wait dur1 dur2 dur3 dur4 numframes tpf vol];

%% Choose one of the following 4 options to specify applied AC Field
% option 3: (voltage sweep)

% 'Voltage (1-6)'   : Voltages for each step
% 'Duration'        : Time between electric field is switched [ms] (i.e. 1sec/frequency/2)
% '#Frames'         : Number frames for each voltage
% 'Wait between'    : Time for AC field off between voltage steps [s]
% 'Time per frame' [ms]; 'Voltage' [V]

input = {'Voltage1:', 'Voltage2:', 'Voltage3:', 'Voltage4:', 'Voltage5:',...
    'Voltage6:', 'Duration', '#Frames:', 'Wait between:', 'Time per Frame:'};
input_default = {'0', '5', '10', '20', '30', '40', '100', '2000', '0', '4'};
tmp = inputdlg(input, 'Parameters', 1, input_default);
 
vol1 = str2double(tmp(1)); 
vol2 = str2double(tmp(2)); 
vol3 = str2double(tmp(3)); 
vol4 = str2double(tmp(4));
vol5 = str2double(tmp(5)); 
vol6 = str2double(tmp(6));
dur = str2double(tmp(7));
numframes = str2double(tmp(8));
waitb = str2double(tmp(9)); 
tpf = str2double(tmp(10));
dum1 = 0;
dum2 = 0;
dum3 = 0;

eField = [vol1 vol2 vol3 vol4 vol5 vol6 dur numframes waitb tpf dum1 dum2 dum3];

%% Choose one of the following 4 options to specify applied AC Field
% option 4: (field axis rotation)

% 'Wait before'     : Time before AC field is switched on [s]
% 'Min Angle'       : Starting AC field axis angle [°]
% 'Max Angle'       : Maximum AC field axis angle [°]
% 'Degree Step'     : Stepsize between two AC field axis angles [°]
% 'Duration'        : Time between electric field is switched for each angle [ms] (i.e. 1sec/frequency/2)
% 'Cycles'          : Number of AC field cycles
% 'Time per frame' [ms]; 'Voltage' [V], 'Frequency' [Hz]

input = {'Wait before:', 'Min Angle:', 'Max Angle:', 'Degree Step:', 'Duration:' 'Cycles', 'Time per Frame:', 'Voltage:', 'Frequency'};
input_default = {'4', '0', '180', '15', '100', '20', '4', '20', '5'};
tmp = inputdlg(input, 'Parameters', 1, input_default);

wait = str2double(tmp(1)); 
minang = str2double(tmp(2)); 
maxang = str2double(tmp(3)); 
degstep = str2double(tmp(4)); 
dur = str2double(tmp(5)); 
cycles = str2double(tmp(6));
tpf = str2double(tmp(7));
vol = str2double(tmp(8));
freq = str2double(tmp(9));

eField = [wait minang maxang degstep dur cycles tpf vol freq];

%% Go through all fitted spots and select class of respective spot

% in case not all blocks of AC field cycles are along the same axis (chose
% between vertically ('UD') and horizontally ('RL')

if length(eField) == 7
    input = {['Sequence of AC field direction for ' num2str(eField(7))  ' blocks?']};
    input_default = {'UD - UD - UD'};
    seq = inputdlg(input, 'Parameters', 1, input_default);
    save seq.mat seq
end
    
mkdir('spots_analyzed')

% remember: counter-clockwise (-) and clockwise (+) switches

i = 1; % Set starting value here
K = spots;
go_on = 1;

while go_on 
    
    
    
    disp(['To go : ' num2str(K-i)]);
    key = 0;
   
    if exist('seq')
        [class, dwell, pick, Theta_D_shifted, Theta_cum_D, Theta_bet_D, rotation, net_rot, ind, turns, new_mean] = plot_spot_eField_nocb_gyro_col(i, spots, data, eField, seq);    
    else
        [class, dwell, pick, Theta_D_shifted, Theta_cum_D, Theta_bet_D, rotation, net_rot, ind, turns, new_mean] = plot_spot_eField_nocb_gyro_col(i, spots, data, eField);    
    end
    
    button_cw = uicontrol('Style', 'pushbutton', 'String', 'clockwise (+)', 'Position', [70 20 200 35], 'Callback', 'key = 1;'); 
    button_ccw = uicontrol('Style', 'pushbutton', 'String', 'counter-clockwise (-)', 'Position', [270 20 200 35], 'Callback', 'key = 2;'); 
    button_non = uicontrol('Style', 'pushbutton', 'String', 'no net rotation', 'Position', [470 20 200 35], 'Callback', 'key = 3;'); 
    button_change = uicontrol('Style', 'pushbutton', 'String', 'direction change', 'Position', [670 20 200 35], 'Callback', 'key = 5;');
    button_unclear = uicontrol ('Style', 'pushbutton', 'String', 'unclear', 'Position', [970 20 200 35], 'Callback', 'key = 4;');
    
    while waitforbuttonpress
    end


    pause(1.0)
    %disp(key)
    if isempty(key)
        go_on = 0;
    elseif key==0 % Repeat
        disp('Repeat')
        i = i-1;
    else
        disp([num2str(i) ' classified'])
        manual = 0;
        data{i}.type = key;
        data{i}.dwell = dwell;
        data{i}.clusters = pick;
        data{i}.rotation = rotation;
        data{i}.class = class;
        data{i}.theta = Theta_D_shifted;
        data{i}.theta_cum = Theta_cum_D;
        data{i}.theta_bet = Theta_bet_D;
        data{i}.nrot = net_rot;
        data{i}.maxpos = ind;
        data{i}.turns = turns;
        data{i}.center = new_mean;
        movefile(['spot' num2str(i) '_analysis.tiff'], [path_out '/spots_analyzed'])

    end
    

    if i==K
        go_on = 0;
        disp('end')
    else
        i=i+1;
    end

    close all;

end

save 'data_dwell.mat' data
save 'data_info.mat' path_out eField spots frames

%% histogram of net rotation

netrot = zeros(1,0);
for i = 1:spots
    if data{i}.type ~= 4
        netrot = [netrot; data{i}.turns];
    end
end

minrot = min(netrot);
minrot = floor(minrot);
maxrot = max(netrot);
maxrot = ceil(maxrot);
range = [minrot:1:maxrot];
fd = fitdist(netrot, 'Normal');
cutoff = round(1.645*fd.sigma);

figure()
hold all
h=hist(netrot, range);
hist(netrot, range);

kern = fitdist(netrot, 'Kernel');
ykern = pdf(kern, range);
ykern  =ykern*sum(h);
plot([minrot:1:maxrot], ykern, 'r')

title (['histogram of netto rotations with 90% threshold at +/-' num2str(cutoff) ' turns'])

savefig('netrot_hist.fig')
print(gcf, '-dtiff', 'netrot_hist.tiff');

save 'data_turns.mat' cutoff netrot maxrot minrot range ykern

%% sort analysed particles
cd('spots_analyzed')
mkdir('analyzed - cw')
mkdir('analyzed - ccw')
mkdir('analyzed - nnr')
mkdir('analyzed - discarded')
mkdir('analyzed - dirc')

for i=1:spots
    if data{i}.type == 1
        movefile(['spot' num2str(i) '_analysis.tiff'], [path_out '/spots_analyzed/analyzed - cw'])
    elseif data{i}.type == 2
        movefile(['spot' num2str(i) '_analysis.tiff'], [path_out '/spots_analyzed/analyzed - ccw'])
    elseif data{i}.type == 3
        movefile(['spot' num2str(i) '_analysis.tiff'], [path_out '/spots_analyzed/analyzed - nnr'])
    elseif data{i}.type == 4
        movefile(['spot' num2str(i) '_analysis.tiff'], [path_out '/spots_analyzed/analyzed - discarded'])
    elseif data{i}.type == 5
        movefile(['spot' num2str(i) '_analysis.tiff'], [path_out '/spots_analyzed/analyzed - dirc'])
    end
end

%% localize rotation direction in fit overview image (as tiff)

cd(path_out)

f = msgbox('open rendered tiff file with resolution x10');
uiwait(f);

[FileName, PathName] = uigetfile('*.tif', 'open rendered tiff file with resolution x10');
overview = fullfile(PathName,FileName);

imshow(overview); hold on
centers = zeros(spots,2);

for i = 1:spots
    if data{i}.type == 1  % clockwise(+)
        center = data{i}.center;
        center = center*10; %upscale
        plot(center(1), center(2), 'ro', 'MarkerSize', 20)
        centers(i,:) = center;

    elseif data{i}.type == 2  % counter-clockwise (-)
        center = data{i}.center;
        center = center*10; %upscale
        plot(center(1), center(2), 'bo', 'MarkerSize', 20)
        centers(i,:) = center;

    elseif data{i}.type == 3  % no net rotation (0)
        center = data{i}.center;
        center = center*10; %upscale
        plot(center(1), center(2), 'go', 'MarkerSize', 20)
        centers(i,:) = center;

    elseif data{i}.type == 5  % direction change (+/-)
        center = data{i}.center;
        center = center*10; %upscale
        plot(center(1), center(2), 'yo', 'MarkerSize', 20)
        centers(i,:) = center;
    end
end


    savefig('overview_analyzed.fig')
    print(gcf, '-dtiff', 'overview_analyzed.tiff');
    save 'overview.mat' centers overview PathName

close all

%% Calculate number of spots in each class

classes = zeros(4,2);

    for j = 1:spots
            if data{j}.type == 1 
                classes(1,2) = classes(1,2) + 1;
            elseif data{j}.type == 2 
                classes(2,2) = classes(2,2) + 1;
            elseif data{j}.type == 3 
                classes(3,2) = classes(3,2) + 1;
            elseif data{j}.type == 5
                classes(4,2) = classes(4,2) + 1;
            end
    end  
    
for m = 1:3
    classes(m,1) = m;
end 
    classes(4,1) = 5;
    
cd(path_out)
save 'classes.mat' 'classes'

%% Plot overview of classes 
    
if classes(4,2) == 0    
    cur_fig = figure;
    bar(classes(1:3,2))
    set(gca, 'xticklabel', {'clockwise (+)', 'counter-clockwise (-)', 'no net rotation'})
else
    cur_fig = figure;
bar(classes(1:4,2))
set(gca, 'xticklabel', {'clockwise (+)', 'counter-clockwise (-)', 'no net rotation', 'direction change'})
end
savefig('classes.fig')
print(gcf, '-dtiff', 'classes.tiff');

%% Plot overview of cumulative angular displacement for all rotors

theta_cum_mat = zeros(0,length(data{1}.theta_cum));
turns_tot = 0;
counter = 0;

for i = 1:spots
    if data{i}.type ~= 4
        theta_cum_mat = [theta_cum_mat; data{i}.theta_cum];
        turns_tot = turns_tot + data{i}.turns;
        counter = counter + 1;
    end
end

turns_avg = turns_tot/counter;
theta_cum_avg = mean(theta_cum_mat,1);

    

cur_fig2 = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
hold all
    for i=1:spots
        if data{i}.type ~= 4
            plot(data{i}.theta_cum,  'Color', '#787777')
        end
    end
    plot(theta_cum_avg, 'Color', 'k', 'LineWidth', 3)
    title(['theta cumulative with average net rotation of ' num2str(turns_avg) ' turns'])
    savefig('Theta_cum_all.fig')
    print(gcf, '-dtiff', 'Theta_cum_all.tiff');
hold off

cur_fig3 = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);
    plot(theta_cum_avg, 'Color', 'k', 'LineWidth', 3)
    title(['theta cumulative with average net rotation of ' num2str(turns_avg) ' turns'])
    savefig('Theta_cum_avg.fig')
    print(gcf, '-dtiff', 'Theta_cum_avg.tiff');
hold off

close all

save 'theta_cum.mat' turns_avg turns_tot counter theta_cum_avg theta_cum_mat

%% crop individual movies for each rotor
  
f = msgbox('Choose the folder with all .tiff files');
uiwait(f);

pname=uigetdir(pname, 'Choose the folder with all .tiff files.');

files = pickFirstTiffFiles(pname);
N_files = length(files);  

cd(pname)

frames = 0;
info = cell(N_files, 1);
for i = 1:N_files
    info{i} = imfinfo (files{i});
    frames = frames + length(info{i});
end

cd(path_out)
% export all spots as .tif and .fits
tic
mkdir('spots_tmp')
cd('spots_tmp')

frames_start1 = 1:length(info{1}):frames;
frames_stop1 = frames_start1+length(info{1})-1;
frames_stop1(end) = frames;

% go through all movies and export from each movie all spots
for i = 1:N_files % goes through movies part 1 till ...
    
    mov_out = cell(spots,1);
    
    for k = 1:length(info{i}) % goes through frames in movie part i, spot j
        
        tmp_frame = imread(info{i}(1).Filename, k);
        tmp_frame = double(tmp_frame);
        
        for j = 1:spots;     % goes through spot locations in movie part i 
        

        
            if isfile(['temp_part_' num2str(i) 'spot' num2str(j) '.fits'])
                disp('file already in directory')
            else
                ymin = round((centers(j,2))/10-20);
                ymax = round((centers(j,2))/10+20);
                xmin = round((centers(j,1))/10-20);
                xmax = round((centers(j,1))/10+20);
                
                if xmin < 1
                   xmin = 1;
                end
                if xmax > info{i}(1).Width
                   xmax = info{i}(1).Width;
                end
                if ymin < 1
                   ymin = 1;
                end
                if ymax > info{i}(1).Height
                   ymax = info{i}(1).Height;
                end

                mov_out{j, 1}(:, :, k) = tmp_frame(ymin:ymax, xmin:xmax);

            end  
            
        end
        
        k
        clear tmp_frame
        
    end
    
    for j = 1:spots
    
        fitswrite(mov_out{j, 1}, ['temp_part_' num2str(i) 'spot' num2str(j) '.fits'])
        j
    
    end
    
    clear mov_out
   
    i
   
end

for j = 1:spots % goes through spots
    
    tmp_video = zeros(41, 41, frames, 'int16');
    [sizeX, sizeY, dim] = size(fitsread(['temp_part_' num2str(1) 'spot' num2str(j) '.fits']));

    for i = 1:N_files   
        
        tmp_video(1:sizeX,1:sizeY,frames_start1(i):frames_stop1(i)) = fitsread(['temp_part_' num2str(i) 'spot' num2str(j) '.fits']);
        filename = ['temp_part_' num2str(i) 'spot' num2str(j) '.fits'];

        %delete(filename)
        
    end
    
    sum_mov_out = sum(tmp_video,3);
    sum_mov_out = sum_mov_out-min(sum_mov_out(:));
    max_value = max(max(sum_mov_out));
    factor = 65536/max_value;
    sum_mov_corr = uint16(sum_mov_out.*factor);

    fitswrite(tmp_video, ['spot_' num2str(j) '.fits'])
    imwrite(sum_mov_corr, ['spot_' num2str(j) '.tif'])
    j
    clear tmp_video
    
    movefile('spot*', [path_out '/spots'])
end

cd(path_out)
rmdir('spots_tmp')

toc

%% read spots and adjust fit positions

cd('spots')

fnames = dir('*.tif');
fnamesfits = dir('*.fits');


for i=1:spots
    data{i}.newfit = data{i}.fit-(centers(i,:)/10-21);
end

save 'data_new.mat' data
