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
end

%% calculate sum of all frames and stdev + polar coordinates

dev_img_small = zeros(31, 31, numfids);
sum_img_small = zeros(31, 31, numfids);
pol_img = zeros(200, 360, numfids);
vals = zeros(31, 31, numfids);

for K = 1:numfids
    vals = fitsread([fnamesfits(K).folder, '/corrected/' fnamesfits(K).name(1:end-5) '_CORR.fits']);
    sz = size(vals);
    dev_img_small(1:sz(1),1:sz(2),K) = std(vals,0,3);
    sum_img_small(1:sz(1),1:sz(2),K) = sum(vals,3)./max(max(sum(vals,3))); % calculate sum of stacks
    sum_img_small(1:sz(1),1:sz(2),K) = (sum_img_small(1:sz(1),1:sz(2),K)-(min(min(sum_img_small(1:sz(1),1:sz(2),K)))))./(max(max(sum_img_small(1:sz(1),1:sz(2),K)))-(min(min(sum_img_small(1:sz(1),1:sz(2),K))))); %normalize image
    
    %Convert cartesian to polar coordinates
    pol_img(:,:,K) = ImToPolar(sum_img_small(1:sz(1),1:sz(2),K), 0, 0.7, 200, 360);
    pol_img(:,:,K) = (pol_img(:,:,K)-(min(min(pol_img(:,:,K)))))./(max(max(pol_img(:,:,K)))-(min(min(pol_img(:,:,K)))));
    K
end

disp('ImToPol done')
save sum_dev_small.mat dev_img_small sum_img_small pol_img

%% Go through all fitted spots that are classified as switching and see circular histogram


% Pick clusters in counter-clockwise rotation!

i = 1; % Set starting value here

go_on = 1;

while go_on
    
    
    key = 0;
    
    disp(['To go : ' num2str(length(data)-i)]);
    
    
    [class, dwell, pick, Theta_D_shifted, rotation] = plot_spot_vwcm_nocb_gyro_circhist_col(i, N_peaks, data, fnamesfits, dev_img_small);
    
    
    button_switch3 = uicontrol('Style', 'pushbutton', 'String', 'Switching - 3', 'Position', [150 20 100 35], 'Callback', 'print(''-dtiff'', ''-r300'',  [num2str(fnamesfits(i).name(1:end-5)) ''_fig'']); key = 1');
    button_switch2 = uicontrol('Style', 'pushbutton', 'String', 'Switching - 2', 'Position', [300 20 100 35], 'Callback', 'key = 2;');
    button_stationary = uicontrol('Style', 'pushbutton', 'String', 'Stationary', 'Position', [450 20 100 35], 'Callback', 'key = 3;');
    button_rotating = uicontrol('Style', 'pushbutton', 'String', 'Rotating', 'Position', [600 20 100 35], 'Callback', 'key = 4;');
    button_unclear = uicontrol('Style', 'pushbutton', 'String', 'No full Circle', 'Position', [750 20 100 35], 'Callback', 'key = 5;');
    button_maybe = uicontrol('Style', 'pushbutton', 'String', 'Unclear', 'Position', [900 20 100 35], 'Callback', 'key = 6;');
      
    while waitforbuttonpress
    end
    
    pause(0.7)
    
    if isempty(key)
        go_on = 0;
    elseif key==0 % Repeat
        disp('Repeat')
        i = i-1;
    else
        disp([num2str(i) ' classified'])
        data{i}{1,1}.type = key;
        data{i}{1,1}.dwell = dwell;
        data{i}{1,1}.clusters = pick;
        data{i}{1,1}.rotation = rotation; %rotation can be +1, -1 or 0 according to where rotor is rotating; done per each frame
        data{i}{1,1}.class = class;
        data{i}{1,1}.Theta_D_shifted = Theta_D_shifted;
        %counter = counter-1;
    end
    
    if i==numfids
        go_on = 0;
        disp('end')
    else
        i=i+1;
    end
    
    close all;
    
    
end

save 'data_dwell.mat' data

%% Count classes

classes_AS = zeros(6,2);

for j = 1:length(data)
    if data{j}{1,1}.type == 1
        classes_AS(1,2) = classes_AS(1,2) + 1;
    elseif data{j}{1,1}.type == 2
        classes_AS(2,2) = classes_AS(2,2) + 1;
    elseif data{j}{1,1}.type == 3
        classes_AS(3,2) = classes_AS(3,2) + 1;
    elseif data{j}{1,1}.type == 4
        classes_AS(4,2) = classes_AS(4,2) + 1;
    elseif data{j}{1,1}.type == 5
        classes_AS(5,2) = classes_AS(5,2) + 1;
    elseif data{j}{1,1}.type == 6
        classes_AS(6,2) = classes_AS(6,2) + 1;
    end
end

for m = 1:6
    classes_AS(m,1) = m;
end

disp(classes_AS)

%% save name of files sorted
movies_switching3 = cell(1,K);
movies_switching2 = cell(1,K);
movies_rotating = cell(1,K);
movies_nofullcircle = cell(1,K);

for j = 1:K
    if data{j}{1,1}.type == 1
        movies_switching3(j) = ch1{j}.fname;
    elseif data{j}{1,1}.type == 2
        movies_switching2(j) = ch1{j}.fname;
    elseif data{j}{1,1}.type == 4
        movies_rotating(j) = ch1{j}.fname;
    elseif data{j}{1,1}.type == 5
        movies_nofullcircle(j) = ch1{j}.fname;
    end
end

movies_switching3 = movies_switching3';
movies_switching2 = movies_switching2';
movies_rotating = movies_rotating';
movies_nofullcircle = movies_nofullcircle';

disp('Done')

%% Save new data file

%cd(path_out),
data_corr=data;
save 'data_categorized.mat' 'data' 'data_corr' 'path_out'
save 'data_dwell.mat' 'data'
save 'movies_classified.mat' 'movies_switching3' 'movies_switching2' 'movies_rotating' 'movies_nofullcircle'

%%
types=[];
for i = 1:length(fnamesfits)
    data{i}{1,1}.name = getfield(fnamesfits(i),'name');
    types(i) = data{i}{1,1}.type;
end

types = types';

mkdir analysis
cd analysis
save 'data_names.mat' data


%% save names in a file
fileID = fopen('Names.txt','w');
% %fprintf(fileID,'%s\t %i\n', fnamesfits.name, types);
% fprintf(fileID,'%s\t %i\n',fnamesfits.name, types);

for i=1:length(fnamesfits)
    fprintf(fileID, '%s ', fnamesfits(i).name');
    fprintf(fileID,'\t');
    fprintf(fileID,'%i', types(i));
    fprintf(fileID,'\n');
end

fclose(fileID);

%% Get positions from VWCM - ALL PARTICLES

xy=zeros(length(data{i}{1,1}.vwcm.pos),length(fnamesfits));

index = 1;
for i=1:length(fnamesfits)
    xy(:,index) = [data{i}{1,1}.vwcm.pos(:,1)];
    xy(:,index+1) = [data{i}{1,1}.vwcm.pos(:,2)];
    index = index+2;
end

save 'xy_pos.mat' xy
writematrix(xy,'xy.csv','Delimiter','tab')

%% Get histogram theta & RMSD

K = length(fnamesfits);
rad1 = 10;
rad2 = 10;

index = 1;

for i=1:length(fnamesfits)
    if data{i}{1,1}.type == 1 || data{i}{1,1}.type == 2 || data{i}{1,1}.type == 4 || data{i}{1,1}.type == 5
        
        xydata = data{i}{1,1}.vwcm.pos;
        pick = data{i}{1,1}.clusters;
        new_mean = mean(pick);
        
        
        Theta_cum(1) = 0;
        
        normal = pick(1,:) - new_mean;
        normal = normal./norm(normal);
        [Theta_norm,Rho] = cart2pol(normal(1,1),normal(1,2));
        Theta_norm = wrapTo2Pi(Theta_norm);
        Theta_norm_shifted = Theta_norm + (60*pi/180);
        
        Theta_norm_D = Theta_norm*180/pi;
        Theta_norm_D_shifted = Theta_norm_shifted*180/pi;
        R = [cosd(Theta_norm_D) -sind(Theta_norm_D); sind(Theta_norm_D) cosd(Theta_norm_D)];
        R_shifted = [cosd(Theta_norm_D_shifted) -sind(Theta_norm_D_shifted); sind(Theta_norm_D_shifted) cosd(Theta_norm_D_shifted)];
        
        
        for k = 1:(size(data{i}{1,1}.vwcm.pos,1)-1) % Goes through all frames
            
            vector = data{i}{1,1}.vwcm.pos(k,:) - new_mean;
            vector = vector./norm(vector);
            vector = vector*R;
            vector_shifted = vector*R_shifted;
            
            vector2 = data{i}{1,1}.vwcm.pos(k+1,:) - new_mean;
            vector2 = vector2./norm(vector2);
            vector2 = vector2*R;
            vector2_shifted = vector2*R_shifted;
            
            [Theta(k),Rho] = cart2pol(vector(1,1),vector(1,2));
            [Theta_shifted(k),Rho] = cart2pol(vector_shifted(1,1),vector_shifted(1,2));
            
            [Theta_2(k),Rho] = cart2pol(vector2(1,1),vector2(1,2));
            
            
            Theta(k) = wrapTo2Pi(Theta(k));
            Theta_shifted(k) = wrapTo2Pi(Theta_shifted(k));
            Theta_2(k) = wrapTo2Pi(Theta_2(k));
            
            %Theta inbetween
            Theta_bet(k) = Theta_2(k) - Theta(k);
            
            %ensure shortest distence has been chosen (i.e. forbid theta_bet>180)
            if abs(Theta_bet(k)) < pi
                Theta_bet(k) = Theta_bet(k);
            else
                if Theta_bet(k) < 0
                    Theta_bet(k) = Theta_bet(k) + 2*pi;
                else
                    Theta_bet(k) = Theta_bet(k) - 2*pi;
                end
            end
            
            Theta_cum(k+1) = Theta_bet(k) + Theta_cum(k);
            
        end
        
        
        Theta_D = Theta*180/pi;
        Theta_D_shifted = Theta_shifted*180/pi;
        Theta_2_D = Theta_2*180/pi;
        Theta_bet_D = Theta_bet*180/pi;
        Theta_cum_D = Theta_cum*180/pi;
        Theta_cum_D=Theta_cum_D';
        
        %root mean square deviation
        Theta_D0 = Theta_cum_D(1,1);
        RMSD_sq = zeros(1,length(Theta_cum_D));
        RMSD_sq(1,1)=(Theta_cum_D(1,1)-Theta_D0)^2;
        
        for l=2:length(Theta_cum_D)
            RMSD_sq(l) = ((Theta_cum_D(l)-Theta_D0)^2+RMSD_sq(l-1));
        end
        
        RMSD_sq_temp = zeros(1,length(RMSD_sq));
        for l=1:length(RMSD_sq_temp)
            RMSD_sq_temp(l) = RMSD_sq(l)/l;
        end
        
        RMSD_1 = sqrt(RMSD_sq_temp);
      
        
        RMSD_particle(:,index) = [RMSD_1];
       
        
        
        %angular velocity
        w = zeros(1,length(Theta_bet));
        w(1,1) = 0;
        
       
        
        for l=2:length(Theta_bet)
            w(l) = Theta_bet(l)/1;
        end
        
        w_particle(:,index) = [w];
        %h_w = histogram(w_particle,'Normalization','probability');
        %w_hist(:,index) = [h_w.Values];
        
        space_particle(index) = sum(abs(Theta_bet_D));
        
        %angular hist
        h = histogram(Theta_D, 360);
        [max_h, ind] = max(h.Values);
        h_bins = h.BinEdges;
        shift = ind - 60;
        Theta_D_60 = Theta_D - shift;
        Theta_D_60_rad = Theta_D_60 *pi/180;
        Theta_D_60_rad = wrapTo2Pi(Theta_D_60_rad);
        Theta_D_60_D = Theta_D_60_rad * 180/pi;
        
        h_bin1 = histogram(Theta_D_60_D, 360);
        
        Theta_cum_particle(:,index) = [Theta_cum_D];
        Theta_particle(:,index) = [Theta_D_shifted];
        Theta_hist(:,index) = [h_bin1.Values];
        
        
        
        sum_hist_shift = sum(Theta_hist,2);
        avg_hist_shift = sum_hist_shift/(index);
        
        Theta_hist_norm = Theta_hist./sum(Theta_hist);
        avg_hist_shift_norm = sum(Theta_hist_norm,2)/(index);
        
        h_bin3 = histogram(Theta_D_60_D, 72);
        Theta_hist_3(:,index) = [h_bin3.Values];
        Theta_hist_norm_3 = Theta_hist_3./sum(Theta_hist_3);
        
        avg_hist_shift_norm_3 = sum(Theta_hist_norm_3,2)/(index);
        
        
        
        
        index = index+1;
        
    end
end

RMSD_avg = sum(RMSD_particle,2)/(index-1);
w_particle_all = reshape(w_particle,[],1);
space_particle = space_particle';

save 'theta.mat' Theta_particle Theta_hist sum_hist_shift avg_hist_shift Theta_hist_norm Theta_hist_norm_3 avg_hist_shift_norm avg_hist_shift_norm_3
save 'RMSD-omega.mat' Theta_cum_particle RMSD_particle RMSD_avg w_particle w_particle_all space_particle%w_hist
writematrix(Theta_particle,'Theta_particle.csv','Delimiter','tab')
writematrix(Theta_hist,'Theta_hist.csv','Delimiter','tab')

disp('Done')

%% Get rotations
rotations=[];
index2 = 1;
for i=1:length(data)
    if data{i}{1,1}.type == 1 || data{i}{1,1}.type == 2 || data{i}{1,1}.type == 4 || data{i}{1,1}.type == 5
        
        rotations(index2,:) = [data{i}{1,1}.rotation];
        index2 = index2+1;
    end
end

rotations = rotations';
save 'rotations.mat' rotations
writematrix(rotations,'rotations.csv','Delimiter','tab')

%% superimpose vwcm location per frame on movies of perfect particles
mkdir('VWCM_superimposed')
cd('VWCM_superimposed')


for i=1:length(data)
    
    movie_out = VideoWriter([fnamesfits(i).name(1:end-5) '_vwcm_superimp.avi'], 'Uncompressed AVI');
    open(movie_out)
    movie_in = fitsread([fnamesfits(i).folder, '/corrected/' fnamesfits(i).name(1:end-5) '_CORR.fits']);
    movie_in = movie_in(:,:,2:end);
    crosscolor = {'blue' 'red' 'green'};
    for j = 1:length(movie_in)
        frame_in = (movie_in(:,:,j)-min(min(movie_in(:,:,j))))./(max(max(movie_in(:,:,j)))-min(min(movie_in(:,:,j))));
        frame_out = insertMarker(frame_in, data{i}{1,1}.vwcm.pos(j,:), '+', 'color', crosscolor{data{i}{1,1}.class(j)}, 'size', 1);
        writeVideo(movie_out, frame_out);
    end
    close(movie_out)
    i
end


%% Distance from centre

index = 1; 

%for i=1:length(fnamesfits)
for i=1:fnamesfits
    if data{i}{1,1}.type == 1 || data{i}{1,1}.type == 2 %|| data{i}{1,1}.type == 4 || data{i}{1,1}.type == 5
        
        xydata = data{i}{1,1}.vwcm.pos;
        pick = data{i}{1,1}.clusters;
        new_mean = mean(pick);
        
        for j = [1:length(data{i}{1,1}.vwcm.pos(:,1))]
            
            if (data{i}{1,1}.vwcm.pos(j,1) == 0)
                d(j) = 0;
            else
                d(j) = sqrt((data{i}{1,1}.vwcm.pos(j,1)-new_mean(1)).^2+(data{i}{1,1}.vwcm.pos(j,2)-new_mean(2)).^2);
            end
        end
        %plot(d)
        mean_d(i) = mean(d(d~=0));
        %ylim ([0,10])
        %title(['Normalized distance from mean value, mean is ' num2str(mean_d)]);
        d_particle(:,index) = [d];
        index = index+1;
    end
    
end

save 'Distance_particles_2+3.mat' d_particle