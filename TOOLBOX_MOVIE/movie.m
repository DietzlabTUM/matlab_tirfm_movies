classdef movie < handle
    %TRACE_SELECTION Summary of this class goes here
    %   Detailed explanation goes here
     
    properties
        frames; %all images to be read
        N_read = 300; % number of images to read in one portion.
        counter = 1; % internal counter for reading the movie
         
        sequence; % sequence to be read, e.g. 101010
        first; % first image to be read
        last; % last image to be read
         
        pname; %pathname of file location
        fname; %filename of file
         
        sizeX; % number of pixel in X-dimension
        sizeY; % number of pixel in Y-dimension
        mov_length; % number of frames in thw whole movue
         
        info; %fits info
        h_min; %minimal heigth for peak fidning
         
        input; % 0=fits, 1=tiff-stack
        fnames; % cell with all filenames, only for tiff-stack
         
        N_frame_per_fits; % stores the number of frames in one fits-file
         
        drift; % stores displacement in x and y over whole movie through drift.
    end
     
    methods
        %constructor
        function obj = movie(pname, fname, first, last, sequence) % fname is the filename of the first fits-file
            obj.sequence = sequence;
            obj.pname = pname;
            obj.fname = cell(1,1);
            obj.fname{1} = fname;
             
            if strcmp(fname(end-2:end), 'tif')
                obj.input = 1; % read tiff data
            else
                obj.input = 0; % read fits data
            end
             
            if obj.input == 1 % tiff-stack
                tmp = dir([pname filesep '*.tif']);
                obj.fnames = {tmp.name};
                obj.info = 'Tiff-Stack info';
                obj.sizeX = size(imread([pname filesep obj.fnames{1}]),2);
                obj.sizeY = size(imread([pname filesep obj.fnames{1}]),1);
                obj.mov_length = length(obj.fnames);
            else % fits
                obj.info = cell(1,1);
                obj.info{1} = fitsinfo([obj.pname filesep obj.fname{1}]);
                 
                obj.N_frame_per_fits = 4095; %obj.info{1}.PrimaryData.Size(3);
                 
                obj.sizeX = obj.info{1}.PrimaryData.Size(1); 
                obj.sizeY = obj.info{1}.PrimaryData.Size(2);
                 
                tmp = dir([pname filesep fname(1:end-5) '_X*.fits']); %returns additional  change * to wildcard for 1-2 character/integers
                [~,idx] = sort([tmp.datenum]);
                tmp = tmp(idx);
                for i=1:length(tmp)
                    obj.fname = [obj.fname tmp(i).name];
                    obj.info = [obj.info fitsinfo([obj.pname filesep tmp(i).name])];
                end
                 
                f_tot = 0; % calculate total number of frames
                for i=1:length(obj.fname)
                    f_tot = f_tot + obj.info{i}.PrimaryData.Size(3);           
                end
                obj.mov_length = f_tot;
            end
             
             
            obj.first = first;
            if last == -1
                obj.last = obj.mov_length;
            else
                obj.last = last;
            end
            obj.frames = obj.getFrames(obj.sequence, obj.first, obj.last);
             
            % set drift to zero (= default value)
            obj.drift = zeros(length(obj.frames),2);
             
        end
         
        %generate a list of images to be read from the movie
         function frames = getFrames(obj, sequence, first, last)
            frames = [];
            tmp = first:last;
            for i=1:length(tmp)
               if(  sequence(  mod(i-1,size(sequence,2))+1 )  )
                  frames = [frames tmp(i) ];       
               end    
            end
         end
          
         %reads one frame
         function [img] = readFrame(obj,framenumber)
              
              if obj.input == 1 % tif
                    img = double(imread([obj.pname filesep obj.fnames{framenumber}]));
              else
                   i_fits = ceil(framenumber/obj.N_frame_per_fits);    % index of fits file
                   framenumber_effektive = mod(framenumber-1, obj.N_frame_per_fits) +1;
                   img = fitsread([obj.pname filesep obj.fname{i_fits}],  'Info', obj.info{i_fits}, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [framenumber_effektive framenumber_effektive] });  % fits
              end
         end
          
         % initialize the counter
         function  initRead(obj)
              obj.counter = 1;
         end
         
         %generate a list of images to be read from long movies
         function [img] = readLongFrame(obj,framenumber)
              
              if obj.input == 1 % tif
                    img = double(imread([obj.pname filesep obj.fnames{framenumber}]));
              else
                   framenumber_effektive = framenumber;
                   img = fitsread([obj.pname filesep obj.fname{1}],  'Info', obj.info{1}, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [framenumber_effektive framenumber_effektive] });  % fits
              end
         end
         
         %reads the next N_read frames
         function [tmp, frames_out, go_on] = readNext(obj)
            read_stop = obj.counter+obj.N_read-1; % last frame to read
            go_on = 1;
            if read_stop >= length(obj.frames)
                read_stop = length(obj.frames);
                go_on = 0;
            end
             
             
            tmp = zeros(obj.sizeX, obj.sizeY, read_stop-obj.counter+1);
            frames_out = obj.frames(obj.counter:read_stop);
             
            display(['Reading frame ' num2str(frames_out(1)) ' to ' num2str(frames_out(end))])
            for i=1:length(frames_out)
                tmp(:,:,i) = obj.readFrame(frames_out(i));
            end
             
            obj.counter = obj.counter + length(frames_out); 
         end

         %reads the next N_read frames for long movies
         function [tmp, frames_out, go_on] = readLongNext(obj)
            read_stop = obj.counter+obj.N_read-1; % last frame to read
            go_on = 1;
            if read_stop >= length(obj.frames)
                read_stop = length(obj.frames);
                go_on = 0;
            end
             
             
            tmp = zeros(obj.sizeX, obj.sizeY, read_stop-obj.counter+1);
            frames_out = obj.frames(obj.counter:read_stop);
             
            display(['Reading frame ' num2str(frames_out(1)) ' to ' num2str(frames_out(end))])
            for i=1:length(frames_out)
                tmp(:,:,i) = obj.readLongFrame(frames_out(i));
            end
             
            obj.counter = obj.counter + length(frames_out); 
         end
          
          
          
         % trace the movie 
         function [ traces, itraces, avg_frame ] = trace_movie(obj, h_min, r_find, r_integrate, min_length )
            traces = cell(0,1);
            itraces = cell(0,1);
            avg_frame = zeros(obj.sizeX, obj.sizeY);
             
            go_on = 1;
            obj.initRead;
            N = 0;
                      
            %display('Tracing movie... please wait')
            while go_on
                [movie, frames, go_on]  = obj.readNext;
                [traces, itraces] = append_traces(movie, traces, itraces, frames, h_min, r_find, r_integrate, min_length);
             
                avg_frame = avg_frame + sum(movie,3);
                N = N + length(frames);
            end
            avg_frame = avg_frame ./ N; 
            %display('Done tracing movie.')
         end
          
         % integrate specific regions in movie
         function [itraces ] = traces_movie_position(obj, positions, r_integrate )
            itraces = cell(0,1);
            go_on = 1;
            obj.initRead;
            while go_on
                [movie, frames, go_on]  = obj.readNext;
                itraces = append_traces_to_position(movie, itraces, frames, positions, r_integrate);
            end
         end
          
         % integrate intensities over specific areas in one frame (suitable for parfor)
         % includes coarse background subtraction 
         function [ints] = int_spots_in_frame(obj, frame_idx, spots_pos, d_int)
             %display(['Tracing frame #' num2str(frame_num)])
             ints = zeros(1,size(spots_pos,1));
             img = obj.readLongFrame(obj.frames(frame_idx));
             spots_pos = round(spots_pos + ones(size(spots_pos))*diag(obj.drift(frame_idx,:)));
             for i=1:size(spots_pos,1)
                 sub_img = img(max(1,spots_pos(i,2)-d_int):min(obj.sizeY,spots_pos(i,2)+d_int),...
                     max(1,spots_pos(i,1)-d_int):min(obj.sizeX,spots_pos(i,1)+d_int));
                 %sub_img = sub_img - min(sub_img(:));
                 ints(i) = sum(sub_img(:));
             end
         end
          
         % integrate intensities over specific areas in many frames;
         % return intensity trace
         function [itraces] = int_spots_in_frames(obj, frame_nums, spots_pos, d_int)
             itraces = cell(size(spots_pos,1),1);
             tmp = zeros(length(frame_nums),size(spots_pos,1));
             for j = 1:length(frame_nums)
                 tmp(j,:) = obj.int_spots_in_frame(frame_nums(j), spots_pos, d_int);
             end
             for s = 1:size(itraces,1)
                 itraces{s} = [reshape(frame_nums,length(frame_nums),1) zeros(length(frame_nums),3)];
                 itraces{s}(:,2:3) = ones(length(frame_nums),2)*diag(spots_pos(s,:)) + obj.drift;
                 itraces{s}(:,4) = tmp(:,s);
             end
         end
              
         % determine peak-finding threshholds
         function [h_min, p_out] = get_h_min(obj, r_find, N_img)
             
             if exist('N_img', 'var') % use average frame
                 img = obj.average_image(N_img); % average first N_img images
             else % use first frame if N_img is not specified
                if obj.input == 1 % tiff
                    img = double(imread([obj.pname filesep obj.fnames{obj.frames(1)}]));
                else % fits
                    img = fitsread([obj.pname filesep obj.fname{1}],  'Info', obj.info{1}, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [obj.frames(1) obj.frames(1)] }); % read first frame                
                end
             end
            p = find_peaks2d(img, r_find, 0, 0); % finding all possible peaks p has x, y, height, height-bg, I, I-I_bg
             
             
            close all
            figure('units','normalized','outerposition',[0 0 1 1])
            img_mean = mean(img(:));
            img_std = std(img(:));            
             
            p_7std = p(find(p(:,4)>=7*img_std), :); % this is just an estimate, #of peaks found may vary since peak_find algorithm return height as int not double
            p_5std = p(find(p(:,4)>=5*img_std), :);
            p_3std = p(find(p(:,4)>=3*img_std), :);
 
            %plot
            subplot(1, 2, 1)
            imagesc(img), colorbar, axis image, colormap gray, hold on
            if size(p_3std,1)>0
               h(1) = plot(p_3std(:,1)+1, p_3std(:,2)+1, 'ro');
            end
            if size(p_5std,1)>0
               h(2) = plot(p_5std(:,1)+1, p_5std(:,2)+1, 'go');
            end
            if size(p_7std,1)>0
               h(3) = plot(p_7std(:,1)+1, p_7std(:,2)+1, 'bo');
            end
            legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
                   
             
            subplot(1, 2, 2)
            xhist = min(p(:,4)):5:max(p(:,4));
             n = hist(p(:,4), xhist);
            semilogy(xhist, sum(n)-cumsum(n)), hold on
            h(1) = vline(3*img_std, 'r');
            h(2) = vline(5*img_std, 'g');
            h(3) = vline(7*img_std, 'b');
            legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
            set(gca, 'XLim', [0 xhist(end)])
            xlabel('Minimal height'), ylabel('# of peaks found')
            axis square
              
              
            % promp
             options.WindowStyle='normal';
             prompt={'Enter min heigth (default=5*sigma):'};
             def={num2str(round(5*img_std))};
             threshold = inputdlg(prompt, strcat('Enter threshold:'), 1, def, options);
             %threshold= round(5*img_std);
            h_min = str2double(threshold(1));
           close all
             
            obj.h_min = h_min;
                          
            p_out = p(find(p(:,4)>=h_min),:);
         end
          
         % determine peak-finding at concret pos
         function [p_out] = get_peaks_spec(obj, r_find, N_img)
             
             if exist('N_img', 'var') % use average frame
                 img = obj.average_image(N_img); % average first N_img images
             else % use first frame if N_img is not specified
                if obj.input == 1 % tiff
                    img = double(imread([obj.pname filesep obj.fnames{obj.frames(1)}]));
                else % fits
                    img = fitsread([obj.pname filesep obj.fname{1}],  'Info', obj.info{1}, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [obj.frames(1) obj.frames(1)] }); % read first frame                
                end
             end
            p = find_peaks2d(img, r_find, 0, 0); % finding all possible peaks p has x, y, height, height-bg, I, I-I_bg
             
            for i = 1:size(p,1)
                row_p = p(i,:);
                if (20 <= row_p(1)) && (row_p(1) <= 22) && (20 <= row_p(2)) && (row_p(2) <= 22)
                    p_out = row_p;
                end
            end
         end
         
         %select center of image manually
         function [p_out] = get_peaks_center(obj, r_find, N_img)
            img = fitsread([obj.pname filesep obj.fname{1}],  'Info', obj.info{1}, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [obj.frames(1) obj.frames(N_img)] }); % read first N_img frame                
            stdev_img = std(img,0,3);
            
            f1 = figure('units','normalized','outerposition',[0 0 1 1]);
                imshow(stdev_img,[min(min(stdev_img)) max(max(stdev_img))], 'InitialMagnification', 'fit'); 
                hold on
                [x,y] = ginput(1);
                close gcf
                p_out = [round(x), round(y), 0, 0, 0, 0];
                
         end        
                
         % determine peak-finding threshold and return found peaks for supplied image
         % (more generic than get_h_min)
         function [threshold, p_out] = get_threshold(obj, r_find, img)
             
            p = find_peaks2d(img, r_find, 0, 0); % finding all possible peaks p has x, y, height, height-bg, I, I-I_bg
             
            close all
            figure('units','normalized','outerposition',[0 0 1 1])
            img_mean = mean(img(:));
            img_std = std(img(:));            
             
            p_7std = p(find(p(:,4)>=7*img_std), :); % this is just an estimate, #of peaks found may vary since peak_find algorithm return height as int not double
            p_5std = p(find(p(:,4)>=5*img_std), :);
            p_3std = p(find(p(:,4)>=3*img_std), :);
 
            %plot
            subplot(1, 2, 1)
            imagesc(img), colorbar, axis image, colormap gray, hold on
            if size(p_3std,1)>0
                h(1) = plot(p_3std(:,1)+1, p_3std(:,2)+1, 'ro');
            end
            if size(p_5std,1)>0
                h(2) = plot(p_5std(:,1)+1, p_5std(:,2)+1, 'go');
            end
            if size(p_7std,1)>0
                h(3) = plot(p_7std(:,1)+1, p_7std(:,2)+1, 'bo');
            end
            legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
                   
             
            subplot(1, 2, 2)
            xhist = min(p(:,4)):5:max(p(:,4));
            n = hist(p(:,4), xhist);
            semilogy(xhist, sum(n)-cumsum(n)), hold on
            h(1) = vline(3*img_std, 'r');
            h(2) = vline(5*img_std, 'g');
            h(3) = vline(7*img_std, 'b');
            legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
            set(gca, 'XLim', [0 xhist(end)])
            xlabel('Minimal height'), ylabel('# of peaks found')
            axis square
             
             
            % promp
            options.WindowStyle='normal';
            prompt={'Enter min heigth (default=5*sigma):'};
            def={num2str(round(5*img_std))};
            threshold_dlg = inputdlg(prompt, strcat('Enter threshold:'), 1, def, options);
            threshold = str2double(threshold_dlg(1));
            close all
                          
            p_out = p(find(p(:,4)>=threshold),:);
         end       
          
          
         % generate average image, starting from first frame until N_max
         function [ avg_frame ] = average_image(obj, N_max )
 
            if N_max <= 0 % adjust to full movie lenght
                N_max = obj.mov_length; 
            end
            avg_frame = zeros(obj.sizeX, obj.sizeY);
             
            go_on = 1;
            obj.initRead;
            N = 0;
            while go_on
                [movie, frames, go_on]  = obj.readNext;
             
                if N+length(frames) < N_max
                    avg_frame = avg_frame + sum(movie,3);
                    N = N + length(frames);
                else
                    k = min(N_max-N, length(frames));
                    avg_frame = avg_frame + sum(movie(:,:,1:k),3);
                    N = N + k;
                    go_on = 0;
                end
            end
            avg_frame = avg_frame ./ N; 
         end
          
          % generate average image of last N_avg frames
         function [ avg_frame ] = average_image_last(obj, N_avg)
 
            N_max = obj.mov_length;
             
            avg_frame = zeros(obj.sizeX, obj.sizeY);
             
            go_on = 1;
            obj.counter = length(obj.frames) - N_avg + 1; %start frame for averaging 
            N = 0;
            while go_on
                [movie, frames, go_on]  = obj.readNext;
             
                if N+length(frames) < N_max
                    avg_frame = avg_frame + sum(movie,3);
                    N = N + length(frames);
                else
                    k = min(N_max-N, length(frames));
                    avg_frame = avg_frame + sum(movie(:,:,1:k),3);
                    N = N + k;
                    go_on = 0;
                end
            end
            avg_frame = avg_frame ./ N; 
         end
          
     % generate average image of last N_avg frames for long movies
         function [ avg_frame ] = average_image_lastLong(obj, N_avg)
 
            N_max = obj.mov_length;
             
            avg_frame = zeros(obj.sizeX, obj.sizeY);
             
            go_on = 1;
            obj.counter = length(obj.frames) - N_avg + 1; %start frame for averaging 
            N = 0;
            while go_on
                [movie, frames, go_on]  = obj.readLongNext;
             
                if N+length(frames) < N_max
                    avg_frame = avg_frame + sum(movie,3);
                    N = N + length(frames);
                else
                    k = min(N_max-N, length(frames));
                    avg_frame = avg_frame + sum(movie(:,:,1:k),3);
                    N = N + k;
                    go_on = 0;
                end
            end
            avg_frame = avg_frame ./ N; 
         end
         
         
         % Fits a PSF to one spot in each frame in fit_frames, at given
        % positions fit_pos, using fit parameter sigma
         
        function [ pos_trace ] = fit_psf_to_movie(obj, fit_frames, fit_pos, sigma)
         
            w_fit = 3*sigma;
            s_x = sigma;
            s_y = s_x;
            pos_trace = [];
            counter = 0;
            i=1;
             
            while counter<20 && i<=size(fit_frames,1)
                cur_frame = obj.readFrame(fit_frames(i));     
                x_0 = fit_pos(i,1)+1;
                y_0 = fit_pos(i,2)+1;
                A = cur_frame(y_0,x_0);
                [X,Y] = meshgrid(max(x_0-w_fit,1):min(x_0+w_fit,512),max(y_0-w_fit,1):min(y_0+w_fit,512));
                Z = cur_frame(max(y_0-w_fit,1):min(y_0+w_fit,512),max(x_0-w_fit,1):min(x_0+w_fit,512));
                xyz_data=zeros(size(X,1)*size(X,2),3);
                xyz_data(:,1) = reshape(X, size(X,1)*size(X,2), 1); %make a vector out of matrix X
                xyz_data(:,2) = reshape(Y, size(Y,1)*size(Y,2), 1); %make a vector out of matrix X
                xyz_data(:,3) = reshape(Z, size(Z,1)*size(Z,2), 1); %make a vector out of matrix X
                bg = mean(xyz_data(:,3));
 
                param_init = [x_0 y_0 s_x s_y A-bg bg];
                options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',50000,...
                    'TolFun',1e-9,'MaxIter',1000, 'TolX', 1e-9); 
                [param, chi2, residual, exitflag, output] = lsqcurvefit(@gauss2d_bg, param_init,xyz_data(:,1:2), xyz_data(:,3),[],[], options);
                %display([num2str(output.iterations) ', ' num2str(output.funcCount)])
                if exitflag <= 0
                        display(['WARNING: Fitting gaussian failed. Exitflag: ' num2str(exitflag)])
                        counter = counter + 1;
                        if counter == 20
                            display('Moving on to next spot')
                        end
                elseif counter > 0
                   counter = 0;
                end
                 
            pos_trace = [pos_trace ; [param chi2]];
            i=i+1;
            end
        end
         
      % Fits a PSF to one spot in each frame in fit_frames, at given
      % positions fit_pos, using fit parameter sigma; difference to
      % fit_psf_to movie: uses same value for sigma_y AND sigma_y
      % -> symmetric PSF, one fit parameter less.
         
        function [ pos_trace ] = fit_sym_psf_to_movie(obj, fit_frames, fit_pos, sigma)
         
            w_fit = 3*sigma;
            S = sigma;
            pos_trace = [];
            counter = 0;
            i=1;
             
            while counter<20 && i<=size(fit_frames,1)
                cur_frame = obj.readFrame(fit_frames(i));     
                x_0 = fit_pos(i,1)+1;
                y_0 = fit_pos(i,2)+1;
                A = cur_frame(y_0,x_0);
                [X,Y] = meshgrid(max(x_0-w_fit,1):min(x_0+w_fit,512),max(y_0-w_fit,1):min(y_0+w_fit,512));
                Z = cur_frame(max(y_0-w_fit,1):min(y_0+w_fit,512),max(x_0-w_fit,1):min(x_0+w_fit,512));
                xyz_data=zeros(size(X,1)*size(X,2),3);
                xyz_data(:,1) = reshape(X, size(X,1)*size(X,2), 1); %make a vector out of matrix X
                xyz_data(:,2) = reshape(Y, size(Y,1)*size(Y,2), 1); %make a vector out of matrix X
                xyz_data(:,3) = reshape(Z, size(Z,1)*size(Z,2), 1); %make a vector out of matrix X
                bg = mean(xyz_data(:,3));
 
                param_init = [x_0 y_0 S A-bg bg];
                options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',50000,...
                    'TolFun',1e-9,'MaxIter',1000, 'TolX', 1e-9); 
                [param, chi2, residual, exitflag, output] = lsqcurvefit(@gauss2d_bg_sym, param_init,xyz_data(:,1:2), xyz_data(:,3),[],[], options);
                %display([num2str(output.iterations) ', ' num2str(output.funcCount)])
                if exitflag <= 0
                        display(['WARNING: Fitting gaussian failed. Exitflag: ' num2str(exitflag)])
                        counter = counter + 1;
                        if counter == 20
                            display('Moving on to next spot')
                        end
                elseif counter > 0
                   counter = 0;
                end
                 
            pos_trace = [pos_trace ; [param chi2]];
            i=i+1;
            end
        end
         
      % Fits symmetric PSFs to all spots occuring in one frame in fit_frames,
      % starting at given positions fit_pos, using fit parameter sigma;
      % difference to fit_sym_psf_to movie: processing order is reversed:
      % frame over spot.
      % uses same value for sigma_y AND sigma_y -> symmetric PSF, one fit
      % parameter less than in e.g. fit_psf_to_movie.
      % counter keeping track of misfit frames for a certain spot
      % is abolished here (no sense).
         
        function [ fit_output ] = fit_sym_psfs_to_frame(obj, frame_num, spots_pos, sigma)
         
            w_fit = 3*sigma;
            S = sigma;
            fit_output = zeros(size(spots_pos,1),6);
            frame = obj.readFrame(frame_num);
             
            for i = 1:size(spots_pos,1) 
                if spots_pos(i,1)>=1 && spots_pos(i,2)>=1 && spots_pos(i,1) <= obj.sizeX && spots_pos(i,2)<=obj.sizeY
                x_0 = spots_pos(i,1); % already includes drift correction
                y_0 = spots_pos(i,2); % already includes drift correction
                A = frame(round(y_0),round(x_0));
                [X,Y] = meshgrid(max(round(x_0)-w_fit,1):min(round(x_0)+w_fit,512),max(round(y_0)-w_fit,1):min(round(y_0)+w_fit,512));
                Z = frame(max(round(y_0)-w_fit,1):min(round(y_0)+w_fit,512),max(round(x_0)-w_fit,1):min(round(x_0)+w_fit,512));
                xyz_data=zeros(size(X,1)*size(X,2),3);
                xyz_data(:,1) = reshape(X, size(X,1)*size(X,2), 1); %make a vector out of matrix X
                xyz_data(:,2) = reshape(Y, size(Y,1)*size(Y,2), 1); %make a vector out of matrix Y
                xyz_data(:,3) = reshape(Z, size(Z,1)*size(Z,2), 1); %make a vector out of matrix Z
                bg = mean(xyz_data(:,3));
 
                param_init = [x_0 y_0 S A-bg bg];
                options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',50000,...
                    'TolFun',1e-9,'MaxIter',1000, 'TolX', 1e-9); 
                [param, chi2, residual, exitflag] = lsqcurvefit(@gauss2d_bg_sym, param_init,xyz_data(:,1:2), xyz_data(:,3),[],[], options);
                %display([num2str(output.iterations) ', ' num2str(output.funcCount)])
                if exitflag <= 0
                        display(['WARNING: Fitting gaussian failed. Exitflag: ' num2str(exitflag)])
                end
                 
                fit_output(i,1:6)= [param chi2];
                end
            end
        end
         
         
        % Fits a PSF to one spot in each frame in fit_frames, at given
        % positions fit_pos, using fit parameter sigma; difference to
        % fit_psf_to_movie: uses MATLAB function "fit" with lower/upper
        % bounds for parameters
         
        function [ fit_data ] = fit_psf_2_movie(obj, fit_frames, fit_pos, sigma)
         
            w_fit = 3*sigma;
            s_x = sigma;
            s_y = s_x;
            fit_data = cell(size(fit_frames,1),3);
            counter = 0;
            i=1;
             
            while counter<20 && i<=size(fit_frames,1)
                cur_frame = obj.readFrame(fit_frames(i));     
                x_0 = fit_pos(i,1);
                y_0 = fit_pos(i,2);
                A = cur_frame(y_0,x_0);
                [X,Y] = meshgrid(max(x_0-w_fit,1):min(x_0+w_fit,512),max(y_0-w_fit,1):min(y_0+w_fit,512));
                Z = cur_frame(max(y_0-w_fit,1):min(y_0+w_fit,512),max(x_0-w_fit,1):min(x_0+w_fit,512));
                xyz_data=zeros(size(X,1)*size(X,2),3);
                xyz_data(:,1) = reshape(X, size(X,1)*size(X,2), 1); %make a vector out of matrix X
                xyz_data(:,2) = reshape(Y, size(Y,1)*size(Y,2), 1); %make a vector out of matrix X
                xyz_data(:,3) = reshape(Z, size(Z,1)*size(Z,2), 1); %make a vector out of matrix X
                bg = mean(xyz_data(:,3));
 
                param_init = [x_0 y_0 s_x s_y A bg];
                %options = ['Algorithm','Levenberg-Marquardt','Display','off', 'MaxFunEvals',50000,...
                 %   'TolFun',1e-9,'MaxIter',1000, 'TolX', 1e-9];
                ft = fittype(@(x_0,y_0,s_x,s_y,A,bg,x1,x2)gauss2d_bg([x_0,y_0,s_x,s_y,A,bg], [x1,x2]), ...
                    'coefficient', {'x_0','y_0','s_x','s_y','A','bg'}, 'independent', {'x1','x2'});
                [func, gof, output] = fit(xyz_data(:,1:2), xyz_data(:,3), ft,'StartPoint', param_init,...
                    'Lower', [max(x_0-10,1) max(y_0-10,1) 0 0 0 0], 'Upper', [min(x_0+7,obj.sizeX), min(y_0+7,obj.sizeY), 5 5 inf inf],...
                    'Algorithm','Trust-Region','Display','off',...
                    'MaxFunEvals',50000, 'TolFun',1e-9,'MaxIter',1000, 'TolX', 1e-9);
                %display([num2str(output.iterations) ', ' num2str(output.funcCount)])
                if output.exitflag <= 0
                        display(['WARNING: Fitting gaussian failed. Exitflag: ' num2str(output.exitflag)])
                        counter = counter + 1;
                        if counter == 20
                            display('Moving on to next spot')
                        end
                elseif counter > 0
                   counter = 0;
                end
                 
            fit_data{i,1} = func;
                        fit_data{i,2} = gof;
                                    fit_data{i,3} = output;
            i=i+1;
            end
        end
         
        % VIRTUAL WINDOW CENTER of MASS (VWCM) position estimator per frame
         
        function [pos, deltas_end, Ns, pos_max, v_max, stdevs] = vwcm_in_frame(obj, frame_num, spots_pos, w, epsilon, N_max)
            pos = zeros(size(spots_pos));
            deltas_end = zeros(size(spots_pos,1),1);
            Ns = zeros(size(deltas_end));
            pos_max = zeros(size(pos));
            v_max = zeros(size(deltas_end));
            stdevs = zeros(size(deltas_end));
            img = obj.readLongFrame(frame_num);
             
            for s = 1:size(spots_pos,1)
                if spots_pos(s,1) >=1 && spots_pos(s,2) >=1 && spots_pos(s,1) <= obj.sizeX && spots_pos(s,2) <= obj.sizeY
                % start values
                x = spots_pos(s,1);
                y = spots_pos(s,2);
                 
                % position of maximum pixel and standard deviation              
                X = max(1,round(x)-w):min(obj.sizeX,round(x)+w);
                Y = max(1,round(y)-w):min(obj.sizeY,round(y)+w);
                W = img(Y,X);
                W = W - min(W(:)); %coarse background correction
                 
                [v_max(s), ind] = max(W(:));
                [a, b] = ind2sub(size(W),ind);
                pos_max(s,:) = [b a];
                stdevs(s) = std(W(:));
                 
                % vwcm position estimation algorithm
                cm = [x y];
                delta = epsilon;
                n = 1;
 
                % begin iterations
                while n<=N_max && delta>=epsilon
                     
                X = max(1,round(x)-w):min(obj.sizeX,round(x)+w);
                Y = max(1,round(y)-w):min(obj.sizeY,round(y)+w);
                W = img(Y,X);
                W = W - min(W(:)); %coarse background correction              
 
                %x
                dx = x-round(x);
 
                Sx = sum(W,1);                
                Sx(1) = Sx(1)*((X(1)==1)*0.5+(X(1)>1)*(0.5-dx));
                Sx(end) = Sx(end)*((X(end)==obj.sizeX)*0.5+(X(end)<obj.sizeX)*(0.5+dx));
                 
                X(1) = max(1.25 , mean([x-w X(1)+0.5]));
                X(end) = min(obj.sizeX-0.25 , mean([x+w X(end)-0.5]));
 
                x = (Sx*X')/sum(Sx);
 
                %y
                dy = y-round(y);
 
                Sy = sum(W,2);
                Sy(1) = Sy(1)*((Y(1)==1)*0.5+(Y(1)>1)*(0.5-dy));
                Sy(end) = Sy(end)*((Y(end)==obj.sizeY)*0.5+(Y(end)<obj.sizeY)*(0.5+dy));
 
                Y(1) = max(1.25, mean([y-w Y(1)+0.5]));
                Y(end) = min(obj.sizeY-0.25, mean([y+w Y(end)-0.5]));
 
                y = (Y*Sy)/sum(Sy);
 
                delta = norm([x y] - cm);
                cm = [x y];
                n=n+1;
                end
                 
            deltas_end(s) = delta;
            pos(s,:) = cm;
            Ns(s) = n-1;
            end
            end
        end
    end
     
     
end