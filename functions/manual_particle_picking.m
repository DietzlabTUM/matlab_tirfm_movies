function [pos_ch1, N_peaks] = manual_particle_picking(i, stdev_img, mag)

width = size(stdev_img{i},1)/mag;
height = size(stdev_img{i},2)/mag;

f1 = figure('units','normalized','outerposition',[0 0 1 1]);
n = 1;

for j = 1:mag
    for k = 1:mag
        %show image and set limits according to magnification
        hold on
        imshow(-stdev_img{i},[min(min(-stdev_img{i}))*0.4 max(max(-stdev_img{i}))], 'InitialMagnification', 'fit'); hold on
        xlim([(j-1)*width+1 j*width])
        ylim([(k-1)*height+1 k*height])
        
        %select rotors, press return when all are selected
        while true
            [x, y, button] = ginput(1);
            if isempty(x);
                break;
            end
            
            if button == 3;
                x_n(n-1) = [];
                y_n(n-1) = [];
                n = length(x_n)+1;
                hold off
                imshow(-stdev_img{i},[min(min(-stdev_img{i}))*0.4 max(max(-stdev_img{i}))], 'InitialMagnification', 'fit');
                xlim([(j-1)*width+1 j*width])
                ylim([(k-1)*height+1 k*height])
                hold on
            end
            
            if button == 1;
                x_n(n) = x; % save all points you continue getting
                y_n(n) = y;
                n = n+1;
            end
            
            hold on
            plot(x_n, y_n, 'or', 'MarkerSize', mag*20)
            drawnow
        end
    end
end

pos_ch1 = zeros(length(x_n), 5);
pos_ch1(:,1:2) = [x_n', y_n'];
N_peaks = n;

hold off
close gcf