function [ criteria ] = filter_spots( sigma, ratio_threshold, sigma_filter)
% sigma: sigma_x and sigma_y of spots
% ratio_threshold: ratio min, max
% 

spotsize_mean = mean(sqrt( sigma(:,1).^2 +  sigma(:,2).^2));
spotsize_std = std(sqrt( sigma(:,1).^2 +  sigma(:,2).^2));

% determine spotsize threshold
close all

xhist = 0:0.5:5; %pixel
n = hist(sqrt( sigma(:,1).^2 +  sigma(:,2).^2), xhist);
bar(xhist, n );hold on
xlabel('$\sqrt{(\Delta x)^2 + (\Delta y)^2}$ in pixel', 'Interpreter', 'LaTex')
ylabel('Frequency')


spotsize_thresh = [spotsize_mean-sigma_filter(1)*spotsize_std  spotsize_mean+sigma_filter(2)*spotsize_std];
vline(spotsize_thresh(1), 'g');
vline(spotsize_thresh(2), 'g');



criteria = ones(size(sigma,1),2 );
criteria(:,1) = [sigma(:,1)./sigma(:,2) < ratio_threshold(2) & sigma(:,1)./sigma(:,2) > ratio_threshold(1)];
criteria(:,2) = [sqrt( sigma(:,1).^2 +  sigma(:,2).^2 ) > spotsize_thresh(1) & sqrt( sigma(:,1).^2 +  sigma(:,2).^2 ) < spotsize_thresh(2)];


end

