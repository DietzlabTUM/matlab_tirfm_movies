% This and other related scripts were written by Dr. Matthias Schickinger, Dr. Philip Ketterer, Dr. Jonas Funke,
% Anna-Katharina Pumm, Dr. Wouter Engelen and Eva Bertosin

%% load spots

%change to correct path

load('/Users/username/Documents/green01/20xx-xx-xx_analysis/data.mat')
cd '/Users/username/Documents/green01/20xx-xx-xx_analysis/spots'

%% calculate std dev

fnames = dir('*.tif');
fnamesfits = dir('*.fits');
numfids = length(fnames);
dev_img = zeros (41, 41, numfids);
vals = cell(41, 41, numfids);
vals_norm = cell(1,numfids);
N_peaks=numfids;

for K = 1:numfids
    vals = fitsread(fnamesfits(K).name);
    sz = size(vals);
    dev_img(1:sz(1),1:sz(2),K) = std(vals,0,3);
    K
end

save dev_img
disp('Std dev done')
