%%% Function that scales all values of an NxN matrix in reference of an
%%% upper and lower bound value. Useful for scaling greyscale values of a
%%% number of images individually.


function [scmat] = scalematrix(mat, lb, ub);

scmat = double(mat-min(mat(:)));
scmat = scmat./max(scmat(:)).*(ub-lb);
scmat = floor(scmat+lb);

end
    