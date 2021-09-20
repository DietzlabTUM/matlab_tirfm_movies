function [ result, err, ci, area ] = fit_gauss2d_mainaxis_bg( x0, y0, sigma, w_fit, img ) % (x0, y0, s, w, image)
%UNTITLED2 Summary of this function goes here
%   c_init = [x_0 y_0 s_x s_y A-bg bg];
    alpha = 0;
    options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',50000,'TolFun',1e-9,'MaxIter',50000, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',

    area = [  max(x0-w_fit,1) , min(x0+w_fit,size(img,1)) , max(y0-w_fit,1) , min(y0+w_fit,size(img,2)) ];
    [X,Y] = meshgrid(area(1):area(2) , area(3):area(4) );
    Z = img(area(3):area(4), area(1):area(2));
    xyz_data=zeros(size(X,1)*size(X,2),3);
    xyz_data(:,1) = reshape(X, size(X,1)*size(X,2), 1); %make a vector out of matrix X
    xyz_data(:,2) = reshape(Y, size(Y,1)*size(Y,2), 1); %make a vector out of matrix X
    xyz_data(:,3) = reshape(Z, size(Z,1)*size(Z,2), 1); %make a vector out of matrix X
    
    bg = mean(xyz_data(:,3));
    
    a = 1./(2*sigma.^2); %circle
    b = a;  %circle
    c = 0; %circle
    
    c_init = [x0 y0 a b c img(y0,x0)-bg bg];
    
    
    
    % fit 
    [coef,r,J,COVB,mse] = nlinfit(  xyz_data(:,1:2)  ,  xyz_data(:,3), @gauss2d_bg_rotated, c_init, options);

    coef(3) = abs(coef(3));
    coef(4) = abs(coef(4));

    % calculate errors
    ci = nlparci(coef,r,'Jacobian',J);
    
    % caculate c_err
    err = coef-ci(:,1)'+ci(:,2)'-coef;
    
    % fitted parameters
    a = coef(3); b = coef(4); c = coef(5);
    
    
    % determine the width and angle
    [V, D] = eig([a c; c b]);
    
    sigma1 = 1 / sqrt(2*abs(D(1,1))); % lambda = 1 / 2 sigma^2 => sigma = 1 / sqrt...
    sigma2 = 1 / sqrt(2*abs(D(2,2)));
    theta = acos(V(1,1));
    result = [coef(1:2) sigma1 sigma2 theta coef(6:7)];
    
    % caculate errors with gaussian error propagation
    tmp = sqrt((a+b).^2/4 - a.*b+c.^2);
    
    %lambda1
    a1 = 0.5*(1-(a-b)/tmp); % dlambda1 / da
    a2 = 0.5*(1-(b-a)/tmp); % dlambda1 / db
    a3 = -c / tmp; % dlambda1 / d
    dlambda1 = sqrt( a1.^2 * err(3).^2 + a2.^2 * err(4).^2 + a3.^2 * err(5).^2 );
    
    %lambda2
    a1 = 0.5*(1+(a-b)/tmp); % dlambda1 / da
    a2 = 0.5*(1+(b-a)/tmp); % dlambda1 / db
    a3 = +c / tmp; % dlambda1 / dc
    dlambda2 = sqrt( a1.^2 * err(3).^2 + a2.^2 * err(4).^2 + a3.^2 * err(5).^2 );
    
    %theta
    dtheta = 0;
    
    err(3:5) = [dlambda1 dlambda2 dtheta];
    
    
end

