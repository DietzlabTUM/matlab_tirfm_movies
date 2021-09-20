function [ result, err, ci, area, residuals ] = fit_gauss2d_bg( x0, y0, sigma, w_fit, img ) % (x0, y0, s, w, image)
%UNTITLED2 Summary of this function goes here
%   c_init = [x_0 y_0 s A-bg bg];
    
    options = optimset('Algorithm','levenberg-marquardt','Display','off', 'MaxFunEvals',50000,'TolFun',1e-9,'MaxIter',200, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',
       %options = optimset('TolFun',1e-9,'MaxIter',200, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',

    x0 = round(x0);
    y0 = round(y0);
   
    if x0 < 1 
        x0 = 1;
    end
    if x0 > 512
        x0 = 512;
    end
    if y0 < 1 
        y0 = 1;
    end
    if y0 > 512
        y0 = 512;
    end    
   
    
    area = [  max(x0-w_fit,1) , min(x0+w_fit,size(img,1)) , max(y0-w_fit,1) , min(y0+w_fit,size(img,2)) ];
    [X,Y] = meshgrid(area(1):area(2) , area(3):area(4) );
    Z = img(area(3):area(4), area(1):area(2));
    xyz_data=zeros(size(X,1)*size(X,2),3);
    xyz_data(:,1) = reshape(X, size(X,1)*size(X,2), 1); %make a vector out of matrix X
    xyz_data(:,2) = reshape(Y, size(Y,1)*size(Y,2), 1); %make a vector out of matrix X
    xyz_data(:,3) = reshape(Z, size(Z,1)*size(Z,2), 1); %make a vector out of matrix X
    
    
    % initial parameters
    bg = mean(xyz_data(:,3));
 
    c_init = [x0 y0 sigma img(y0,x0)-bg  bg];
    
  
    % fit the parameters
    [coef,residuals,J,COVB,mse, ErrorModelInfo] = nlinfit(  xyz_data(:,1:2),  xyz_data(:,3), @gauss2d_samesigma_bg, c_init, options);
    
    
    coef(3) = abs(coef(3));

    % calculate errors
    ci = nlparci(coef,residuals,'Jacobian',J);
    
    % caculate c_err
    err = coef-ci(:,1)'+ci(:,2)'-coef;

    result = coef(1:5);
    
    %disp(ErrorModelInfo)
    
end