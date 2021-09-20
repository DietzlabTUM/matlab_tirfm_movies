function [ z ] = gauss2d_bg_sym(parameter, x) %x = (x_coor, y_coor)
%gauss2d_bg_sym computes a vector z of each high such as x, y, z triples
%sigma_x is equal to sigma_y, hence symm(etric)
    x_0 = parameter(1);
    y_0 = parameter(2);
    S = parameter(3);
    A = parameter(4);
    bg = parameter(5);
    z = A * exp(- ( (x(:,1)-x_0).^2./(2*S^2) + (x(:,2)-y_0).^2./(2*S^2) ) ) + bg;

end

