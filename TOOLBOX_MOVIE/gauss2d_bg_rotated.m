function [ z ] = gauss2d_bg_rotated(parameter, x) %x = (x_coor, y_coor)
%gauss2d_bg computes a vector z of each high such as x, y, z triples
    x_0 = parameter(1);
    y_0 = parameter(2);
    a = parameter(3);
    b = parameter(4);
    c = parameter(5);
    A = parameter(6);
    bg = parameter(7);
    
    
    z = A * exp(- ( abs(a).*(x(:,1)-x_0).^2 + abs(b).*(x(:,2)-y_0).^2 + 2.*c.* (x(:,1)-x_0).*(x(:,2)-y_0) ) ) + bg;
    

end

