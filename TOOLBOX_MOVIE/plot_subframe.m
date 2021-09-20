function [] = plot_subframe( img, x, y, w )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %check if x and y are inside of img
    if x < 1
        display(['Warning: x changed from ' num2str(x) ' to 1.'])
        x=1;
    end
    if x > size(img,1)
        display(['Warning: x changed from ' num2str(x) ' to ' num2str(size(img,1))])
        x=size(img,1);
    end
    if y < 1
        display(['Warning: y changed from ' num2str(y) ' to 1.'])
        y=1;
    end
    if y > size(img,2)
        display(['Warning: y changed from ' num2str(y) ' to ' num2str(size(img,2))])
        y=size(img,1);
    end
    colormap gray
    corner = ([max(1, x-w) min( x+w , size(img,1) ) max(1, y-w) min( y+w , size(img,2) )   ]);
    corner_int = int16([max(1, x-w) min( x+w , size(img,1) ) max(1, y-w) min( y+w , size(img,2) )   ])  ;
    scale = [min(min(img(corner_int(3):corner_int(4), corner_int(1):corner_int(2)))) max(max(img(corner_int(3):corner_int(4), corner_int(1):corner_int(2)))) ];
    imagesc(img, [scale(1) scale(2)]), colorbar
    set(gca, 'XLim', corner(1:2))
    set(gca, 'YLim', corner(3:4))
    
end

