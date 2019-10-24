function [ rows, cols ] = rect2sub( rect_output )
%rect2sub takes a rect in pixels in the form [xmin ymin xwidth yheight] and
%gives the subscript of the matrix image.
%   Example: 
%   imshow(IM,[]);
%   rect_output = getrect(gcf);
%   imshow(IM(rect2sub(rect_output))); %shows the roi

rect_output = round(rect_output);

xmin    = rect_output(1);
ymin    = rect_output(2);
xwidth  = rect_output(3);
yheight = rect_output(4);
rows = (ymin : ymin+yheight);
cols = (xmin : xmin+xwidth);

end

