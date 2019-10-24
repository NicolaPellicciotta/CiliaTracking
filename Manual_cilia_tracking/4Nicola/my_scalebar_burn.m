function[image] = my_scalebar_burn(image, real_length, px2realunits, position, grey_color)
%this function actually replaces pixels of the image with the scalebar.
%Only to use when making videos

%% parsing input

if nargin < 4
    position = [0.9 0.1];
    grey_color = 1-1e-15;
end

if nargin < 5
    grey_color = 'w';
end

if ~strcmp(grey_color,'w') && ~strcmp(grey_color,'k')
    error('Color can only be ''w'' or ''k''.');
end

%% checking input

if ~any(strcmp({'uint8';'uint16';'uint32';'uint64'}, class(image)))
    error('Image has to belong to an unsigned integer class.');
end

YPosPx = position(1)*size(image,1);
XPosPx = position(2)*size(image,2);

LengthPx = real_length/px2realunits;
ThicknessPx = 0.01 * size(image,1);

x = round(XPosPx : XPosPx+LengthPx); %at a certain point there was a minus instead of a plus ??
x = x(:)'; %easier to think of it as a line vector
y = round(YPosPx : YPosPx+ThicknessPx);

if strcmp(grey_color,'w')
    image(y,x) = max(image(:));
else
    image(y,x) = min(image(:));
end



