function [x0, y0] = find_new_trajectory_centre(x,y, flag_debug)
% find_new_trajectory_centre finds if possible a point so that the angles
% off the line with it are on a monotonous function

if nargin < 3 || isempty(flag_debug)
    flag_debug = false;
end

x = x(:);
y = y(:);

% put x and y in new units, everything has to be positive and multiplied by
% 100;
sf = 100;                   % scale factor
x_off = min(x.*sf) - 1;      % x offset
y_off = min(y.*sf) - 1;      % y offsset

% new coordinates
xx = x.*sf - x_off;
yy = y.*sf - y_off;

% create mask full within the curve
mask = poly2mask(xx,yy,ceil(max(yy)),ceil(max(xx)));
mask = imerode(mask,strel('disk',5));

% create grid to cover curve
[x_grid, y_grid] = meshgrid(1:ceil(max(xx)), 1:ceil(max(yy)));

x_grid = x_grid(mask);
y_grid = y_grid(mask);


% could do a for loop but very slow i think. matrix much better

% matrix of point - centre, size is n_centres x n_points_curve
dx = xx' - x_grid;
dy = yy' - y_grid;

% find angles between all points along the line and the centres
phi = unwrap(atan2(dy,dx),pi/2,2);

% now differentiate the angle along the line and take the sign of the
% derivative
dphi = diff(phi,1,2);

% find good centres: derivative has to ba always positive or always negative
idx_good_centres = all(dphi > 0,2) | all(dphi < 0, 2);
goodness_centres = min(abs(dphi),[],2);

% if there are many good centres, find the one where the minimum of the
% |derivative| is highest
if sum(idx_good_centres) > 1
    goodness_centres(~idx_good_centres) = 0;
    [~,idx_best_centre] = max(goodness_centres);
else
    idx_best_centre = find(idx_good_centres);
end


% so this is the best point
xx0 = x_grid(idx_best_centre);
yy0 = y_grid(idx_best_centre);

% now let's put it back in normal coordinates
x0 = (xx0 + x_off) ./ sf;
y0 = (yy0 + y_off) ./ sf;



if flag_debug
    % plot of curves
    hf = figure; clf
    hf.PaperPositionMode = 'auto';
    hf.Units = 'centimeters';
    hf.Position([3 4]) = [16 6];
    
    subplot(121)
    plot(x,y)
    hold on
    plot(x0,y0,'ro');
    axis tight
    
    subplot(122)
    imagesc(mask)
    hold on
    plot(xx,yy)
    plot(x_grid(~idx_good_centres),y_grid(~idx_good_centres),'r.')
    plot(x_grid(idx_good_centres),y_grid(idx_good_centres),'y.')
    plot(xx0,yy0,'go')
    set(gca,'YDir','normal')
    
end

