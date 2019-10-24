function [ baseline ] = find_baseline( Intersections, Stroke )
%find_baseline Set a frame of reference parrallel/perpendicular to the cell's surface
%   The idea is to use the PIV measuring lines and go back towards the
%   cell's surface, then find the points of the surface closest to where
%   the cilium is and use them to set the frame of reference

% spacing between piv measuring lines

try
    % this only works if all the lines have the same length (rare)
    lines_dx = mean(diff([Intersections.line_x],1,2),2);
    lines_dy = mean(diff([Intersections.line_y],1,2),2);
catch
    
    % find correspondences (minimum distance) between inmost and second-to-inmost lines
    [idx,~] = ...
        knnsearch([Intersections(2).line_x, Intersections(2).line_y],...
        [Intersections(1).line_x, Intersections(1).line_y]);
    lines_dx = diff([Intersections(1).line_x,Intersections(2).line_x(idx)],1,2);
    lines_dy = diff([Intersections(1).line_y,Intersections(2).line_y(idx)],1,2);
    
    
end %try


% assuming that all videos had 1um spacing between lines
baseline.x = Intersections(1).line_x - lines_dx .* Intersections(1).height_mum;
baseline.y = Intersections(1).line_y - lines_dy .* Intersections(1).height_mum;

% now in general we can not assume that baseline.x and y will be sorted,
% because with the "regression" of the piv lines we might have some points
% out of order. So parametrise the x and y
baseline.t = cumsum(ones(size(baseline.x)));

% then fit them with a polynomial, then use that as cell surface
baseline.x = polyval(polyfit(baseline.t, baseline.x, 4),baseline.t);
baseline.y = polyval(polyfit(baseline.t, baseline.y, 4),baseline.t);

% now find points within 2um (radius) of the average base
baseline_dist_from_base = hypot( baseline.x - Stroke.basal_avg_xx, baseline.y - Stroke.basal_avg_yy) .* Stroke.px2mum;
idx_baseline_close_cilium = baseline_dist_from_base <= 2;
baseline.x_close = baseline.x(idx_baseline_close_cilium);
baseline.y_close = baseline.y(idx_baseline_close_cilium);

% interpolate them with a straight line
baseline.P = polyfit(baseline.x_close, baseline.y_close, 1);

% and find component of the tangential (parallel) versor
baseline.dx = diff(baseline.x_close([1,end]));
baseline.dy = diff(polyval(baseline.P,baseline.x_close([1,end])));
baseline.parav_x = baseline.dx / hypot(baseline.dx, baseline.dy);
baseline.parav_y = baseline.dy / hypot(baseline.dx, baseline.dy);

% now I have to find the perpendicular versor. I want the normal to be
% always defined going "up" from the cell's surface

% first attempt:
baseline.perpv_x = -baseline.parav_y;
baseline.perpv_y =  baseline.parav_x;

 
% apply transformation to average cilium base and intersection with highest
% line
[basal_avg_para, basal_avg_perp] = xy2paraperp(baseline, Stroke.basal_avg_xx, Stroke.basal_avg_yy);
[intersections_para, intersection_perp] = xy2paraperp(baseline, Intersections(end).inter_xx, Intersections(end).inter_yy);

% now check whether TM is the good transformation or not
flag_good_transform = mean(intersection_perp) > basal_avg_perp;

% change transformation if needed
if ~flag_good_transform
    baseline.perpv_x =  baseline.parav_y;
    baseline.perpv_y = -baseline.parav_x;
end %if


end

