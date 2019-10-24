function [para, perp] = xy2paraperp(baseline, x, y)
%xy2parperp takes x,y coordinates of vectors and outputs them in
%parallel/perpendicular components

para = x.*baseline.parav_x + y.*baseline.parav_y;
perp = x.*baseline.perpv_x + y.*baseline.perpv_y;


end

