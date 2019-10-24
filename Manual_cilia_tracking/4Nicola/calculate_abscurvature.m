function [k] = calculate_abscurvature(x,y,t)
% curvature of parametric curve, equation from wikipedia
    
dx = gradient(x)./gradient(t);
ddx = gradient(dx)./gradient(t);

dy = gradient(y)./gradient(t);
ddy = gradient(dy)./gradient(t);

num = abs( dx.*ddy - ddx.*dy );
den = (dx.*dx + dy.*dy);
den = sqrt(den);
den = den .* den .* den;

k = num ./ den;


end