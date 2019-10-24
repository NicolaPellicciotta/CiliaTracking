function [hs] = plotc(x,y,c, LineWidth, flag_flat)

if nargin < 5 || isempty(flag_flat)
    flag_flat = false;
end

if nargin < 4 || isempty(LineWidth)
    LineWidth = 2;
end


if nargin < 3 || isempty(c)
    c = cumsum(ones(size(y)));
end

if flag_flat
    edgecol = 'flat';
else
    edgecol = 'interp';
end

z = zeros(size(x));

x = x(:)';
y = y(:)';
z = z(:)';
c = c(:)';


hs = surface([x;x],[y;y],[z;z],[c;c],...
    'facecol','no',...
    'edgecol',edgecol,...
    'linew',LineWidth);

end