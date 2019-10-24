function [ xx, yy ] = transf2D( TM, x, y )
%transf2D Applyes the affine transformation matrix to the vectors x and y


if ~isvector(x) || ~isvector(y)
    error('x and y have to be vectors');
end

if numel(x) ~= numel(y)
    error('x and y have to be of the same size');
end

if ~ismatrix(TM)
    error('TM has to be a 2D matrix');
end

if size(TM,1) ~= size(TM,2)
    error('TM has to be a square matrix');
end

if any(size(TM) > 3) || any(size(TM) < 2)
    error('TM can aeither be a 2-by-2 or 3-by-3 matrix');
end


flag_wascolumn_x = iscolumn(x);
flag_wascolumn_y = iscolumn(y);

% put in row if needed
if flag_wascolumn_x, x = x(:)'; end
if flag_wascolumn_y, y = y(:)'; end

if size(TM,1) == 3 % affine transformation
    
    IM = vertcat(x, y, ones(1,numel(x))); %input matrix
    OM = TM * IM;   % output matrix
    
    xx = OM(1,:);
    yy = OM(2,:);
    
    
elseif size(TM,1) == 2 % only simmetries and rotations about origin, no shifting or scaling allowed
    
    IM = vertcat(x, y); %input matrix
    OM = TM * IM;   % output matrix
    
    xx = OM(1,:);
    yy = OM(2,:);
    
else
    error('TM is of the wrong size')
end %if


% return proper size
if flag_wascolumn_x, xx = xx(:); end
if flag_wascolumn_y, yy = yy(:); end


end

