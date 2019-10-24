function [ordered_x, ordered_y, idx] = sortbypath_pdist(x,y)

if ~isvector(x) || ~isvector(y)
    error('function only works with vectors');
end

x = x(:);
y = y(:);

% figure;
% plot(x,y,'o-')

data = [x,y];

% create matrix of all pairwise distances
D = squareform(pdist(data,'euclidean'));

% find the points the most far from all the others
dist = sum(D,2);
[~,ext] = max(dist);

% initialise result vector
idx = nan(numel(x),1);

% first point is the "extreme" point found a few lines above
idx(1) = ext;

% delete from the distance matrix the column corresponding to the first
% point
D(:,idx(1)) = NaN;

% greedy nearest-neighbour algorithm

for i = 2:numel(x)
    
    % find smaller distance among points not yet touched
    [~,idx(i)] = min(D(idx(i-1),:));
    
    % hence delete corresponding column from distance matrix
    D(:,idx(i)) = NaN;
    
end


% sort vecors according to new order (stored in idx)
ordered_x = x(idx);
ordered_y = y(idx);

end



%% old version, using knnsearch. It works, but knnsearch is order O(N^2)
% 
% function [ordered_x, ordered_y, idx] = sortbypath(x,y)
% 
% if ~isvector(x) || ~isvector(y)
%     error('function only works with vectors');
% end
% 
% x = x(:);
% y = y(:);
% 
% % figure;
% % plot(x,y,'o-')
% 
% data = [x,y];
% [Idx,D] = knnsearch(data,data,'K',numel(x));
% dist = sum(D,2);
% [~,ext] = max(dist);
% 
% idx = nan(numel(x),1);
% 
% idx(1) = ext;
% idx(2) = Idx(ext,2);
% 
% for i=3:numel(x)
%     
%     % while loop to find the closest neighbour not yet in the path
%     j = 2;
%     while any( Idx(idx(i-1),j) == idx(1:i-1) )
%         j = j+1;
%     end
%     idx(i) = Idx(idx(i-1),j);
%     
% end
% 
% 
% ordered_x = x(idx);
% ordered_y = y(idx);
% 
% end