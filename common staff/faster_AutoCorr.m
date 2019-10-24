function [res, lags] = faster_AutoCorr(cc, ll)
% only works on column vectors now
%%

% cc = rand(10,1)
% ll=5

ind_M = [0:ll-1] + cumsum(ones(size(cc)));
ccounts = ind_M <= numel(cc);
ind_M = mod(ind_M-1, numel(cc))+1;


M = cc(ind_M);
M = M.*ccounts;

res = cc .* M ;
res = sum(res)./sum(ccounts) - mean(cc)^2;
res = res ./ var(cc);

lags = 0:ll-1;

