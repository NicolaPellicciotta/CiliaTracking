function [m, sm] = weightednanmean(x,s,dim)



w = ones(size(s))./(s.^2);
w(s == 0) = 0;    %metto a zero il peso dei valori con errori 0 (i.e. non mediati, altrimenti un errore ce l'avrebbero, a meno che sia media di due cose uguali)
nans = isnan(x);
x(nans) = 0;
w(nans) = 0;

% n = sum(~nans,dim);
% w(n==0) = NaN; % prevent divideByZero warnings
% Sum up non-NaNs, and divide by the number of non-NaNs.
sumweights = sum(w,dim);
m = sum(w.*x,dim) ./ sumweights;
sm = ones(size(m))./ sqrt(sumweights);
end