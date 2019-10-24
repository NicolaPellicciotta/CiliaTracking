function [normChi2] = lognormalChi2 (data, bins, plotflag)
% LOGNORMALCHI2 fits the logarithm of "data" with a normal ditribution with
% "bins" bins (or using directly "bins" if it's ana array) normalize the 
%histogram and calculates the Chi2, then divides it by the number of bins
%and returns the normalized value (Chi2/n)

if nargin<3
    plotflag = false;
end

if min(size(data)) > 1
    disp('Data has to be an array. Returning')
    return
end

if size(data,2) ~= 1
    data = data';
end

logdata = log(data);
[oc,xc] = hist(logdata,bins);
pd = fitdist(logdata,'normal');
ec = pdf('normal',xc,pd.mu,pd.sigma);
nec = ec.*trapz(xc,oc);

Chi2 = sum(((oc-nec).^2)./ec);
normChi2 = Chi2/length(xc);

if plotflag == true
    f = figure;
    p1 = plot(xc,oc,'o');
    hold on
    p2 = plot(xc,nec,'r');
    set(p2,'LineWidth',1);
    ylabel('Counts, normalized')
    legend('log(data) distr, normalized',['fit, \chi^{2}/n = ',num2str(normChi2,'%0.3e')],'Location','NorthEast');
end

end