function [cc, xb, fit, pd, hb, hp] = my_lognormalfit(data, bins, flag_show_figure, ha, color_histogram, color_fit)

%bins can be either a scalar or an array
if nargin < 6, color_fit = 'r'; end
if nargin < 5, color_histogram = 'b'; end
if nargin < 4, hf = figure; ha = gca; end
if nargin < 3, flag_show_figure = 0; end
if nargin < 2, bins = 10; end

[cc,xb] = hist(data,bins);
parmhat = lognfit(data);
pd = ProbDistUnivParam('lognormal',parmhat);
fit = pdf(pd,xb);
%fit = fit.*trapz(xb,cc);
% ncc = cc./trapz(xb,cc);
ncc = cc./(sum(cc)*mean(diff(xb)));

if flag_show_figure == 1
    hold(ha,'on');
    hb = bar(ha,xb,ncc,1,color_histogram);
    set(get(hb,'children'),'FaceAlpha',.5,'EdgeAlpha',0);

    hp = plot(ha,xb,fit,'-','LineWidth',2,'Color',color_fit);

end

end
