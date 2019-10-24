function [] = clclk_curv_plot(clclk)


ind_clicked = [clclk.clicked_frames.frame_number];

normr_array = [clclk.points(ind_clicked).normr];
fluct_normr_array = normr_array - reshape(smooth(normr_array,5),1,[]);
mm = median(fluct_normr_array);
ss = iqr(fluct_normr_array);
ind_clicked(fluct_normr_array > mm + 3*ss) = [];

nplots = 2;

hf = figure;
clf
hf.Position = [120 300 560*nplots 420];

ls = 0.4;
is = (1-nplots*ls)/(nplots+1);

for i = 1:nplots
    
    ha(i) = axes;
    ha(i).Units = 'normalized';
    ha(i).Box = 'on';
    ha(i).Position = [i*is + (i-1)*ls, 0.15, ls, 0.8];
    
end %for


%% cilium shape plot

axes(ha(1))
ha(1).NextPlot = 'add';
pbar = ha(1).PlotBoxAspectRatio;
ha(1).DataAspectRatio = [1 1 1];
ha(1).PlotBoxAspectRatio = pbar;
ha(1).YDir = 'reverse';

for i = ind_clicked
    plotc(clclk.points(i).cilium_xx, clclk.points(i).cilium_yy, clclk.points(i).curvature ./ clclk.px2mum, 3);
end %for

ha(1).XTick = [];
ha(1).YTick = [];

ha(1).XLim = [ min(arrayfun(@(i)min(clclk.points(i).cilium_xx), ind_clicked)),...
    max(arrayfun(@(i)max(clclk.points(i).cilium_xx), ind_clicked)) ];

ha(1).XLim = ha(1).XLim + .05 * diff(ha(1).XLim) .* [-2 1]; %this gives nice white border before cilia hit colorbar

if min(arrayfun(@(i)min(clclk.points(i).cilium_yy), ind_clicked)) < min(ha(1).YLim) ||...
        max(arrayfun(@(i)max(clclk.points(i).cilium_yy), ind_clicked)) > max(ha(1).YLim)
    
    ha(1).YLim = [ min(arrayfun(@(i)min(clclk.points(i).cilium_yy), ind_clicked)),...
        max(arrayfun(@(i)max(clclk.points(i).cilium_yy), ind_clicked)) ];
    
    ha(1).YLim = ha(1).YLim + .05 * diff(ha(1).YLim) .* [-1 1];
    
    
end %if

ha(1).XLimMode = 'auto';
ha(1).PlotBoxAspectRatio = pbar;


hc(1) = colorbar;
hc(1).Location = 'manual';
hc(1).Position = ha(1).Position .* [1 1 .03 1];
hc(1).Position(1) = ha(1).Position(1) - hc(1).Position(3);
hc(1).Label.String = 'Curvature, [\mum^{-1}]';
hc(1).Label.FontSize = 16;


%% straightened cilium plot

axes(ha(2))
ha(2).NextPlot = 'add';
pbar = ha(2).PlotBoxAspectRatio;


for i = ind_clicked
    plotc(clclk.clicked_frames(i).timestamp .* ones(size(clclk.points(i).curvature)),...
        clclk.points(i).tt .* clclk.px2mum,...
        clclk.points(i).curvature ./ clclk.px2mum, 8);
end %for
% hs = scatter(straightened_ii,straightened_tt,30,straightened_kk,'filled');

xlabel(ha(2),'Time, [s]','FontSize',16);
ylabel(ha(2),'Arclength, [\mum]','FontSize',16);


%%%%%%%%%%%% new figure: autocorrelation of curvature %%%%%%%%%%%%


hf2 = figure;

ha2 = axes;

% imagesc(dt, ds, curv_xcorr2);
hs = surf(clclk.Curv.dt, clclk.Curv.ds, clclk.Curv.curv_xcorr2);
hs.EdgeColor = 'none';
xlabel(ha2, '\tau, [s]','FontSize',16)
ylabel(ha2, '\Deltas, [\mum]','FontSize',16)

hc2 = colorbar;
hc2.Label.String = 'C_\kappa(\Deltas,\tau)';
hc2.Label.FontSize = 16;

hf3 = figure;

ha3 = axes;

% imagesc(dt, ds, curv_xcorr2);
hs = surf(clclk.Curv.dt_rel, clclk.Curv.ds_rel, clclk.Curv.curv_xcorr2);
hs.EdgeColor = 'none';
xlabel(ha3, '\tau, [cycles]','FontSize',16)
ylabel(ha3, '\Deltas, [cilium length]','FontSize',16)

hc2 = colorbar;
hc2.Label.String = 'C_\kappa(\Deltas,\tau)';
hc2.Label.FontSize = 16;




hf5 = figure;
ha(5) = gca;
him5 = imagesc(minmax([clclk.clicked_frames.timestamp]), minmax(clclk.Curv.ntt) * clclk.px2mum,...
    clclk.Curv.int_curv_mat_1omum);
him5.AlphaData = ~isnan(clclk.Curv.int_curv_mat_1omum);
set(ha(5),'YDir','normal')
xlabel(ha(5),'Time, [s]','FontSize',16);
ylabel(ha(5),'Arclength, [\mum]','FontSize',16);
ha(5).YLim = ha(2).YLim;
ha(5).CLim = ha(1).CLim;
hold on
plot([clclk.clicked_frames(ind_clicked).timestamp], max(ha(5).YLim) .* ones(size(ind_clicked)),...
    'rv','MarkerSize',4,'MarkerFacecolor','r');

hc5 = colorbar;
hc5.Label.String = 'Curvature, [\mum^{-1}]';
hc5.Label.FontSize = 16;


hf6 = figure;
ha(6) = gca;
him6 = imagesc(minmax([clclk.clicked_frames.timestamp]), minmax(clclk.Curv.ntt_cl) * clclk.px2mum,...
    clclk.Curv.int_curv_mat_commonlength_1omum);
him6.AlphaData = ~isnan(clclk.Curv.int_curv_mat_commonlength_1omum);
set(ha(6),'YDir','normal')
xlabel(ha(6),'Time, [s]','FontSize',16);
ylabel(ha(6),'Arclength, [\mum]','FontSize',16);
ha(6).CLim = ha(1).CLim;
hold on
plot([clclk.clicked_frames(ind_clicked).timestamp], max(ha(6).YLim) .* ones(size(ind_clicked)),...
    'rv','MarkerSize',4,'MarkerFacecolor','r');

hc6 = colorbar;
hc6.Label.String = 'Curvature, [\mum^{-1}]';
hc6.Label.FontSize = 16;

end