function [clclk] = clclk_curvature_map(clclk_filename, flag_plot, flag_force_userIDbasetip, flag_save, flag_debug)

if nargin < 1 || isempty(clclk_filename)
    % if isempty(clclkfilename)
    [clclk_filename, clclk_pathname] = uigetfile('*.clclk','Select a .clclk file.');
else
    [clclk_filename, clclk_pathname] = parse_filename(clclk_filename);
end

if nargin < 2 || isempty(flag_plot)
    flag_plot = 1;
end %if

if nargin < 3 || isempty(flag_force_userIDbasetip)
    flag_force_userIDbasetip = 1;
end %if

if nargin < 4 || isempty(flag_save)
    flag_save = true;
end

if nargin < 5 || isempty(flag_debug)
    flag_debug = false;
end

% %% debug inputs
% cch
% warning('I''m debugging');
% clclk_filename = 'test.22Jul2015_16.14.46.clclk';
% % clclk_filename = 'test.22Jul2015_15.58.09.clclk';
% [clclk_filename, clclk_pathname] = parse_filename(clclk_filename);

%% read clclk file and order and interpolate clicked points

[clclk] = clclk_reader(fullfile(clclk_pathname, clclk_filename),0,'blue_gradient',flag_force_userIDbasetip);
% close all
uiwait(gcf)
if ~flag_plot
    close(gcf)
end

%% remove misclicking (because this is only done as a display thing in clclk_reader)

ind_clicked = [clclk.clicked_frames.frame_number];

normr_array = [clclk.points(ind_clicked).normr];
fluct_normr_array = normr_array - reshape(smooth(normr_array,5),1,[]);
mm = median(fluct_normr_array);
ss = iqr(fluct_normr_array);
ind_clicked(fluct_normr_array > mm + 3*ss) = [];


%% straightened cilia with color-coded curvature

straightened_xx = [];
straightened_yy = [];
straightened_kk = [];
straightened_tt = [];
straightened_ii = [];

for i = ind_clicked
    
    straightened_xx = vertcat(straightened_xx,clclk.points(i).cilium_xx(:));
    straightened_yy = vertcat(straightened_yy,clclk.points(i).cilium_yy(:));
    straightened_tt = vertcat(straightened_tt,clclk.points(i).tt(:));
    straightened_kk = vertcat(straightened_kk,clclk.points(i).curvature(:));
    
    straightened_ii = vertcat(straightened_ii,i(ones(size(clclk.points(i).curvature(:)))));
end %for

% figure
% % colormap redblue
% scatter(straightened_ii,straightened_tt,30,straightened_kk,'filled')
% shg


%% interpolate data for curvature map


% first create a common tt (as long as the longest snapshot
max_tt = arrayfun(@(i)clclk.points(i).tt(end),ind_clicked); % end arclength of cilia snapshots, need for interpolating the end of the missing data
[~,wh] = max(max_tt);
ntt = clclk.points(ind_clicked(wh)).tt;    % common tt, will use for curvature matrix


% index on all the points between first and last clicked
int_ind_clicked = linspace(min(ind_clicked), max(ind_clicked),range(ind_clicked)+1);

% interpolate to find endpoint of missing snapshots
int_max_tt = interp1(ind_clicked, max_tt, int_ind_clicked ,'cubic' );

% Interpolate Px and Py to "fix" missing data
matPx = vertcat(clclk.points(ind_clicked).cilium_Px)';
matPy = vertcat(clclk.points(ind_clicked).cilium_Py)';
intPx = interp1(ind_clicked, matPx', int_ind_clicked ,'cubic' )';
intPy = interp1(ind_clicked, matPy', int_ind_clicked ,'cubic' )';

% now create the equally-spaced (time-wise) curvature matrix, using
% reconstructed snapshots
int_curv_mat = nan(numel(ntt), numel(int_ind_clicked));     % preallocation

if flag_plot
    figure;
    hold on
end %if

for i = 1:numel(int_ind_clicked)
    
    t = ntt(ntt <= int_max_tt(i));
    x = polyval(intPx(:,i), t);
    y = polyval(intPy(:,i), t);
    
    if flag_plot,  plot(x,y); end
    
    %     int_curv_mat(1:numel(t),i) = calculate_signedcurvature(x, y, t);
    int_curv_mat(1:numel(t),i) = calculate_abscurvature(x, y, t);
end %for

int_curv_mat_1omum = int_curv_mat ./ clclk.px2mum; %go in real world units


% ask the user how much of a cycle this is
answer = Inf;
while ~isinrange(answer, [0, 10]) 
    answer = inputdlg('what fraction of a cycle is this?', 'fraction of cycle');
    answer = str2num(answer{1});
end
cycle_completion = answer;

% crop the matrix at shortest cilium length
int_curv_mat_commonlength = int_curv_mat(~any(isnan(int_curv_mat),2),:);
int_curv_mat_commonlength_1omum = int_curv_mat_commonlength ./ clclk.px2mum; %go in real world units

% calculate autocorrelation (at positive lags)
curv_xcorr2 = my_normxcorr2(int_curv_mat_commonlength_1omum);
curv_xcorr2 = curv_xcorr2(1:floor(end/2), 1:floor(end/2));

% calculate lags
ds = ntt(1:size(curv_xcorr2,1)) * clclk.px2mum;     % space lag
dt = interp1(ind_clicked, [clclk.clicked_frames(ind_clicked).timestamp], int_ind_clicked, 'linear'); % lag time
dt = dt(1:size(curv_xcorr2, 2)) - dt(1);

% calculate lags in units of cycles and cilium length
ds_rel = ntt(1:size(curv_xcorr2,1)) ./ ntt(end); %relative
dt_rel = interp1(ind_clicked, [clclk.clicked_frames(ind_clicked).timestamp], int_ind_clicked, 'linear'); % lag time
dt_rel = dt_rel - dt_rel(1);
dt_rel = dt_rel ./ dt_rel(end) .* cycle_completion;
dt_rel = dt_rel(1:size(curv_xcorr2, 2));

% get cross section of autocorrelation down the line of slowest decay
[cdt_rel, cds_rel, curv_xcorr2_cross_section, curv_xcorr2_theta] = autocorr_cross_section(curv_xcorr2,...
    dt_rel(:), ds_rel(:), flag_debug);


% now look at average properties:

% calculate mean and std of curvature over time
mean_int_curv_1omum = nanmean(int_curv_mat_1omum,2);
std_int_curv_1omum = nanstd(int_curv_mat_1omum,[],2);
N_int_curv_1omum = sum(~isnan(int_curv_mat_1omum),2);

mean_int_curv_commonlength_1omum = nanmean(int_curv_mat_commonlength_1omum,2);
std_int_curv_commonlength_1omum = nanstd(int_curv_mat_commonlength_1omum,[],2);
N_int_curv_commonlength_1omum = sum(~isnan(int_curv_mat_commonlength_1omum),2);

%% return

% interpolated timepoints
clclk.Curv.int_ind_clicked      = int_ind_clicked;
clclk.Curv.int_timestamp      = [clclk.clicked_frames(int_ind_clicked).timestamp];

% interpolated cilia polynomial fitting
clclk.Curv.ntt          = ntt;          % arclength
clclk.Curv.ntt_cl       = ntt(1:size(int_curv_mat_commonlength,1));         % arclength
clclk.Curv.int_max_tt   = int_max_tt;   % where to stopalong the cilium for each timepoint
clclk.Curv.intPx        = intPx;        % interpolated polyfit
clclk.Curv.intPy        = intPy;

% curvature
clclk.Curv.int_curv_mat                     = int_curv_mat;                     % curvature(arclength, time)
clclk.Curv.int_curv_mat_1omum               = int_curv_mat_1omum;               % curvature in SI units
clclk.Curv.int_curv_mat_commonlength        = int_curv_mat_commonlength;        % curv choopped at min common length
clclk.Curv.int_curv_mat_commonlength_1omum  = int_curv_mat_commonlength_1omum;  % same but SI units
clclk.Curv.mean_int_curv_1omum              = mean_int_curv_1omum;              % curvature averaged over time, SI units
clclk.Curv.std_int_curv_1omum               = std_int_curv_1omum;               % std of curv over time, SI units
clclk.Curv.N_int_curv_1omum                 = N_int_curv_1omum;                 % non nans data points for std error
clclk.Curv.mean_int_curv_commonlength_1omum = mean_int_curv_commonlength_1omum; % curvature averaged over time, SI units
clclk.Curv.std_int_curv_commonlength_1omum  = std_int_curv_commonlength_1omum;  % std of curv over time, SI units
clclk.Curv.N_int_curv_commonlength_1omum    = N_int_curv_commonlength_1omum;    % non nans data points for std error

% curvature autocorrelation
clclk.Curv.ds               = ds;           % lag space
clclk.Curv.dt               = dt;           % lag time
clclk.Curv.curv_xcorr2      = curv_xcorr2;  % autocorrelation
clclk.Curv.ds_rel           = ds_rel;       % lag space, units of cilum_length
clclk.Curv.dt_rel           = dt_rel;       % lag time, units of cycle
clclk.Curv.cycle_completion = cycle_completion; % how much of a cycle is this

% autocorrelation cross section
clclk.Curv.cdt_rel                      = cdt_rel;                      % x coordinate of the cross section line
clclk.Curv.cds_rel                      = cds_rel;                      % y coordinate of the cross section line
clclk.Curv.curv_xcorr2_cross_section    = curv_xcorr2_cross_section;    % cross section at x,y
clclk.Curv.curv_xcorr2_theta            = curv_xcorr2_theta;            % angle of the cross section plane

clclk.old_savename = clclk.savename;
clclk.savename = clclk_filename;

if flag_save
    save([clclk.savename,'_curv'], '-struct', 'clclk');
end %if

if ~flag_plot
    return
end

%%%%%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%

%% figure, general properties

nplots = 2;

hf = figure(201);
clf
hf. Position = [120 300 560*nplots 420];

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
hs = surf(dt, ds, curv_xcorr2);
hs.EdgeColor = 'none';
xlabel(ha2, '\tau, [s]','FontSize',16)
ylabel(ha2, '\Deltas, [\mum]','FontSize',16)

hc2 = colorbar;
hc2.Label.String = 'C_\kappa(\Deltas,\tau)';
hc2.Label.FontSize = 16;

hf3 = figure;

ha3 = axes;

% imagesc(dt, ds, curv_xcorr2);
hs = surf(dt_rel, ds_rel, curv_xcorr2);
hs.EdgeColor = 'none';
xlabel(ha3, '\tau, [cycles]','FontSize',16)
ylabel(ha3, '\Deltas, [cilium length]','FontSize',16)

hc2 = colorbar;
hc2.Label.String = 'C_\kappa(\Deltas,\tau)';
hc2.Label.FontSize = 16;

% 
% hf3 = figure;
% ha3 = axes;
% hp3 = plot(ntt * clclk.px2mum, std_int_curv_1omum);
% xlabel(ha3,'Arclength, [\mum]','FontSize',16);
% ylabel(ha3,'\sigma(\kappa), [\mum^{-1}]','FontSize',16);
% 
% 
% hf4 = figure;
% ha4 = axes;
% hp4 = shadedErrorBar(ntt * clclk.px2mum, mean_int_curv_1omum, std_int_curv_1omum);
% xlabel(ha4,'Arclength, [\mum]','FontSize',16);
% ylabel(ha4,'<\kappa>, [\mum^{-1}]','FontSize',16);
%%
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



% try alternative method, finding the angle that maximises the integral of the curve
function [cx, cy, cross_section, theta_star] = autocorr_cross_section(IM, x, y, flag_debug)


if nargin < 4 || isempty(flag_debug)
    flag_debug = false;
else
    flag_debug = true;
end
% threshold for calculating the integral
r_star = 0.2; %this will be in either cilium length, cycle, or a combination of the 2

% rename fors quicker coding
% IM = clclk.Curv.curv_xcorr2;
% x = clclk.Curv.dt_rel(:);
% y = clclk.Curv.ds_rel(:);

% get initial estimate of propagation angle
nangles = 91;
theta_array = linspace(0, pi/2, nangles);
for i = numel(theta_array):-1:1
    theta = theta_array(i);
    cross_section_integral(i) = integral_of_autocorr_Xsection(IM, x, y, r_star, theta);
end %for
[~,wh] = max(cross_section_integral);
theta0 = theta_array(wh);

% create anonymous function to hide parameters and finnd maximising angle
% since using fminsearch, put a minus sign here
intAutocorrXsect = @(th)(-1)*integral_of_autocorr_Xsection(IM, x, y, r_star, th);
theta_star = fminsearch(intAutocorrXsect,theta0);

hf = figure;
clf
centred_resized_figure(hf,[1 2]);

ha1 = axes('Position',[0.1 0.15 0.35 0.8]);
ha2 = axes('Position',[0.6 0.15 0.35 0.8]);
ha2.NextPlot = 'add';
ha2.Box = 'on';

% plot autocorr first
% colormap(ha1,'gray');
surf(ha1,x, y, IM, 'edgecolor', 'none');
xlabel(ha1, '\tau, [cycles]','FontSize',16)
ylabel(ha1, '\Deltas, [cilium length]','FontSize',16)
% set(ha1,'ydir','normal');
ha1.NextPlot = 'add';

% get cross section
[cx,cy,cross_section] = improfile(x([1 end]), y([1 end]), IM,...
    [0, 1]*cos(theta_star), [0, 1]*sin(theta_star) , 'bilinear');

cx(isnan(cross_section)) = [];
cy(isnan(cross_section)) = [];
cross_section(isnan(cross_section)) = [];

% plot cross section on right axes
plot(ha2, hypot(cx,cy), cross_section,'.')

% plot the direction of maximum cross_section_integral on the
% autocorrelation map
plot3(ha1, cx, cy, cross_section, 'color', 'r');


if flag_debug
    
    nangles = 91;
    hf = figure;
    clf
    centred_resized_figure(hf,[2 2]);
    cmap = jet(nangles);
    
    ha1 = axes('Position',[0.1 0.15 0.35 0.35]);
    ha2 = axes('Position',[0.6 0.15 0.35 0.35]);
    ha2.NextPlot = 'add';
    ha2.Box = 'on';
    
    ha3 = axes('Position',[0.1 0.6 0.35 0.35]);
    ha3.NextPlot = 'add';
    ha3.Box = 'on';
    
    ha4 = axes('Position',[0.6 0.6 0.35 0.35]);
    

    % plot autocorr first
    % colormap(ha1,'gray');
    imagesc(ha1,x([1 end]), y([1 end]), IM);
    set(ha1,'ydir','normal');
    ha1.NextPlot = 'add';
    
    % plot decay curve at each angle
    theta_array = linspace(0, pi/2, nangles);
    for i = numel(theta_array):-1:1
        
        theta = theta_array(i);
        
        % get cross section
        [cx_,cy_,cross_section_] = improfile(x([1 end]), y([1 end]), IM,...
            [0, r_star]*cos(theta), [0, r_star]*sin(theta) , 'bilinear');
        
        % plot direction on left axes
        plot(ha1, [0, r_star]*cos(theta), [0, r_star]*sin(theta), 'color', cmap(i,:));
        
        % plot cross section on right axes
        plot(ha2, hypot(cx_,cy_), cross_section_, '.', 'color', cmap(i,:));
        
        % grow vector with value of integral of cross section
        cross_section_integral(i) = integral_of_autocorr_Xsection(IM, x, y, r_star, theta);
        
    end %for
    
    % plot cross_section_integral vs angle
    plot(ha3, theta_array, cross_section_integral)
    
    % then find the maximum
    [~,wh] = max(cross_section_integral);
    
    % plot the direction of maximum cross_section_integral on the
    % autocorrelation map
    imagesc(ha4,x([1 end]), y([1 end]), IM);
    set(ha4,'ydir','normal');
    ha4.NextPlot = 'add';
    plot(ha4, [0, r_star]*cos(theta_star), [0, r_star]*sin(theta_star), 'color', 'r');
    plot(ha4, [0, r_star]*cos(theta_array(wh)), [0, r_star]*sin(theta_array(wh)), 'color', 'k');
    
end %if

end


function [autocorr_Xsection_integral] = integral_of_autocorr_Xsection(IM, x, y, r_star, theta)
% measures the integral of the cross section at angle theta

theta = max(theta, 0);
theta = min(theta, pi/2);

% take Xsection at desired angle
[cx,cy,cross_section] = improfile(x([1 end]), y([1 end]), IM,...
    [0, r_star]*cos(theta), [0, r_star]*sin(theta) , 'bilinear');

% measure integral
autocorr_Xsection_integral = trapz(hypot(cx,cy), cross_section);

end %function



%% curvature matrix NOT PLOTTING THIS BECAUSE TT HAS GOT DFFERENT SPACING FOR EACH FRAME, SO CAN'T PUT IT AS A MATRIX
%
% axes(ha(3))
% ha(3).NextPlot = 'add';
% pbar = ha(3).PlotBoxAspectRatio;
%
% hs = pcolor(allcurvmat_XX,allcurvmat_YY,allcurvmat);
% hs.EdgeColor = 'none';



