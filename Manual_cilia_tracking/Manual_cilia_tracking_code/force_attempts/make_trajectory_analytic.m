clear all
close all

cd E:\Data\Cilia_Profile\Clicked_Cilia\N_clclk

fl = dir('*.clclk_force');

% fc = 4
fc = 7 % very good

% load file
clclk = load(fl(fc).name,'-mat');


%% change names of centre-of-drag coordinates for quicker coding. x == para, y == perp. In um

flag_cc = 1; %flag common_cylinders

if flag_cc
    cod_x = [clclk.Force.displ.commoncyls_cod_cpara]' .* clclk.Stroke.px2mum;
    cod_y = [clclk.Force.displ.commoncyls_cod_cperp]' .* clclk.Stroke.px2mum;
else
    cod_x = [clclk.Force.displ.cod_cpara]' .* clclk.Stroke.px2mum;
    cod_y = [clclk.Force.displ.cod_cperp]' .* clclk.Stroke.px2mum;
end %if
ts = [clclk.Force.displ.timestamp_s]';

% rescale ts so that it's between 0 and 1
T_s = (ts(end) + median(diff(ts)) - ts(1)); % "measured" period of a cilium as per clicking, in seconds
ts =(ts - ts(1)) ./ T_s;


% vectors of "velocities" of the centre of drag
% x velocity u, y velocity v (PIV notation)

cod_u = diff(cod_x) ./ diff(ts);
cod_v = diff(cod_y) ./ diff(ts);

% force
if flag_cc
    force_x = [clclk.Force.displ.commoncyls_F_para_pN]';
    force_y = [clclk.Force.displ.commoncyls_F_perp_pN]';
else
    force_x = [clclk.Force.displ.tot_F_para_pN]';
    force_y = [clclk.Force.displ.tot_F_perp_pN]';
end

% force angle (good to understand pow/rec)
force_angle = atan2(force_y, force_x);
    

%% plots

close all

% maxiplot
displ_plot(clclk,flag_cc)

% ------- angle plot -------%

figure(2); clf
hold on

% plot(ts,unwrap(force_angle,pi),'o')
ts_itp = interp(ts,4);                                                     % interpolated time
fa_itp = interp(unwrap(force_angle,pi),4);                                 % interpolated force_angle
sfx_itp = interp1( ts, sign(force_x) ,ts_itp);                             % interpolate sign of force at the interpolated times
comb_itp = abs(diff(sfx_itp)) .* abs(diff(fa_itp));                        % multiply |d/dt of sign| and |d/dt of angle| to highlight switching time

[~,ts_switch] = findpeaks(comb_itp, ts_itp(1:end-1));                      % find switching time

% plot(xx,yy,'.')
hold on
plot(ts(1:end-1),abs(diff(unwrap(force_angle,pi))),'o')                    % |d/dt of force_angle|
plot(ts_itp(1:end-1),abs(diff(fa_itp)),'.')                                % |d/dt interpolated force| this is nice

plot(ts_itp(1:end-1),abs(diff(sfx_itp)),'.')                               % |d/dt interpolated sign|

plot(ts_itp(1:end-1), comb_itp,'-k')                                       % multiplication |d/dt sign| and |d/dt angle|

plot(ts_switch*[1 1],get(gca,'YLim'),'g')                                  % switching time as found by findpeaks




%% maybe a good idea to take fft on x and y, and then only keep good modes?


% THIS IS VERY GOOD AND WORKS KEEP WORKING ON THIS



MT = 1; % this is supposed to be the period, but I have rescaled ts

% find centre
c_cod_x = mean(cod_x);
c_cod_y = mean(cod_y);


%---------------- fourier transform
ofxx = fft(cod_x,numel(cod_x));
ofyy = fft(cod_y,numel(cod_y));


% only keep a few modes
modes_to_keep = 4; % offset, 1, 2, 3
fxx = ofxx;
fyy = ofyy;
% fxx(modes_to_keep+1:end) = 0;
% fyy(modes_to_keep+1:end) = 0;
fxx(modes_to_keep+1:end-(modes_to_keep-1)) = 0;
fyy(modes_to_keep+1:end-(modes_to_keep-1)) = 0;

% delete first mode to centre the trajectory
fxx(1) = 0;
fyy(1) = 0;
                         


%-------------------- transform into analytical

% trajectories
fun_x_old = @(t) 2./length(cod_x) .* sum( arrayfun(@(m,tt) ...
    real(fxx(m)) * cos( 2*pi *(m-1) * tt/MT) - imag(fxx(m)) * sin( 2*pi *(m-1) * tt/MT),...
    repmat((1:modes_to_keep)',1,numel(t)),...
    repmat(t,modes_to_keep,1)));
fun_y_old = @(t) 2./length(cod_y) .* sum( arrayfun(@(m,tt) ...
    real(fyy(m)) * cos( 2*pi *(m-1) * tt/MT) - imag(fyy(m)) * sin( 2*pi *(m-1) * tt/MT),...
    repmat((1:modes_to_keep)',1,numel(t)),...
    repmat(t,modes_to_keep,1)));



% fun_y = @(t) 1./length(cod_y) .* (fyy(2) * sin(2*pi*t) + fyy(3) * sin(2*pi*2*t));

% create "analytical" time and angle
t_an = linspace(0,1,1e3);
phi_an = linspace(0,2*pi,1e3);


%------------------ find better centre
% as is, the centre is in the mean x and y position of fun_*(t_an). This
% though means that fun_phi(t_an) can be non-monotonic, and that is an issue
% when it comes to inverting it
[x0,y0] = find_new_trajectory_centre(fun_x_old(t_an), fun_y_old(t_an),1);
% x0=0;
% y0=0;

% now write this into fxx, fyy for ifft purposes
fxx(1) = -x0*length(cod_x);
fyy(1) = -y0*length(cod_y);

% inverse transform
x_plot = real(ifft(fxx,numel(cod_x),'symmetric'));
y_plot = real(ifft(fyy,numel(cod_y),'symmetric'));

% interpolation for plots, not really needed though
% ts_itp = interp(ts,4);                                                     % interpolated time
% xx_plot = interp1( ts, x_plot ,ts_itp);                            
% yy_plot = interp1( ts, y_plot ,ts_itp);  

% now write the better centre into fxx, fyy for analytical purposes (there's a factor of
% 2 )
fxx(1) = -x0*length(cod_x)/2;
fyy(1) = -y0*length(cod_y)/2;

% write the coefficient so that there is no need for an additional
% normalisation factor but the trajectories are just a sum of sines and
% cosines
x_traj_cos_coeff = 2 ./ length(cod_x) .*   real(fxx(1:modes_to_keep)) ;
x_traj_sin_coeff = 2 ./ length(cod_x) .* (-imag(fxx(1:modes_to_keep)));
y_traj_cos_coeff = 2 ./ length(cod_y) .*   real(fyy(1:modes_to_keep)) ;
y_traj_sin_coeff = 2 ./ length(cod_y) .* (-imag(fyy(1:modes_to_keep)));


% and redefine the trajectories
fun_x = @(t) sum( arrayfun(@(m,tt) ...
    x_traj_cos_coeff(m) * cos( 2*pi *(m-1) * tt/MT) + x_traj_sin_coeff(m) * sin( 2*pi *(m-1) * tt/MT),...
    repmat((1:modes_to_keep)',1,numel(t)),...
    repmat(t,modes_to_keep,1)));
fun_y = @(t) sum( arrayfun(@(m,tt) ...
    y_traj_cos_coeff(m) * cos( 2*pi *(m-1) * tt/MT) + y_traj_sin_coeff(m) * sin( 2*pi *(m-1) * tt/MT),...
    repmat((1:modes_to_keep)',1,numel(t)),...
    repmat(t,modes_to_keep,1)));


% velocities
% fun_u = @(t) diff(fun_x(t))./diff(t); % this doesn't work for single-point measurement
% fun_v = @(t) diff(fun_y(t))./diff(t);

fun_u = @(t) (fun_x(t + 1e-3) - fun_x(t)) .* 1e3;
fun_v = @(t) (fun_y(t + 1e-3) - fun_y(t)) .* 1e3;

fun_u_an = @(t) sum( arrayfun(@(m,tt) ...
    -(2*pi *(m-1)) * x_traj_cos_coeff(m) * sin( 2*pi *(m-1) * tt/MT) +...
     (2*pi *(m-1)) * x_traj_sin_coeff(m) * cos( 2*pi *(m-1) * tt/MT),...
    repmat((1:modes_to_keep)',1,numel(t)),...
    repmat(t,modes_to_keep,1)));
fun_v_an = @(t) sum( arrayfun(@(m,tt) ...
    -(2*pi *(m-1)) * y_traj_cos_coeff(m) * sin( 2*pi *(m-1) * tt/MT) +...
     (2*pi *(m-1)) * y_traj_sin_coeff(m) * cos( 2*pi *(m-1) * tt/MT),...
    repmat((1:modes_to_keep)',1,numel(t)),...
    repmat(t,modes_to_keep,1)));

% polar
fun_r = @(t) hypot(fun_y(t), fun_x(t));
fun_phi = @(t) unwrap( atan2(fun_y(t), fun_x(t)) );

% create t = f(phi) as lookup table 
% angle_for_lookuptable = mod(fun_phi(t_an(1:end-1)),2*pi);
size_lookup_old = 4096;
t_an2 = linspace(-1, 2, size_lookup_old+1);
angle_for_lookuptable_old = fun_phi(t_an2(1:end-1)) + 4*pi;
time_for_lookuptable_old = t_an2(1:end-1);

% sort so that angle is increasing (useful for interpolation, harmless here I think) 
[angle_for_lookuptable_old, sortingorder_old] = sort(angle_for_lookuptable_old);
time_for_lookuptable_old = time_for_lookuptable_old(sortingorder_old);


% THAT would work but also slow as lookup table for about three periods. We
% only need enough to have 0:2pi on theta OR 0:1 on tau
padding = 0.05;
idx_needed_range = isinrange(time_for_lookuptable_old, 0-padding, 1+padding) |...
    isinrange(angle_for_lookuptable_old, 0 - padding, 2*pi + padding);

% now find the first and last element of time within found range
time_extremes = minmax(time_for_lookuptable_old(idx_needed_range));

% and use it to create a second lookup table
size_lookup = 1024;
time_for_lookuptable = linspace(time_extremes(1), time_extremes(2), size_lookup+1);
time_for_lookuptable = time_for_lookuptable(1:end-1);
angle_for_lookuptable = fun_phi(time_for_lookuptable);

% now have to understand if angle_for_lookuptable is mostly negative or not 
while nnz(angle_for_lookuptable < 0) > nnz(angle_for_lookuptable > 0)
    angle_for_lookuptable = angle_for_lookuptable + 2*pi;
end

% sort so that angle is increasing (useful for interpolation, harmless here I think) 
[angle_for_lookuptable, sortingorder] = sort(angle_for_lookuptable);
time_for_lookuptable = time_for_lookuptable(sortingorder);

% now we are good, create a lookuptable for x and y (will speed up things)
x_for_lookuptable = fun_x(time_for_lookuptable);
y_for_lookuptable = fun_y(time_for_lookuptable);




% now this is good and can be written to file


absk = min(abs(angle_for_lookuptable));
% fun_t_of_phi = @(phi) interp1(angle_for_lookuptable, time_for_lookuptable, mod(phi+absk,2*pi)-2*pi-absk,'linear','extrap');
fun_t_of_phi = @(phi) interp1(angle_for_lookuptable, time_for_lookuptable, mod(phi,2*pi),'linear','extrap');
% fun_t_of_phi = @(phi) interp1(angle_for_lookuptable, time_for_lookuptable, mod(phi,2*pi),'spline','extrap');

% redefine experimental data based on latest translations, so that it is
% quicker to plot them

cod_xx = cod_x - c_cod_x - x0;
cod_yy = cod_y - c_cod_y - y0;
cod_rr = hypot(cod_xx,cod_yy);
cod_pphi = unwrap(atan2(cod_yy, cod_xx),pi/2);
cod_uu = diff(cod_xx)./diff(ts);
cod_vv = diff(cod_yy)./diff(ts);

%%
% save lookuptable and trajectory parameters in a file
[~,savename,~] = fileparts(fl(fc).name);
savename = strcat(savename, '.clclk_traj_par');
save(savename, 'x_traj_cos_coeff',...
    'x_traj_sin_coeff',...
    'y_traj_cos_coeff',...
    'y_traj_sin_coeff',...
    'angle_for_lookuptable',...
    'time_for_lookuptable');

%-------------------- write 2 files for the lookptable -------------%

% % names first
% [~,datsavename,~] = fileparts(fl(fc).name);
% angle_datsavename = strcat(datsavename, '_anglelookup.dat');
% time_datsavename  = strcat(datsavename, '_timelookup.dat');
% 
% % write angle
% fid = fopen(angle_datsavename, 'w');
% wc = fwrite(fid, angle_for_lookuptable, 'double');
% ws = fclose(fid);
% if ws == 0
%     disp('angles file closed ok')
%     disp([num2str(wc), ' entries written.']);
% else
%     disp('close failed')
%     disp([num2str(wc), ' entries written.']);
% end
% 
% 
% % write time
% fid = fopen(time_datsavename, 'w');
% wc = fwrite(fid, time_for_lookuptable, 'double');
% ws = fclose(fid);
% if ws == 0
%     disp('times file closed ok')
%     disp([num2str(wc), ' entries written.']);
% else
%     disp('close failed')
%     disp([num2str(wc), ' entries written.']);
% end

% that worked fine. Maybe it is better to write all parameters in a single
% file

% --------------------- write all parameters in a single file ------- %

% 1 double  => modes_to_keep, i.e. how many elements in *_traj*coeff
% 4 doubles => x_traj_cos_coeff
% 4 doubles => x_traj_sin_coeff
% 4 doubles => y_traj_cos_coeff
% 4 doubles => y_traj_sin_coeff
% 1 double  => N_lines of lookup table
% N_lines double => angle_for_lookuptable
% N_lines double => time_for_lookuptable
% N_lines double => x_for_lookuptable
% N_lines double => y_for_lookuptable

[~,datsavename,~] = fileparts(fl(fc).name);
trajparam_datsavename = strcat(datsavename, '_rotors_trajparam.dat');

% open
fid = fopen(trajparam_datsavename, 'w');

% write

% 1 double
wc = fwrite(fid, modes_to_keep, 'double');
if wc ~= 1, error('didn''t write all!'); end

% 4 x modes_to_keep doubles
wc = fwrite(fid, x_traj_cos_coeff, 'double');
if wc ~= modes_to_keep, error('didn''t write all!'); end
wc = fwrite(fid, x_traj_sin_coeff, 'double');
if wc ~= modes_to_keep, error('didn''t write all!'); end
wc = fwrite(fid, y_traj_cos_coeff, 'double');
if wc ~= modes_to_keep, error('didn''t write all!'); end
wc = fwrite(fid, y_traj_sin_coeff, 'double');
if wc ~= modes_to_keep, error('didn''t write all!'); end

wc = fwrite(fid, size_lookup, 'double');
if wc ~= 1, error('didn''t write all!'); end

wc = fwrite(fid, angle_for_lookuptable, 'double');
if wc ~= size_lookup, error('didn''t write all!'); end
wc = fwrite(fid, time_for_lookuptable, 'double');
if wc ~= size_lookup, error('didn''t write all!'); end

wc = fwrite(fid, x_for_lookuptable, 'double');
if wc ~= size_lookup, error('didn''t write all!'); end
wc = fwrite(fid, y_for_lookuptable, 'double');
if wc ~= size_lookup, error('didn''t write all!'); end

% close
ws = fclose(fid);

if ws == 0
    disp('param file closed ok');
end

return
%% --------------------------- plots

% close all



% -------------- positions
hf = figure(101); clf
hf.Units = 'centimeters';
hf.Position([3 4]) = [16 5];
hf.Units = 'normalized';
% centred_resized_figure(hf,[1 3]);

% trajectory
ha1 = axes('Position',[0.01+0.065 0.26 0.25 0.7]);
colormap(ha1,'jet');
hold on;
box on;
cmap = parula(numel(cod_xx));
scatter(cod_xx, cod_yy,20,cmap,'o','filled')
% plot(xx_plot,yy_plot,'.-')
cmapvel = jet(numel(fun_x(t_an)));
scatter(fun_x(t_an),fun_y(t_an),5,hypot(fun_u(t_an), fun_v(t_an)),'o','filled')
plot(0,0,'g*');
xlim(minmax(cod_xx)*1.2);
xlabel('x_{cod} - <x_{cod}>_t, [\mum]','FontSize',12);
ha1.XLabel.Units = 'normalized';
ha1.XLabel.Position(2) = -0.16;
ylabel('y_{cod} - <y_{cod}>_t, [\mum]','FontSize',12);
ha1.YLabel.Units = 'normalized';
ha1.YLabel.Position(1) = -0.1;
ha1.YLabel.Position(2) = 0.35;

% x
ha2 = axes('Position',[0.01+0.065*2+4/15 0.26 0.25 0.7]);
hold on;
box on;
plot(ts,cod_xx,'ob','MarkerSize',6)
plot(ts,x_plot,'.b','MarkerSize',6)
plot(t_an, fun_x(t_an),'b--')
plot(ts_itp, fun_x(ts_itp(:)'),'b-')
xlabel('t/T','FontSize',12);
ha2.XLabel.Units = 'normalized';
ha2.XLabel.Position(2) = -0.16;
ylabel('x_{cod} - <x_{cod}>_t, [\mum]','FontSize',12);
ha2.YLabel.Units = 'normalized';
ha2.YLabel.Position(1) = -0.1;
ha2.YLabel.Position(2) = 0.35;

% x spectrum
ha2in = axes('Position',[0.01+0.065*2+4/15 0.26 0.08 0.15]);
ha2in.FontSize = 6;
hold on;
box on;
plot(abs(ofxx(1:ceil(end/2))),'b--');
plot(abs(fxx(1:ceil(end/2))),'b-');
ha2in.XAxisLocation = 'top';
ha2in.YAxisLocation = 'right';
ylim(1.5*minmax(abs(fxx(1:ceil(end/2)))))
xlim([1,numel(fxx(1:ceil(end/2)))])
xlabel('modes')
ha2in.XLabel.Units = 'normalized';
ha2in.XLabel.Position(2) = 1.35;


% y
ha3 = axes('Position',[0.01+0.065*3+8/15 0.26 0.25 0.7]);
hold on
box on
plot(ts,cod_yy,'or','MarkerSize',6)
plot(ts,y_plot,'.r','MarkerSize',6)
plot(t_an, fun_y(t_an),'r--')
plot(ts_itp, fun_y(ts_itp(:)'),'r-')
xlabel('t/T','FontSize',12);
ha3.XLabel.Units = 'normalized';
ha3.XLabel.Position(2) = -0.16;
ylabel('y_{cod} - <y_{cod}>_t, [\mum]','FontSize',12);
ha3.YLabel.Units = 'normalized';
ha3.YLabel.Position(1) = -0.1;
ha3.YLabel.Position(2) = 0.35;

% y spectrum
ha3in = axes('Position',[0.01+0.065*3+8/15 0.26 0.08 0.15]);
ha3in.FontSize = 6;
hold on;
box on;
plot(abs(ofyy(1:ceil(end/2))),'r--');
plot(abs(fyy(1:ceil(end/2))),'r-');
ha3in.XAxisLocation = 'top';
ha3in.YAxisLocation = 'right';
ylim(1.5*minmax(abs(fyy(1:ceil(end/2)))))
xlim([1,numel(fyy(1:ceil(end/2)))])
xlabel('modes')
ha3in.XLabel.Units = 'normalized';
ha3in.XLabel.Position(2) = 1.35;

% print2svg(hf,'D:\Dropbox\Universita\PhD\Cilia_to_rotors\Figures\xy_trajectory')


%% ------------- velocities

hf = figure(102); clf
hf.Units = 'centimeters';
hf.Position([3 4]) = [16 5];
hf.Units = 'normalized';
% centred_resized_figure(hf,[1 3]);

%v(u)
ha1 = axes('Position',[0.01+0.065 0.26 0.25 0.7]);
hold on;
box on;
cmap = parula(numel(cod_uu));
scatter(cod_uu, cod_vv,20,cmap,'o','filled')
plot(fun_u(t_an), fun_v(t_an))
xlabel('u_{cod}, [\mum]','FontSize',12);
ylabel('v_{cod}, [\mum]','FontSize',12);
ha1.YLabel.Units = 'normalized';
ha1.YLabel.Position(1) = -0.12;

% u(t)
ha2 = axes('Position',[0.01+0.065*2+4/15 0.26 0.25 0.7]);
hold on
box on
plot(ts(1:end-1),cod_uu,'bo','MarkerSize',6);
plot(t_an, fun_u(t_an),'b')
xlabel('t/T','FontSize',12);
ylabel('u_{cod}, [\mum]','FontSize',12);
ha2.YLabel.Units = 'normalized';
ha2.YLabel.Position(1) = -0.12;

%v(t)
ha3 = axes('Position',[0.01+0.065*3+8/15 0.26 0.25 0.7]);
hold on
box on
plot(ts(1:end-1),cod_vv,'ro','MarkerSize',6);
% plot(t_an(1:end-1), diff(fun_y(t_an))./diff(t_an),'r')
plot(t_an, fun_v(t_an),'r')
xlabel('t/T','FontSize',12);
ylabel('v_{cod}, [\mum]','FontSize',12);
ha3.YLabel.Units = 'normalized';
ha3.YLabel.Position(1) = -0.12;

% print2svg(hf,'D:\Dropbox\Universita\PhD\Cilia_to_rotors\Figures\uv_velocity')


%% -------------- polar coordinates

hf = figure(103); clf
hf.Units = 'centimeters';
hf.Position([3 4]) = [16 5];
hf.Units = 'normalized';

% coordinates
cmap = parula(numel(cod_xx));
polarscatter(cod_pphi, cod_rr,20,cmap,'o','filled')
ha1 = gca;
ha1.Position = [0.06 0.15 0.25 0.7];
hold on
box on
colormap(ha1,'jet');
polarscatter(fun_phi(t_an), fun_r(t_an),...
    5,hypot(fun_u(t_an), fun_v(t_an)),'o','filled')


% r
ha2 = axes('Position',[0.01+0.065*2+4/15 0.26 0.25 0.7]);
hold on
box on
plot(ts, cod_rr, 'bo','MarkerSize',6);
plot(t_an, fun_r(t_an),'b')
xlabel('t/T','FontSize',12);
ylabel('r_{cod} - <r_{cod}>_t, [\mum]','FontSize',12);
ha2.YLabel.Units = 'normalized';
ha2.YLabel.Position(1) = -0.1;
ha2.YLabel.Position(2) = 0.35;

% phi
ha3 = axes('Position',[0.01+0.065*3+8/15 0.26 0.25 0.7]);
hold on
box on
plot(ts, cod_pphi, 'ro','MarkerSize',6);
plot(t_an, fun_phi(t_an),'r')
xlabel('t/T','FontSize',12);
ylabel('\phi, [rad]','FontSize',12);
ylim([-2.5*pi,pi/2]);
ha3.YLabel.Units = 'normalized';
ha3.YLabel.Position(1) = -0.18;
ha3.YLabel.Position(2) = 0.35;
setpiyticklabels(ha3,1)


% print2svg(hf,'D:\Dropbox\Universita\PhD\Cilia_to_rotors\Figures\polar_trajectory')


%% coordinates vs angle

% t_an2 = linspace(0,2,2e2);
t_an2 = t_an;
% xX = [-0.5 0.5] + minmax(fun_phi(t_an2));
xX = [-0.1 0.1] + [0 2*pi];
phi_an2 = linspace(0,2*pi,50);


hf = figure(104); clf;
hf.Units = 'centimeters';
hf.Position([3 4]) = [16 5];

ha1 = axes('Position',[0.005+0.065 0.26 0.25 0.7]);
hold on
box on
scatter(mod(fun_phi(t_an2),2*pi), fun_x(t_an2),20,t_an2,'filled');
scatter(phi_an2, fun_x(fun_t_of_phi(phi_an2)),10,'rx');
plot(mod(cod_pphi,2*pi), cod_xx, 'bo','MarkerSize',6);
xlim(xX)
xlabel('\phi, [rad]','FontSize',12);
ylabel('x_{cod}, [\mum]','FontSize',12);
ha1.YLabel.Units = 'normalized';
ha1.YLabel.Position(1) = -0.1;
setpixticklabels(ha1)

ha2 = axes('Position',[0.01+0.065*2+4/15 0.26 0.25 0.7]);
hold on
box on
scatter(mod(fun_phi(t_an2),2*pi), fun_y(t_an2),20,t_an2,'filled');
scatter(phi_an2, fun_y(fun_t_of_phi(phi_an2)),10,'rx');
plot(mod(cod_pphi,2*pi), cod_yy, 'ro','MarkerSize',6);
xlim(xX)
xlabel('\phi, [rad]','FontSize',12);
ylabel('y_{cod}, [\mum]','FontSize',12);
ha2.YLabel.Units = 'normalized';
ha2.YLabel.Position(1) = -0.12;
setpixticklabels(ha2)

ha3 = axes('Position',[0.01+0.065*3+8/15 0.26 0.25 0.7]);
hold on
box on
scatter(mod(fun_phi(t_an2),2*pi), fun_r(t_an2),20,t_an2,'filled');
scatter(phi_an2, fun_r(fun_t_of_phi(phi_an2)),10,'rx');
plot(mod(cod_pphi,2*pi), cod_rr,...
    'o','MarkerSize',6);
xlim(xX)
xlabel('\phi, [rad]','FontSize',12);
ylabel('r_{cod}, [\mum]','FontSize',12);
ha3.YLabel.Units = 'normalized';
ha3.YLabel.Position(1) = -0.1;
setpixticklabels(ha3)

% print2svg(hf,'D:\Dropbox\Universita\PhD\Cilia_to_rotors\Figures\xy_trajectory_vs_angle.svg');



%% velocities vs angle

% t_an2 = linspace(0,2,2e4);
t_an2 = t_an;
phi_an2 = linspace(0,2*pi,50);

% xX = [-0.5 0.5] + minmax(fun_phi(t_an2));
xX = [-0.1 0.1] + [0 2*pi];


hf = figure(105); clf;
hf.Units = 'centimeters';
hf.Position([3 4]) = [16 5];

ha1 = axes('Position',[0.015+0.065 0.26 0.25 0.7]);
hold on
box on
scatter(mod(fun_phi(t_an2),2*pi), fun_u(t_an2),20,t_an2,'filled');
scatter(phi_an2, fun_u(fun_t_of_phi(phi_an2)),10,'rx');
% plot(cod_pphi(1:end-1),cod_uu,'bo','MarkerSize',6);
xlim(xX)
xlabel('\phi, [rad]','FontSize',12);
ylabel('u_{cod}, [\mum]','FontSize',12);
ha1.YLabel.Units = 'normalized';
ha1.YLabel.Position(1) = -0.15;
ha1.YLabel.Position(2) = 0.4;
setpixticklabels(ha1)


ha2 = axes('Position',[0.005+0.065*2+4/15 0.26 0.25 0.7]);
hold on
box on
scatter(mod(fun_phi(t_an2),2*pi), fun_v(t_an2),20,t_an2,'filled');
scatter(phi_an2, fun_v(fun_t_of_phi(phi_an2)),10,'rx');
% plot(cod_pphi(1:end-1),cod_vv,'bo','MarkerSize',6);
xlim(xX)
xlabel('\phi, [rad]','FontSize',12);
ylabel('v_{cod}, [\mum]','FontSize',12);
ha2.YLabel.Units = 'normalized';
ha2.YLabel.Position(1) = -0.08;
ha2.YLabel.Position(2) = 0.4;
setpixticklabels(ha2)

ha3 = axes('Position',[0.01+0.065*3+8/15 0.26 0.25 0.7]);
hold on
box on
scatter(mod(fun_phi(t_an2),2*pi), hypot(fun_u(t_an2), fun_v(t_an2)),20,t_an2,'filled');
scatter(phi_an2, hypot(fun_u(fun_t_of_phi(phi_an2)),fun_v(fun_t_of_phi(phi_an2))),10,'rx');
% plot(cod_pphi(1:end-1),hypot(cod_uu,cod_vv),'o','MarkerSize',6);
xlim(xX)
xlabel('\phi, [rad]','FontSize',12);
ylabel('(u^2+v^2)^{1/2}_{cod}, [\mum]','FontSize',12);
ha3.YLabel.Units = 'normalized';
ha3.YLabel.Position(1) = -0.12;
ha3.YLabel.Position(2) = 0.4;
setpixticklabels(ha3)

% print2svg(hf,'D:\Dropbox\Universita\PhD\Cilia_to_rotors\Figures\uv_velocity_vs_angle.svg');


%% time vs angle


% t_an2 = linspace(0,2,2e4);
t_an2 = t_an;
phi_an2 = linspace(0, 4*pi, 1e3);
phi_an2(end) = [];




hf = figure(106); clf;
hf.Units = 'centimeters';
hf.Position([3 4]) = [16 5];

ha1 = axes('Position',[0.015+0.065 0.26 0.25 0.7]);
hold on
box on
plot(t_an2,fun_phi(t_an2) ,'LineWidth',1.2)
plot(ts, cod_pphi, 'o','MarkerSize',6)
ylim([-0.5 0.5] + minmax(fun_phi(t_an2)));
ylabel('\phi_{cod}, [rad]','FontSize',12);
xlabel('t/T','FontSize',12);
ha1.YLabel.Units = 'normalized';
ha1.YLabel.Position([1 2]) = [-0.13 0.45];
yticks([-2*pi:pi:0]);
setpiyticklabels(ha1)

ha2 = axes('Position',[0.005+0.065*2+4/15 0.26 0.25 0.7]);
hold on
box on
plot(angle_for_lookuptable, time_for_lookuptable);
xlim([-2*pi-0.5 4*pi+0.5]);
xticks([-2*pi:2*pi:4*pi]);
xlabel('angle lookuptable','FontSize',12);
ylabel('time lookuptable','FontSize',12);
ha2.YLabel.Units = 'normalized';
ha2.YLabel.Position([1 2]) = [-0.1 0.45];
setpixticklabels(ha2)


ha3 = axes('Position',[0.01+0.065*3+8/15 0.26 0.25 0.7]);
hold on
box on
plot(phi_an2, fun_t_of_phi(phi_an2),'.','MarkerSize',1)
% plot(phi_an2(1:30:end), fun_t_of_phi(phi_an2(1:30:end)), 'r.')
xlim([-0.5 0.5] + minmax(phi_an2) );
ylim([-0.05 0.05] + minmax(fun_t_of_phi(phi_an2)));
xlabel('\phi, [rad]','FontSize',12);
ylabel('t/T = f(\phi)','FontSize',12);
ha3.YLabel.Units = 'normalized';
ha3.YLabel.Position([1 2]) = [-0.18 0.45];
xticks(0:pi:4*pi);
setpixticklabels(ha3)

% print2svg(hf,'D:\Dropbox\Universita\PhD\Cilia_to_rotors\Figures\angle_lookup.svg');

%% compare performances in vector lookup vs function lookup

tic; fun_t_of_phi(phi_an);toc
tic; for i=1:numel(phi_an), fun_t_of_phi(phi_an(i));end; toc




%% animated plot showing local "tangent" velocity

hf = figure(107);
clf;
axes;

scatter(fun_x(fun_t_of_phi(phi_an)),fun_y(fun_t_of_phi(phi_an)),2,...
    hypot(fun_u(fun_t_of_phi(phi_an)),fun_v(fun_t_of_phi(phi_an))));
hold on;
axis image
xlim([-0.4 0.4] + minmax(fun_x(t_an)));
ylim([-0.4 0.4] + minmax(fun_y(t_an)));


for phi = phi_an
    
    tt = fun_t_of_phi(phi);
    hq = quiver(fun_x(tt), fun_y(tt), 0.05*fun_u(tt), 0.05*fun_v(tt),'bo');
    hp = plot([0 fun_x(tt)],[0 fun_y(tt)],'r');
    drawnow;
    pause(0.05)
    delete(hq);
    delete(hp);
    
end

%% 
figure(108);
clf
% hold on

N_angles = 60;
phi_an2 = linspace(0,2*pi-2*pi/N_angles,N_angles);

polarplot(fun_phi(fun_t_of_phi(phi_an2)),fun_r(fun_t_of_phi(phi_an2)),'.')
hold on


for i = phi_an2
polarplot(fun_phi(fun_t_of_phi(i)).*[1 1],fun_r(fun_t_of_phi(i)).*[0 1],'b');
end





