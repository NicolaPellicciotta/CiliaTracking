function [ hf, fmat, vmat ] = plot_force_along_cilium( clclk, flag_plotlog )
%plot_force_along_cilium plots force along cilium as calculated by
%clclk_calculate_force

%% input check (copied from displ_plot)

if nargin < 2 || isempty(flag_plotlog)
    flag_plotlog = false;
end

if isfield(clclk, 'displ') % then input argument is actually clclk.Force
    dummy = clclk;
    clear clclk;
    clclk.Force = dummy;
    clear dummy;
end


if isfield(clclk,'Force') % then function invoked as intended
    
    displ = clclk.Force.displ;
    
    if isfield(clclk.Force,'baseline')
        baseline = clclk.Force.baseline;
    else
        baseline.parav_x = 1;
        baseline.parav_y = 0;
        baseline.perpv_x = 0;
        baseline.perpv_y = 1;
    end %if
    
end %if


% flag for displ being checked against experimental error in clicking
% length
% flag is true by default
% if nargin < 2 || isempty(flag_commoncyls)
%     flag_commoncyls = true;
% end %if

% but if flag is true and no common cylinders, prevent error and warn user
% if flag_commoncyls
%     if isfield(displ,'commoncyls_ll_um')
%         flag_commoncyls = true;
%     else
%         warning('Commoncyls properties not present.');
%         flag_commoncyls = false;
%     end
% end
% keyboard

%% figure


hf = figure;
hf.Color = 'w';

hapos = [0.1 0.15 0.75 0.8];

haf = axes;
haf.Position = hapos;
haf.Box = 'on';
haf.NextPlot = 'add';

% at each click
if isvector(displ(1).cil(1).tr_um)
    disp('legacy')
    for i = 1:numel(displ)
        
        % find the approximate arclength of the centre of the cylinders
        cl = cumsum([0;displ(i).cyl_al_um] );
        
        % plot it "straightened" at the right time
        
        hp(i) = plotc(displ(i).timestamp_s .* ones(size(cl)),...
            cl, ...
            [displ(i).cyl_Fampl_pN./displ(i).cyl_al_um;0], 8, 1);
    end %for
else %new version
    
    % define common arclength for all plotting
    % I know that cilia are at most long 12um, let's use 15um as an upper
    % limit just in case
    cl = linspace(0,15,150); %in um
    cl = cl(:);
    
    % initialise matrix to store results
    % find how many frames is the difference between the first and last
    % clicked cilium
    fmat = nan(numel(cl), range(clclk.idx_clicked_frames)+1); 
    vmat = nan(numel(cl), range(clclk.idx_clicked_frames)+1); 
    
    for i = 1:numel(displ)
        
%         % define an arclength for plotting purposes
%         cl = linspace(min(displ(i).cil(1).tr_um(:)),max(displ(i).cil(1).tr_um(:)),100);
%         cl = cl(:);
        
        % now find which cilinders span the point identified by cl
        idx = cl >= displ(i).cil(1).tr_um(:,1)' & cl <= displ(i).cil(1).tr_um(:,2)';
        ff = repmat(displ(i).cyl_Fampl_pN(:)'./displ(i).cyl_al_um(:)', numel(cl), 1 ); %force per unit length
        ff(~idx) = nan;
        F = nanmean(ff,2);
        
        displ(i).cyl_vel_umos = hypot(displ(i).cyl_displ_x_um, displ(i).cyl_displ_y_um)./...
            displ(i).time_interval_s;
        vv = repmat(displ(i).cyl_vel_umos(:)', numel(cl), 1 ); %velocity
        vv(~idx) = nan;
        V = nanmean(vv,2);
        
        if flag_plotlog
            Fcolors = log10(F);
        else
            Fcolors = F;
        end
        
        hp(i) = plotc(displ(i).timestamp_s .* ones(size(cl)),...
            cl, Fcolors, 8, 1);
        
        % write in matrix
        ind_in_mat = ...
            clclk.idx_clicked_frames( clclk.Force.displ(i).cil(1).timestamp==[clclk.clicked_frames(clclk.idx_clicked_frames).timestamp]);
        ind_in_mat = ind_in_mat - min(clclk.idx_clicked_frames) + 1;
        fmat(:,ind_in_mat) = F;
        fmat(end, ind_in_mat) = 1; %this says that this entry was real, not interpolated!
        
        vmat(:, ind_in_mat) = V;
        
    end %for
    
    
    
end %if


haf.XLabel.String = 'Time, [s]';
haf.XLabel.FontSize = 14;
haf.YLabel.String = 'Arclength, [\mum]';
haf.YLabel.FontSize = 14;


hc = colorbar;
haf.Position = hapos;

hc.Units = 'Normalized';
hc.Position = [sum(hapos([1 3])) + 0.01, hapos(2), 0.03,  hapos(4)];
if flag_plotlog
hc.Label.String = 'log_{10}(Force per unit length / (pN/\mum))';
else
hc.Label.String = 'Force per unit length, [pN/\mum]';
end
hc.Label.FontSize = 14;


figure
imagesc(fmat)

figure
imagesc(vmat)

end

