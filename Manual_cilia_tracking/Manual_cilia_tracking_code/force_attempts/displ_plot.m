function [ hf, hf2, hf3 ] = displ_plot( clclk, flag_commoncyls, flag_cylforces )
%displ_plot Plots the force of the cilium on the fluid


%% input check


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
if nargin < 2 || isempty(flag_commoncyls)
    flag_commoncyls = true;
end %if

if nargin < 3 || isempty(flag_cylforces)
    flag_cylforces = false;
end

% but if flag is true and no common cylinders, prevent error and warn user
if flag_commoncyls
    if isfield(displ,'commoncyls_ll_um')
        flag_commoncyls = true;
    else
        warning('Commoncyls properties not present.');
        flag_commoncyls = false;
    end
end




%% figure

hf = figure;
hf.Position(4) = hf.Position(4) * 1.5;
hf.Position(2) = hf.Position(2) * 0.5;
hf.Position(3) = hf.Position(3) * 1.5;
hf.Color = 'w';

%----------------------------------- axes with forces ---------------------------------%

haf = axes;
haf.Position = [0.1 0.08 0.57 0.35];
hold on
haf.Box = 'on';
haf.XGrid = 'on';
haf.GridAlpha = 0.08;

% plot force over time
hp_par = plot([displ.timestamp_s], [displ.tot_F_para_pN],'g.--','LineWidth',1.2);
hp_per = plot([displ.timestamp_s], [displ.tot_F_perp_pN],'r.--','LineWidth',1.2);
hp_tot = plot([displ.timestamp_s], hypot([displ.tot_F_para_pN],[displ.tot_F_perp_pN]),'k.-','LineWidth',1.2);

% fix axes
haf.XLim = minmax([displ.timestamp_s]) + mean(diff([displ.timestamp_s])) .* [-1.5, +1.5];
haf.YLim = haf.YLim + [0 mean(diff(haf.YTick))/3];

% legend
hleg = legend([hp_par, hp_per, hp_tot],'F_{//}','F_{\perp}','F_{tot}');
hleg.Orientation = 'Horizontal';
hleg.Location = 'SouthEast';
hleg.Box = 'off';
hleg.FontSize = 12;
hleg.Position(2) = haf.Position(2) + 0.01 ;

% labels
xlabel('Time, [s]','FontSize',14);
ylabel('Force, [pN]','FontSize',14);

%------------------------------------ axes with force phase ----------------------------%

haph = axes;
haph.Position = haf.Position;
haph.Position(2) = sum(haph.Position([2 4]));
haph.Position(4) = 0.2;
haph.Box = 'on';
hold on;
haph.XGrid = 'on';
haph.YGrid = 'on';
haph.GridAlpha = 0.08;

% plot phase of the force
hp_ph = plot([displ.timestamp_s], atan2([displ.tot_F_perp_pN],[displ.tot_F_para_pN]),...
    'o','MarkerSize',6,'Color',[0.8500 0.3250 0.098],'LineWidth',1.2);

% and the one on the adjusted number of cylinders
if flag_commoncyls
    hp_ph_cc = plot([displ.timestamp_s], atan2([displ.commoncyls_F_perp_pN],[displ.commoncyls_F_para_pN]),...
        'o','MarkerSize',6,'Color',[0 0.447 0.741],'LineWidth',1.2);
end %if

% fix x axis limit
haph.XLim = haf.XLim;
haph.XTickLabel = [];

% put y axis in units of pi
haph.YLim = [-pi-0.2,pi+0.2];
haph.YTick = -pi:pi/2:pi;
haph.YTickLabels = {'-\pi'; '-\pi/2'; '0'; '\pi/2'; '\pi'};


% labels
ylabel('Angle \phi, [rad]','FontSize',14);


%------------------------------------- axes with cilium length -------------------------%

hal = axes;
hal.Position = haph.Position;
hal.Position(2) = sum(hal.Position([2 4]));
hal.Position(4) = 0.2;
hal.Box = 'on';
hal.XGrid = 'on';
hal.GridAlpha = 0.08;


% plot cilium length as used to calculate the force
yyaxis(hal,'left')
hold on
hp_tcl = plot([displ.timestamp_s],[displ.tot_ll_um],'o','MarkerSize',6,'LineWidth',1.2);
hp_mcl = plot([displ.timestamp_s],mean([displ.tot_ll_um]).*ones(size([displ.timestamp_s])),'-k','LineWidth',1.2);

% fix axes
hal.XLim = haf.XLim;
hal.XTickLabels = [];
hal.YLim = minmax([displ.tot_ll_um]) + max(diff([displ.tot_ll_um])).*[-1 1];
hal.YMinorTick = 'on';
hal.YColor = 'k';

% labels
ylabel('\Sigma_i{\deltal_i}, [\mum]','FontSize',14)

% now put relative scale on the other axis
yyaxis(hal,'right');
plot([displ.timestamp_s],[displ.tot_ll_um]./mean([displ.tot_ll_um]),'-','LineWidth',1.2);
hal.YLim = (minmax([displ.tot_ll_um]) + max(diff([displ.tot_ll_um])).*[-1 1])./mean([displ.tot_ll_um]);


% now, if I have also done the final check
if flag_commoncyls
    
    yyaxis 'left'
    hold on
    
    hp_ccl = plot([displ.timestamp_s],[displ.commoncyls_ll_um],'o-','MarkerSize',6,'LineWidth',1.2);
    
end %if

%----------------------------- axes with cilium shape -----------------------------%

has = axes;
has.Position = hal.Position;
has.Position(2) = sum(hal.Position([2 4]));
has.Position(4) = 0.15;

% apply transformation to cilia
for i = numel(displ):-1:1
    [plstr(i).xpl1, plstr(i).ypl1] = xy2paraperp(baseline, displ(i).cil(1).xx, displ(i).cil(1).yy);
    [plstr(i).xpl2, plstr(i).ypl2] = xy2paraperp(baseline, displ(i).cil(2).xx, displ(i).cil(2).yy);
    plstr(i).meanxpl = mean([plstr(i).xpl1(1), plstr(i).xpl2(1)]);
end

% find how many ums in height we need
minmax_ypl = minmax(vertcat(plstr.ypl1));

% set ylimits accordingly
has.YLim = minmax_ypl + [-2 +2];

% and set xlimits to keep axis "image" like
% has.XLim =
xspan = (hal.Position(3)*hf.Position(3))/(has.Position(4)*hf.Position(4)) * diff(has.YLim);
xlimits = [0, xspan];
has.XLim = xlimits;

% find where to place each base of the cilium to make it correspond to the
% timepoints
xoffsets = ([displ.timestamp_s] - hal.XLim(1))./diff(hal.XLim) .* xspan;

% plot data allowing for xoffsets
hold on
for i = 1:numel(displ)
    
    plstr(i).oxpl2 = plstr(i).xpl2 + xoffsets(i) - plstr(i).meanxpl;
    plstr(i).oxpl1 = plstr(i).xpl1 + xoffsets(i) - plstr(i).meanxpl;
    
    plot(plstr(i).oxpl1, plstr(i).ypl1, 'color', 0.2*[1 1 1],'Clipping','off');
    plot(plstr(i).oxpl2, plstr(i).ypl2,'k','Clipping','off');
    
end

% remove useless ticks
has.YTick = [];
has.XTick = [];
has.XColor = 'none';
has.YColor = 'none';
has.Color = 'none';

% fix data ratio
has.PlotBoxAspectRatioMode = 'manual';
has.DataAspectRatioMode = 'manual';
has.DataAspectRatio = [1 1 1];


%-------------------- axes with force calculated on common clinders -----------------%


if flag_commoncyls
    
    % halve the size of the first axes
    haf.Position(4) = haf.Position(4)/2;
    
    hacf = axes;
    hacf.Position = haf.Position;
    hacf.Position(2) = sum(hacf.Position([2 4]));
    
    hold on
    box on
    
    hacf.XGrid = 'on';
    hacf.GridAlpha = 0.08;
    
    
    % plot common cylinders' forces
    hpc_par = plot([displ.timestamp_s],[displ.commoncyls_F_para_pN],'g.--','LineWidth',1.2);
    hpc_per = plot([displ.timestamp_s],[displ.commoncyls_F_perp_pN],'r.--','LineWidth',1.2);
    hpc_tot = plot([displ.timestamp_s],hypot([displ.commoncyls_F_para_pN], [displ.commoncyls_F_perp_pN]),'k.-','LineWidth',1.2);
    
    
    % legend
    hleg = legend([hpc_par, hpc_per, hpc_tot],'F_{//}^{cc}','F_{\perp}^{cc}','F_{tot}^{cc}');
    hleg.Orientation = 'Horizontal';
    hleg.Location = 'SouthEast';
    hleg.Box = 'off';
    hleg.FontSize = 12;
    hleg.Position(2) = hacf.Position(2) + 0.01 ;
    
    % labels
    ylabel('Force, [pN]','FontSize',14);
    
    % fix axes limits
    hacf.YLim = haf.YLim;
    hacf.XLim = haf.XLim;
    
    % fix axes
    hacf.XTickLabels = [];
    
end %if





%--------------------- draw cylinders on cilium shape ---------------------%

% axes in focus
axes(has);

% apply transformation matrix
for i = 1:numel(displ)
    
    % apply transformation matrix
    [plstr(i).xplr1, plstr(i).yplr1] = xy2paraperp(baseline, displ(i).cil(1).xr, displ(i).cil(1).yr);
    [plstr(i).xplr2, plstr(i).yplr2] = xy2paraperp(baseline, displ(i).cil(2).xr, displ(i).cil(2).yr);
    
    % apply offset
    plstr(i).oxplr2 = plstr(i).xplr2 + xoffsets(i) - plstr(i).meanxpl;
    plstr(i).oxplr1 = plstr(i).xplr1 + xoffsets(i) - plstr(i).meanxpl;
    
    % now plot in 2 colours the cylinders used and not used to calculate
    % the force
    
    if flag_commoncyls
        
        idx_plot = vertcat(true, displ(i).idx_commoncyls);
        
        % cylinders used
        plot(plstr(i).oxplr1(idx_plot), plstr(i).yplr1(idx_plot), 'color', 'g', 'Clipping','off');
        plot(plstr(i).oxplr2(idx_plot), plstr(i).yplr2(idx_plot), 'color', 'g', 'Clipping','off');
        
        
        
    else
        
        plot(plstr(i).oxplr1, plstr(i).yplr1, 'color', [0.85 0.325 0.098],'Clipping','off');
        plot(plstr(i).oxplr2, plstr(i).yplr2, 'color', [0.85 0.325 0.098],'Clipping','off');
        
    end %if
    
end %for





%------------- Axes with cilium shape over time, and tiny forces -------------%

% axes positions
hatf = axes;
hatf.Position(2) = haf.Position(2);
hatf.Position(1) = sum(haf.Position([1 3])) + 0.02;
hatf.Position(3) = 1-hatf.Position(1);
if flag_commoncyls
    hatf.Position(4) = haf.Position(4) + hacf.Position(4);
else
    hatf.Position(4) = haf.Position(4);
end
hatf.Box = 'on';

hatf.XTick = [];
hatf.YTick = [];
hatf.Visible = 'off';

hold on

% for loop on plotting forces
for i = 1:numel(displ)
    
    % first plot cilia
    plot(plstr(i).xpl1, plstr(i).ypl1, 'color',[0.85 0.325 0.098], 'Clipping','off');
    
    % special case: last position
    if i == numel(displ)
        plot(plstr(i).xpl2, plstr(i).ypl2, 'color',[0.85 0.325 0.098], 'Clipping','off');
    end %if
    
    
    %then plot forces of each cylinder
    [plstr(i).cxpl, plstr(i).cypl] = xy2paraperp(baseline, displ(i).cil(1).cyl_cx, displ(i).cil(1).cyl_cy);
    [plstr(i).fxpl, plstr(i).fypl] = xy2paraperp(baseline, displ(i).cyl_F_x_pN, displ(i).cyl_F_y_pN);
    
    if flag_commoncyls
        quiver(plstr(i).cxpl(displ(i).idx_commoncyls), plstr(i).cypl(displ(i).idx_commoncyls),...
            plstr(i).fxpl(displ(i).idx_commoncyls), plstr(i).fypl(displ(i).idx_commoncyls), 'b' );
    else
        quiver(plstr(i).cxpl, plstr(i).cypl,...
            plstr(i).fxpl, plstr(i).fypl, 'b' );
    end
    
end %for

axis image





%------------------ Axes with trajectory of centre of drag ---------------------%

% axes position
hacod = axes;
hacod.Position(3) = hatf.Position(3) .* 0.79;
hacod.Position(1) = hatf.Position(1) + 0.05*hatf.Position(3);
hacod.Position(2) = sum(hatf.Position([2 4])) + 0.05;
hacod.Position(4) = hatf.Position(4);


hacod.Box = 'on';
hold on

% center the CoD positions
if flag_commoncyls
    xpl = [displ.commoncyls_cod_cpara]';
    ypl = [displ.commoncyls_cod_cperp]';
else
    xpl = [displ.cod_cpara]';
    ypl = [displ.cod_cperp]';
end
xpl = xpl - mean(xpl);
ypl = ypl - mean(ypl);

% and put them in microns
px2mum = mean(displ(1).cyl_al_um./displ(1).cyl_al); % retrieve px2mum
xpl_um = xpl * px2mum;
ypl_um = ypl * px2mum;

% plot of center of drag
scatter(xpl_um, ypl_um,[],[displ.timestamp_s],'filled')
% plot forces
if flag_commoncyls
    quiver(xpl_um, ypl_um,...
        [displ.commoncyls_F_para_pN]',[displ.commoncyls_F_perp_pN]');
else
    quiver(xpl_um, ypl_um,...
        [displ.tot_F_para_pN]',[displ.tot_F_perp_pN]');
end

% also plot the basal body
if flag_commoncyls
    basal_xpl_um = (mean(arrayfun(@(i)plstr(i).xpl1(1),1:numel(plstr))) - mean([displ.commoncyls_cod_cpara]) )*px2mum;
    basal_ypl_um = (mean(arrayfun(@(i)plstr(i).ypl1(1),1:numel(plstr))) - mean([displ.commoncyls_cod_cperp]) )*px2mum;
else
    basal_xpl_um = (mean(arrayfun(@(i)plstr(i).xpl1(1),1:numel(plstr))) - mean([displ.cod_cpara]) )*px2mum;
    basal_ypl_um = (mean(arrayfun(@(i)plstr(i).ypl1(1),1:numel(plstr))) - mean([displ.cod_cperp]) )*px2mum;
end
plot( basal_xpl_um, basal_ypl_um, 'o', 'markerSize', 6)
text(basal_xpl_um, basal_ypl_um, '\leftarrow basal body')

% axis image
hacodpos = hacod.Position;

% define axis limits
% allxpl = vertcat(vertcat( plstr.xpl1 ), vertcat(plstr.xpl2));
% allxpl = (allxpl - mean([displ.commoncyls_cod_cpara])) .* px2mum;
% allypl = vertcat(vertcat( plstr.ypl1 ), vertcat(plstr.ypl2));
% allypl = (allypl - mean([displ.commoncyls_cod_cperp])) .* px2mum;
% plot(allxpl,allypl,'.')


% hacod.XLim = minmax( vertcat(allxpl, basal_xpl_um) );
% hacod.YLim = minmax( vertcat(allypl, basal_ypl_um) );
axis image
hacod.XLim = hacod.XLim + diff(minmax(hacod.XLim)) * 0.1 * [-1 +1];
hacod.YLim = hacod.YLim + diff(minmax(hacod.YLim)) * 0.1 * [-1 +1];

hacod.YAxisLocation = 'right';

hacod.XLabel.String = 'x_{cod}, [\mum]';
hacod.XLabel.FontSize = 14;
hacod.XLabel.Units = 'Normalized';
hacod.XLabel.Position(2) = hacod.XLabel.Position(2) + 0.05;

hacod.YLabel.String = 'y_{cod}, [\mum]';
hacod.YLabel.FontSize = 14;
hacod.YLabel.Units = 'Normalized';
hacod.YLabel.Position(1) = hacod.YLabel.Position(1) - 0.05;

% create colorbar with time, for scatterplot
hacodcb = colorbar('NorthOutside');

% put hacod back in correct position
hacod.Position = hacodpos;

% fix colorbar position
hacodcb.Units = 'Normalized';
hacodcb.Position([1 3]) = hacodpos([1 3]);
hacodcb.Position(2) = sum(hacodpos([2 4]));
hacodcb.Position(4) = hacodcb.Position(4)/2;


hacodcb.Label.String = 'Time, [s]';
hacodcb.Label.FontSize = 14;





%------------------- Add Centre of Drag to axes with displacements ---------------%

% axes in focus
axes(has)

if flag_commoncyls
    
    xpl = [displ.commoncyls_cod_cpara]' + xoffsets' - [plstr.meanxpl]';
    ypl = [displ.commoncyls_cod_cperp]';
    
    if ~flag_cylforces
        quiver(xpl,ypl,[displ.commoncyls_F_para_pN]',[displ.commoncyls_F_perp_pN]','Color',has.ColorOrder(1,:),'Clipping','off');
    else
        
        % initialise matrices
        cxpl_mat = nan(max(arrayfun(@(i)numel(plstr(i).cxpl),1:numel(plstr))), numel(plstr));
        cypl_mat = cxpl_mat;
        fxpl_mat = cxpl_mat;
        fypl_mat = cxpl_mat;
        idx_commoncyls_mat = false(size(cxpl_mat));
        
        % fill them
        for i = 1:numel(plstr)
            n = numel(plstr(i).cxpl);
            cxpl_mat(1:n,i) = plstr(i).cxpl;
            cypl_mat(1:n,i) = plstr(i).cypl;
            fxpl_mat(1:n,i) = plstr(i).fxpl;
            fypl_mat(1:n,i) = plstr(i).fypl;
            idx_commoncyls_mat(1:n,i) = displ(i).idx_commoncyls;
            
        end %for
                
        % apply offset 
        cxpl_mat = cxpl_mat + xoffsets - [plstr.meanxpl];
        
        quiver(cxpl_mat(idx_commoncyls_mat), cypl_mat(idx_commoncyls_mat),...
            fxpl_mat(idx_commoncyls_mat), fypl_mat(idx_commoncyls_mat), 'b' );
        
    end %if
else
    xpl = [displ.cod_cpara]' + xoffsets' - [plstr.meanxpl]';
    ypl = [displ.cod_cperp]';
    
    if ~flag_cylforces
        quiver(xpl,ypl,[displ.tot_F_para_pN]',[displ.tot_F_perp_pN]','Color','r','Clipping','off');
    else
        
        % initialise matrices
        cxpl_mat = nan(max(arrayfun(@(i)numel(plstr(i).cxpl),1:numel(plstr))), numel(plstr));
        cypl_mat = cxpl_mat;
        fxpl_mat = cxpl_mat;
        fypl_mat = cxpl_mat;
        
        % fill them
        for i = 1:numel(plstr)
            n = numel(plstr(i).cxpl);
            cxpl_mat(1:n,i) = plstr(i).cxpl;
            cypl_mat(1:n,i) = plstr(i).cypl;
            fxpl_mat(1:n,i) = plstr(i).fxpl;
            fypl_mat(1:n,i) = plstr(i).fypl;
            
        end %for
        
        % apply offset 
        cxpl_mat = cxpl_mat + xoffsets - [plstr.meanxpl];
        
        % plot
        quiver(cxpl_mat, cypl_mat, fxpl_mat, fypl_mat, 'b' );
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Again but in other figure

% keyboard
%%
%------------------ Axes with trajectory of centre of drag ---------------------%

[hf2, hacod] = figure_template;


% plot of center of drag
scatter(xpl_um, ypl_um,[],[displ.timestamp_s],'filled')
% plot forces
if flag_commoncyls
    quiver(xpl_um, ypl_um,...
        [displ.commoncyls_F_para_pN]',[displ.commoncyls_F_perp_pN]');
else
    quiver(xpl_um, ypl_um,...
        [displ.tot_F_para_pN]',[displ.tot_F_perp_pN]');
end

% also plot the basal body
plot( basal_xpl_um, basal_ypl_um, 'o', 'markerSize', 6)
ht = text(basal_xpl_um, basal_ypl_um, 'basal body \rightarrow ');
ht.HorizontalAlignment = 'right';

% axis image
hacodpos = hacod.Position;

% define axis limits
% allxpl = vertcat(vertcat( plstr.xpl1 ), vertcat(plstr.xpl2));
% allxpl = (allxpl - mean([displ.commoncyls_cod_cpara])) .* px2mum;
% allypl = vertcat(vertcat( plstr.ypl1 ), vertcat(plstr.ypl2));
% allypl = (allypl - mean([displ.commoncyls_cod_cperp])) .* px2mum;
% plot(allxpl,allypl,'.')


hacod.XLim = minmax( vertcat(xpl_um, basal_xpl_um) );
hacod.YLim = minmax( vertcat(ypl_um, basal_ypl_um) );
% axis image
hacod.XLim = hacod.XLim + diff(minmax(hacod.XLim)) * 0.1 * [-1 +1];
hacod.YLim = hacod.YLim + diff(minmax(hacod.YLim)) * 0.1 * [-1 +1];

% hacod.YAxisLocation = 'right';

hacod.XLabel.String = 'x_{cod}, [\mum]';
hacod.XLabel.Position(2) = -0.155;
hacod.YLabel.String = 'y_{cod}, [\mum]';

% create colorbar with time, for scatterplot
% hacodcb = colorbar('NorthOutside');
% 
% % put hacod back in correct position
% hacod.Position = hacodpos;
% 
% % fix colorbar position
% hacodcb.Units = 'Normalized';
% hacodcb.Position([1 3]) = hacodpos([1 3]);
% hacodcb.Position(2) = sum(hacodpos([2 4]));
% hacodcb.Position(4) = hacodcb.Position(4)/2;
% 
% 
% hacodcb.Label.String = 'Time, [s]';


%% cilium shape, little forces

% axes positions
[hf3, hatf] = figure_template;


hatf.XTick = [];
hatf.YTick = [];
hatf.Visible = 'off';

hold on

cmap = parula(numel(displ)+1);

% for loop on plotting forces
for i = 1:numel(displ)
    
    % first plot cilia
%     plot(plstr(i).xpl1, plstr(i).ypl1, 'color',[0.85 0.325 0.098], 'Clipping','off');
    plot(plstr(i).xpl1, plstr(i).ypl1,...
        'color',cmap(i,:),...
        'Clipping','off',...
        'LineWidth',2);
    
    % special case: last position
    if i == numel(displ)
        plot(plstr(i).xpl2, plstr(i).ypl2,...
            'color',cmap(end,:),...
            'Clipping','off',...
            'LineWidth',2);
    end %if
    
    
    %then plot forces of each cylinder
    [plstr(i).cxpl, plstr(i).cypl] = xy2paraperp(baseline, displ(i).cil(1).cyl_cx, displ(i).cil(1).cyl_cy);
    [plstr(i).fxpl, plstr(i).fypl] = xy2paraperp(baseline, displ(i).cyl_F_x_pN, displ(i).cyl_F_y_pN);
    
    if flag_commoncyls
        quiver(plstr(i).cxpl(displ(i).idx_commoncyls), plstr(i).cypl(displ(i).idx_commoncyls),...
            10*plstr(i).fxpl(displ(i).idx_commoncyls), 10*plstr(i).fypl(displ(i).idx_commoncyls),...
            0,...
            'color','r');
    else
        quiver(plstr(i).cxpl, plstr(i).cypl,...
            plstr(i).fxpl, plstr(i).fypl ,...
            0,...
            'color',[0 0.447 0.741]);
    end
    
end %for

axis image






end
















