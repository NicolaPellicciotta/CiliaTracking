function [ hf ] = cil_plot( cil, displ, baseline )
%cil_plot quick function to plot the two cilia in the cil structure

if nargin < 3 %use just x,y direction for plotting
    baseline.parav_x = 1;
    baseline.parav_y = 0;
    baseline.perpv_x = 0;
    baseline.perpv_y = 1;
end

if nargin < 2
    displ = [];
end


hf = figure;
hold on;
box on;

axis image

lc = get(gca,'ColorOrder');

ha = gca;


for i = 1:2
    
    % find negative bits
    idx_neg = cil(i).tt < 0;
    idx_pos = ~idx_neg;
    
    % plot them
    [par,per] = xy2paraperp(baseline, cil(i).xx(idx_neg), cil(i).yy(idx_neg));
    plot( par, per, 'k--'  );%,'LineWidth',2);
    
    % plot smooth cilium
%     [par, per] = xy2paraperp(baseline, cil(i).xx(idx_pos), cil(i).yy(idx_pos));
%     hpsc(i) = plot( par,per, 'color', lc(i,:) );%,'LineWidth',2);
    
    % plot markers for cylinder endpoints
    [par, per] = xy2paraperp(baseline, cil(i).xr, cil(i).yr);
    hpsc(i,:) = plot( par', per', '.-', 'color', 0.8 .* lc(i,:),'LineWidth',2);
    
    
    
    % plot markers for cylinders centre points
%     if isfield(cil,'cyl_cx')
%         [par, per] = xy2paraperp(baseline, cil(i).cyl_cx, cil(i).cyl_cy);
%         plot( par, per, 'x', 'color', 0.8 .* lc(i,:),'MarkerSize',4);
%     end %if
    
end %for

% use plot with matrices. plot(X,Y) plots columns of X vs columns of Y. My
% X columns have to be one value of xr on a cilium, and the correspondent
% on the other one. same for Y

% [par, per] = xy2paraperp([cil.xr]', [cil.yr]');
% plot( par, per, 'g')

if ~isempty(displ)
    
    sc = 12;
    [parpos, perpos] = xy2paraperp(baseline, cil(1).cyl_cx, cil(1).cyl_cy);
    [partv, pertv] = xy2paraperp(baseline,  sc*displ.cyl_tv_x, sc*displ.cyl_tv_y );
    [parnv, pernv] = xy2paraperp(baseline,  sc*displ.cyl_nv_x, sc*displ.cyl_nv_y );
    [pards, perds] = xy2paraperp(baseline,  displ.cyl_displ_x, displ.cyl_displ_y ); % displacement
    [paraf, perpf] = xy2paraperp(baseline,  displ.cyl_F_x_pN, displ.cyl_F_y_pN );
%     quiver(parpos, perpos, partv, pertv ,0,'g','LineWidth',1.2);
%     quiver(parpos, perpos, parnv, pernv ,0,'r','LineWidth',1.2);
    quiver(parpos, perpos, pards, perds ,0,'b','LineWidth',1.2);
    quiver(parpos, perpos, sc*paraf, sc*perpf ,0,'k','LineWidth',1.2);
end

hleg = legend(hpsc(:,1),{'1st','2nd'},'Box','off','FontSize',12,'Location','SouthEast')


ha.XTick = [];
ha.YTick = [];
ha.XColor = 'w';
ha.YColor = 'w';

end


