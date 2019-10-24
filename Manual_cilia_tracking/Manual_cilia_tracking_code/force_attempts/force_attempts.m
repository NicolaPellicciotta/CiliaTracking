clear all
close all

% cd E:\Data\Cilia_Profile\Clicked_Cilia\N_clclk
cd D:\Dropbox\Universita\PhD\Manual_cilia_tracking_code\force_attempts

%%
% fl = dir('*.clclk_int');
%
% fc = 5;
% % fc = 8; %file 5 and 8 of the normals look good (aka the base is a single
% % point)
%
%
% load(fl(fc).name, '-mat');

% load('test.22Jul2015_15.58.09.clclk_int','-mat');

% load('E:\Data\Cilia_Profile\Clicked_Cilia\N_clclk\N.22Jul2015_15.33.21.clclk_int','-mat');
% load('E:\Data\Cilia_Profile\Clicked_Cilia\N_clclk\N.22Jul2015_15.43.43.clclk_int','-mat');
% load('E:\Data\Cilia_Profile\Clicked_Cilia\N_clclk\test.22Jul2015_15.54.43.clclk_int','-mat');
% load('E:\Data\Cilia_Profile\Clicked_Cilia\N_clclk\test.22Jul2015_15.54.43b.clclk_int','-mat');
load('E:\Data\Cilia_Profile\Clicked_Cilia\N_clclk\test.22Jul2015_15.58.09.clclk_int','-mat'); %remove idx_clicked_frames(25)
% load('E:\Data\Cilia_Profile\Clicked_Cilia\N_clclk\test.22Jul2015_16.07.35.clclk_int','-mat');
% load('E:\Data\Cilia_Profile\Clicked_Cilia\N_clclk\test.22Jul2015_16.14.46.clclk_int','-mat');
% load('E:\Data\Cilia_Profile\Clicked_Cilia\Nm_clclk\Nm.02Nov2015_09.40.59.clclk_int','-mat')


% profile on

%% parameters

water_viscosity = 1e-3; % Pa*s
cyl_length_um = .8; %this is just a suggestion, will be adjusted on a displ-to-displ basis to prevent excessively short segments

%%

idx_clicked_frames = [clicked_frames.frame_number];


hf = figure;
hold on
set(gca,'ColorOrder',parula(numel(idx_clicked_frames)));
set(gcf,'Color','k');
set(gca,'Color','k')
for i = idx_clicked_frames
    plot(points(i).dw_cilium_xx, points(i).dw_cilium_yy,'LineWidth',1.2);
end

% pause;
close(hf);
% close(fc)

%-----------------------------

%% find parallel and perpendicular component wrt the cell surface

% I have (Fx, Fy) at each time now, but what matters is Ft, Fn wrt to the
% cell surface. Need to find tangent and normal using the lines from the
% clclk_int file

baseline = find_baseline(Intersections, Stroke);


%% for loops on couples of frames

% idx_clicked_frames(25) = [];
warning('frame 112 manually skipped');

for cfc = 1:numel(idx_clicked_frames)-1
    
    % displacement is calculated only between couples of frames. Can use
    % this to decide which "overlapping" region has been clicked" best way
    % is probably to find the offset on one of the two cilia so that max
    % distance travelled by correspondent points is minimised.
    
     
    cil  = prepare_cil_struct( clicked_frames, points, idx_clicked_frames(cfc:cfc+1), Stroke.px2mum );
    
    
    % try and set an offset, have to do it on both cilia
    [ooc, cto] = find_offset_on_cilium( cil );  % offset on cilium, cilium to offset
%     disp([num2str(ooc) ', ' num2str(cto)])
    
    
    % apply the offset on the appropriate cilium
    cil  = apply_offset_to_cilium( cil, ooc, cto );
    
    
    % calculate cylinder's parameters
    cil = calculate_cylinders_parameters( cil, Stroke.px2mum, cyl_length_um );
    
    
    % calculate displacement's properties
    displ(cfc) = calculate_displacements_properties(cil, baseline, water_viscosity, Stroke.px2mum);
    

%     hf = cil_plot(cil,displ(cfc), baseline);
%     pause(0.01);
%     close(hf)
    
end

% profile viewer


%% now a check on displ:
% because of experimental error, cilium's length "changes" overtime.
% Find the minimum length and use that to clip the integral of the force

displ = calculate_commoncyls_properties( displ, baseline );


%% calculate now the centre-of-mass of the cilium, and the drag-centre
% caclulate centre-of-mass of the cilium limited to the commoncylinders,
% and the centre of the force applied (first moment of the force along the
% cilium/common cylinders)


displ = calculate_com_cod(displ, baseline);




%%

figure
hold on

plot([displ.cod_cx],[displ.cod_cy],'r')
plot([displ.commoncyls_cod_cx],[displ.commoncyls_cod_cy],'b')

quiver([displ.cod_cx],[displ.cod_cy], [displ.tot_F_x_pN], [displ.tot_F_y_pN],'r')
quiver([displ.commoncyls_cod_cx],[displ.commoncyls_cod_cy], [displ.commoncyls_F_x_pN], [displ.commoncyls_F_y_pN],'b')


figure
hold on

plot([displ.cod_cpara],[displ.cod_cperp],'r.')
plot([displ.commoncyls_cod_cpara],[displ.commoncyls_cod_cperp],'bo')

quiver([displ.cod_cpara],[displ.cod_cperp], [displ.tot_F_para_pN], [displ.tot_F_perp_pN],'r')
quiver([displ.commoncyls_cod_cpara],[displ.commoncyls_cod_cperp], [displ.commoncyls_F_para_pN], [displ.commoncyls_F_perp_pN],'b')


    

