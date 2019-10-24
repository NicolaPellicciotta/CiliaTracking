function [ clclk ] = clclk_calculate_force( clclk_int_fullfilename, cyl_length_um, fluid_viscosity )
%clclk_calculate_force Calculate the force that a cilium exerts on the
%fluid by using RFT
%   initialy developped asa a script (force_attempts), relies on the
%   existence of a clclk_int file. In case only the clclk_file is present,
%   the user should run the clclk_find_intersections function first


%% input check and parameters

if nargin < 3 || isempty(fluid_viscosity)
    fluid_viscosity = 1e-3; % Pa*s
end


if nargin < 2 || isempty(cyl_length_um)
    cyl_length_um = 0.8; %this is just a suggestion, will be adjusted on a displ-to-displ basis to prevent excessively short segments
end

if nargin < 1 || isempty(clclk_int_fullfilename)
    [clclk_int_filename, clclk_int_pathname] = uigetfile('*.clclk_int');
else
    [clclk_int_filename, clclk_int_pathname] = parse_filename(clclk_int_fullfilename); % this adds the path if a local file was given as an input
end

clclk_int_fullfilename = fullfile(clclk_int_pathname, clclk_int_filename);



%% load data

clclk = load(clclk_int_fullfilename,'-mat'); %remove idx_clicked_frames(25)


%% find parallel and perpendicular component wrt the cell surface

% I have (Fx, Fy) at each time now, but what matters is Ft, Fn wrt to the
% cell surface. Need to find tangent and normal using the lines from the
% clclk_int file

baseline = find_baseline(clclk.Intersections, clclk.Stroke);


%% manually select bad frames

idx_clicked_frames = [clclk.clicked_frames.frame_number];
disp(['There are ',num2str(numel(idx_clicked_frames)),' clicked frames.']);

idx_bad_frames = clclk_scroll_find_misclicking(clclk);
disp(['You selected ',num2str(sum(idx_bad_frames)),' bad frames.']);

idx_clicked_frames(idx_bad_frames) = [];

disp(['Now there are ',num2str(numel(idx_clicked_frames)),' clicked frames.']);


%% for loops on couples of frames


px2mum = clclk.Stroke.px2mum;

parfor cfc = 1:numel(idx_clicked_frames)-1
    
    % displacement is calculated only between couples of frames. Can use
    % this to decide which "overlapping" region has been clicked" best way
    % is probably to find the offset on one of the two cilia so that max
    % distance travelled by correspondent points is minimised.
    
    clicked_frames = clclk.clicked_frames(idx_clicked_frames(cfc:cfc+1));
    points         = clclk.points(idx_clicked_frames(cfc:cfc+1));
    cil = prepare_cil_struct_alt(clicked_frames, points, px2mum);
     
    % this would be same thing, less memory friendly but more
    % "script-friendly":
%     cil  = prepare_cil_struct( clclk.clicked_frames, clclk.points, ...
%         idx_clicked_frames(cfc:cfc+1), px2mum );
    
    
    % try and set an offset, have to do it on both cilia
    [ooc, cto] = find_offset_on_cilium( cil );  % offset on cilium, cilium to offset
    
    
    % apply the offset on the appropriate cilium
    cil  = apply_offset_to_cilium( cil, ooc, cto );
    
    
    % calculate cylinder's parameters
    cil = calculate_cylinders_parameters( cil, px2mum, cyl_length_um );
    
    
    % calculate displacement's properties
    displ(cfc) = calculate_displacements_properties(cil, baseline, fluid_viscosity, px2mum);
    

%     hf = cil_plot(cil,displ(cfc), baseline);
%     pause(0.01);
%     close(hf)
    
end



%% now a check on displ:
% because of experimental error, cilium's length "changes" overtime.
% Find the minimum length and use that to clip the integral of the force

displ = calculate_commoncyls_properties( displ, baseline );


%% calculate now the centre-of-mass of the cilium, and the drag-centre
% caclulate centre-of-mass of the cilium limited to the commoncylinders,
% and the centre of the force applied (first moment of the force along the
% cilium/common cylinders)

displ = calculate_com_cod(displ, baseline);


%% write variables in structure clclk

% Force variables
clclk.Force.baseline            = baseline;
clclk.Force.displ               = displ;
clclk.Force.cyl_length_um       = cyl_length_um;
clclk.Force.fluid_viscosity     = fluid_viscosity;

% name of clclk_int, usually useful to track changes
clclk.clclk_int_fullfilename = clclk_int_fullfilename;

% clicked frames, to make it easier next time I need it
clclk.idx_clicked_frames = idx_clicked_frames;


end %function






%{

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


    
%}
