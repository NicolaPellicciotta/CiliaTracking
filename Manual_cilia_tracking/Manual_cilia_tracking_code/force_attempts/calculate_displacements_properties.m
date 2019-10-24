function [ displ ] = calculate_displacements_properties( cil, baseline, mu, px2mum )
%calculate_displacements_properties given the cil structure, finds the
%displacements between corresponding cylinders, then decompose in (locally)
%tangent and normal components, calculate the force each cylinder exerts on
%the fluid, and then sum up the contributions of all cylinders
%   Detailed explanation goes here

%% input check and parameters

if nargin < 4 || isempty(px2mum)
    warning('um/px will be estimated by existing properties in the cil structure.')
    px2mum = mean(cil(1).cyl_al_um./cil(1).cyl_al);
end

if nargin < 3 || isempty(mu)
    warning('No viscosity was entered. The viscosity of water at 20C will be used.')
    mu = 1e-3; % viscosity of water at 20C, Pa s
end


% radius of a cilium
cil_radius_um = 0.1; % 0.1 um = 100 nm

%% actual function to create displ(acement) structure

% calculate displacements between centers of cylinders
displ.cyl_displ_x = diff([cil.cyl_cx],1,2); % displacement of cylinder's centre, x component
displ.cyl_displ_y = diff([cil.cyl_cy],1,2); % displacement of cylinder's centre, y component
displ.cyl_displ_x_um = displ.cyl_displ_x .* px2mum;
displ.cyl_displ_y_um = displ.cyl_displ_y .* px2mum;


% find time (average between times in cil) and time interval (diff between
% times in cil)
displ.time_interval_s = diff([cil.timestamp]); % dt, time between the two frames in cil
displ.timestamp_s = mean([cil.timestamp]);


% copy here the length of the cylinders (cil(1).cyl_ll_um and
% cil(1).cyl_al_um)
displ.cyl_al = cil(1).cyl_al;
displ.cyl_al_um = cil(1).cyl_al_um;
displ.cyl_ll = cil(1).cyl_ll_um;
displ.cyl_ll_um = cil(1).cyl_ll_um;

% how long was the cilium?
displ.tot_al_um = min( cil(1).tt_um(end), cil(2).tt_um(end) );


% find local tangent and normal, in [x,y] coordinates
% given the approximation of the cilium with cylinders, the tangent will be
% the direction of the segment on the **first** cilium, and the normal the perpendicular to it
displ.cyl_tv_x = diff(cil(1).xr,1,2) ./ cil(1).cyl_ll; % x coordinate of the (local) tangential versor
displ.cyl_tv_y = diff(cil(1).yr,1,2) ./ cil(1).cyl_ll; % x coordinate of the (local) tangential versor
displ.cyl_nv_x = -displ.cyl_tv_y;
displ.cyl_nv_y = displ.cyl_tv_x;


% now write displacement in tangential, normal coordinates
% for each cylinder, the tangential component of the displacement is the
% scalar product between the displacement_xy and the local tangent versor.
% I have to do it for each cylinder though, hence the vectorial (on the
% cylinders) notation
displ.cyl_displ_t = sum([displ.cyl_displ_x,displ.cyl_displ_y] .* [displ.cyl_tv_x,displ.cyl_tv_y], 2);
displ.cyl_displ_n = sum([displ.cyl_displ_x,displ.cyl_displ_y] .* [displ.cyl_nv_x,displ.cyl_nv_y], 2);
displ.cyl_displ_t_um = displ.cyl_displ_t .* px2mum;
displ.cyl_displ_n_um = displ.cyl_displ_n .* px2mum;


% calculate the parallel (tangential) and perpendicular (normal) drag
% coefficients using RFT
displ.cyl_dc_t_Pas = 4*pi*mu ./ ( log(displ.cyl_ll_um.^2./cil_radius_um^2) - 1 ); % cylinder's parallel drag coefficient
displ.cyl_dc_n_Pas = 8*pi*mu ./ ( log(displ.cyl_ll_um.^2./cil_radius_um^2) + 1 ); % cylinder's perpendicular drag coefficient


% calculate the parallel (tangential) and perpendicular (normal) force
% exerted on the fluid by each cylinder
% F (pN) = f*L = drag coefficient (Pa s) * velocity (um/s) * length of cylinder (um)
displ.cyl_F_t_pN = displ.cyl_dc_t_Pas .* displ.cyl_displ_t_um ./ displ.time_interval_s .* displ.cyl_ll_um; % in pN, F = f*L = F/L*L = drag coefficient * velocity * length of cilinder
displ.cyl_F_n_pN = displ.cyl_dc_n_Pas .* displ.cyl_displ_n_um ./ displ.time_interval_s .* displ.cyl_ll_um; % in pN, F = f*L = F/L*L = drag coefficient * velocity * length of cilinder


% each cylinder will exert a force F = F_t \hat{t} + F_n\hat{n}.
% have to convert each (F_t, F_v) in (F_x, F_y) before summing them up

% if matrix to change base is M = [tx nx; ty ny] and (vt,vn) = (vx,vy) * M
% then to go back in xy coordinates you fo (vx, vy) = (vt,vn) * M^-1.
% Pen and paper calculations yield M^-1 = [ny -nx; -ty tx] (det(M) == 1).
% Actually, since nx = -ty and ny = tx, M = [tx -ty; ty tx] and
% M^-1 = M' = [tx ty; -ty tx] (but not using this for sake of clarity

% so for generic vector v = (vt, vn) => vx = vt*ny - vn*ty, vy = -vt*nx + vn*tx
displ.cyl_F_x_pN =  displ.cyl_F_t_pN .* displ.cyl_nv_y - displ.cyl_F_n_pN .*  displ.cyl_tv_y;
displ.cyl_F_y_pN = -displ.cyl_F_t_pN .* displ.cyl_nv_x + displ.cyl_F_n_pN .*  displ.cyl_tv_x;


% now finally "integrate" (aka sum) all the contributions of each cylinder
% take care of overlapping cylinders by dividing the sum of the forces by
% (sum of cylinders lengths/ cilium length)
displ.tot_F_x_pN = sum(displ.cyl_F_x_pN) / sum(displ.cyl_al_um) * displ.tot_al_um;
displ.tot_F_y_pN = sum(displ.cyl_F_y_pN) / sum(displ.cyl_al_um) * displ.tot_al_um;


% now change base so that we have the toal force in the FoR defined by the
% surface of the cell
displ.tot_F_para_pN = displ.tot_F_x_pN * baseline.parav_x + displ.tot_F_y_pN * baseline.parav_y;
displ.tot_F_perp_pN = displ.tot_F_x_pN * baseline.perpv_x + displ.tot_F_y_pN * baseline.perpv_y;


% store for later purposes/normalisations how many cylinders we used for
% calculating the total force (actually, not the number but the sum of the
% lengths of the cylinders we used)
displ.tot_ll_um = sum(displ.cyl_ll_um);


% it is also useful to calculate the amplitude (indep on the base) of the
% force that each cylinder exerts on the fluid (sum of the 2 contributions)
displ.cyl_Fampl_pN = hypot(displ.cyl_F_t_pN, displ.cyl_F_n_pN); % could have done with x,y components too


% copy also the cil structure, so I don't have to keep track of it all the
% time
displ.cil = cil;

end

