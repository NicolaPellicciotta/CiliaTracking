function [ displ_vec ] = calculate_commoncyls_properties( displ_vec, baseline )
%calculate_commoncyls_properties not all ciliasnapshot are of the same
%length. This accounts for it, by finding the cylinders that are common to
%all the cilia (as in: cut off the end on longer ones, do not care about basal part)
%   because of experimental error, cilium's length "changes" overtime.
%   Find the minimum length and use that to clip the integral of the force


min_cil_al_um = min([displ_vec(:).tot_al_um]);
common_limit = ((min_cil_al_um - mean(vertcat(displ_vec.cyl_ll_um)))/mean(diff(displ_vec(1).cil(1).tr_um(:,1))) + 1) *...
    mean(vertcat(displ_vec.cyl_ll_um));
% if you cover something L long with several l long segments every s, the
% cumulative length will be ((L-l)/s + 1)*l


% use it to choose how many cylinders to use
for i = 1:numel(displ_vec)
    
    % find cylinders that fall below limit
    displ_vec(i).idx_commoncyls = cumsum(displ_vec(i).cyl_ll_um) <= common_limit + 0.6;
    
    % re-calculate the integral
    displ_vec(i).commoncyls_F_x_pN = sum(displ_vec(i).cyl_F_x_pN(displ_vec(i).idx_commoncyls)) / ...
        sum(displ_vec(i).cyl_al_um(displ_vec(i).idx_commoncyls)) * min_cil_al_um;
    displ_vec(i).commoncyls_F_y_pN = sum(displ_vec(i).cyl_F_y_pN(displ_vec(i).idx_commoncyls)) / ...
        sum(displ_vec(i).cyl_al_um(displ_vec(i).idx_commoncyls)) * min_cil_al_um;

    
    % now change base so that we have the total force on the common cylinders in the FoR defined by the
    % surface of the cell
    displ_vec(i).commoncyls_F_para_pN = displ_vec(i).commoncyls_F_x_pN * baseline.parav_x + displ_vec(i).commoncyls_F_y_pN * baseline.parav_y;
    displ_vec(i).commoncyls_F_perp_pN = displ_vec(i).commoncyls_F_x_pN * baseline.perpv_x + displ_vec(i).commoncyls_F_y_pN * baseline.perpv_y;
    
    % now I expect these values (sum of the commoncylinders) to fluctuate within 1 um (lenth of
    % cylinders)
    displ_vec(i).commoncyls_ll_um = sum(displ_vec(i).cyl_ll_um(displ_vec(i).idx_commoncyls));
    
end

end

