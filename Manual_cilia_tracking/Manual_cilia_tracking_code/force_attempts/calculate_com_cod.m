function [ displ_vec ] = calculate_com_cod( displ_vec, baseline )
%calculate_com_cod % caclulate centre-of-mass of the cilium limited to the commoncylinders,
% and the centre of the force applied (first moment of the force along the
% cilium/common cylinders)


flag_commoncyls = isfield(displ_vec,'idx_commoncyls');


for cfc = 1:numel(displ_vec)
    
    
    % loop on cilia for cilium-specific properties
    for i = 1:2
        
        % centre of mass of cilium, after applying offset
        displ_vec(cfc).cil(i).com_xx = mean(displ_vec(cfc).cil(i).xx(displ_vec(cfc).cil(i).tt_um > 0));
        displ_vec(cfc).cil(i).com_yy = mean(displ_vec(cfc).cil(i).yy(displ_vec(cfc).cil(i).tt_um > 0));
        
        % centre of mass of cilium, in para/perp coordinate
        [displ_vec(cfc).cil(i).com_para, displ_vec(cfc).cil(i).com_perp] =...
            xy2paraperp(baseline,displ_vec(cfc).cil(i).com_xx, displ_vec(cfc).cil(i).com_yy);
        
        % now calculate the com along the cilium but stopping at the same value
        % we used to stop the commoncylinders
        if flag_commoncyls
            
            % find cutoff length
            cutoff_al_um = sum(displ_vec(cfc).cyl_ll_um(displ_vec(cfc).idx_commoncyls)); %this value is actually common to the 2 cilia snapshot
            
            % find valid arclength values
            idx_tt_um = displ_vec(cfc).cil(i).tt_um > 0 & displ_vec(cfc).cil(i).tt_um < cutoff_al_um;
            
            % find com
            displ_vec(cfc).cil(i).commoncyls_com_xx = mean(displ_vec(cfc).cil(i).xx(idx_tt_um));
            displ_vec(cfc).cil(i).commoncyls_com_yy = mean(displ_vec(cfc).cil(i).yy(idx_tt_um));
            
            % now in para/perp coordinates
            [displ_vec(cfc).cil(i).commoncyls_com_para, displ_vec(cfc).cil(i).commoncyls_com_perp] =...
                xy2paraperp(baseline,displ_vec(cfc).cil(i).commoncyls_com_xx, displ_vec(cfc).cil(i).commoncyls_com_yy);
        end %if
        
    end
    
    
    % now I have to do the same thing for the "centre of drag". cfr
    % Brumley2014 for the definition of centre of drag along each direction
    % as basically the first moment of the force along that direction
    % (sum(xi |fi|)/sum(|fi|)). The amplitude of the force exerted by each
    % cylinder was computed in calculate_displacement_properties. Unlike
    % the cilium centre of mass, as this is defined using the force, it is
    % a property of the displ struct, not disp.cil
    displ_vec(cfc).cod_cx = sum( displ_vec(cfc).cil(1).cyl_cx .* displ_vec(cfc).cyl_Fampl_pN )/...
        sum(displ_vec(cfc).cyl_Fampl_pN);
    displ_vec(cfc).cod_cy = sum( displ_vec(cfc).cil(1).cyl_cy .* displ_vec(cfc).cyl_Fampl_pN )/...
        sum(displ_vec(cfc).cyl_Fampl_pN);
    
    % convert these in para/perp too
    [displ_vec(cfc).cod_cpara, displ_vec(cfc).cod_cperp] = ...
        xy2paraperp(baseline, displ_vec(cfc).cod_cx, displ_vec(cfc).cod_cy);
    
    
    % now again but using only the cylinders in commoncyls
    if flag_commoncyls
        
        displ_vec(cfc).commoncyls_cod_cx = sum( displ_vec(cfc).cil(1).cyl_cx(displ_vec(cfc).idx_commoncyls) .* displ_vec(cfc).cyl_Fampl_pN(displ_vec(cfc).idx_commoncyls) )/...
            sum(displ_vec(cfc).cyl_Fampl_pN(displ_vec(cfc).idx_commoncyls));
        displ_vec(cfc).commoncyls_cod_cy = sum( displ_vec(cfc).cil(1).cyl_cy(displ_vec(cfc).idx_commoncyls) .* displ_vec(cfc).cyl_Fampl_pN(displ_vec(cfc).idx_commoncyls) )/...
            sum(displ_vec(cfc).cyl_Fampl_pN(displ_vec(cfc).idx_commoncyls));
        [displ_vec(cfc).commoncyls_cod_cpara, displ_vec(cfc).commoncyls_cod_cperp] =...
            xy2paraperp(baseline, displ_vec(cfc).commoncyls_cod_cx, displ_vec(cfc).commoncyls_cod_cy);
        
    end %if
    
end %for cfc

end

