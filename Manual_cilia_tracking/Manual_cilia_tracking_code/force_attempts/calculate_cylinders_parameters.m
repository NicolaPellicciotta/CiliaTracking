function [ cil ] = calculate_cylinders_parameters( cil, px2mum, cyl_length_um)
%calculate_cylinders_parameters Divide cilium in cylinders, finds midpoints of tiny cylinders, and
%their actual length
%   so far we've been using cil.xr, cil.yr, cyl.tr, spaced of 1um  along
%   the cilium. This is not important, even if we want to be able to
%   control the size of the cilium, because those point were only used to
%   align the cilia (with the offset). This is the function that actually
%   defines the cylinders, so here's where cyl_length_m matters


% calculate the shortest between the portion of the two cilia
min_arclength_um = min( cil(1).tt_um(end), cil(2).tt_um(end) );

%{
%these lines are not really needed in case of overlapping cylinders:

% find how many integers cyl_length_um fit in min_arclength_um
number_cylinders = floor(min_arclength_um/cyl_length_um);

% now adjust cyl_length_um so that min_arclength_um contains an integer
% numbers of cylinders (this is to avoid having very short cylinders where
% RFT is not valid anymore (drag coefficient diverges)
cyl_length_um = min_arclength_um / number_cylinders;
%} 


% find the (integer multiple of cyl length um) points common between the (portion of) cilium i and
% cilium j

um_start_points_to_find = (0:cyl_length_um/10:min_arclength_um - cyl_length_um)';
um_end_points_to_find = um_start_points_to_find + cyl_length_um;



% for both cilia, 
for i = 1:numel(cil)

    %now find the points whose arclength is the closest to
    % um_points_to_find
    cil(i).idx_umm = horzcat(knnsearch(cil(i).tt_um, um_start_points_to_find),...
        knnsearch(cil(i).tt_um, um_end_points_to_find));
    
    % use the indices just found to find the extremes of the cylinders
    % first column is start points, second is end points
    cil(i).tr_um = reshape(cil(i).tt_um(cil(i).idx_umm(:)),[],2);
    cil(i).tr = reshape(cil(i).tt(cil(i).idx_umm(:)),[],2);
    cil(i).xr = reshape(cil(i).xx(cil(i).idx_umm(:)),[],2);
    cil(i).yr = reshape(cil(i).yy(cil(i).idx_umm(:)),[],2);
    
    
    % centre of cylinders is midpoint
    cil(i).cyl_cx = mean(cil(i).xr, 2); % cylinder's centre x
    cil(i).cyl_cy = mean(cil(i).yr, 2); % cylinder's centre y
    
    % arclength of each "cylinder"
    cil(i).cyl_al = diff(cil(i).tr,1,2);
    cil(i).cyl_al_um = diff(cil(i).tr_um,1,2); % in um
    
    % linear length of each cylinder (using endpoints coordinates)
    cil(i).cyl_ll = hypot(diff(cil(i).xr,1,2), diff(cil(i).yr,1,2));
    cil(i).cyl_ll_um = cil(i).cyl_ll .* px2mum; % in um
    
end %for

end

