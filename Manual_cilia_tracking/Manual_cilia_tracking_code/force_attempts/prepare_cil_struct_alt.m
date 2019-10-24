function [ cil ] = prepare_cil_struct_alt( clicked_frames, points, px2mum )
%prepare_cil_struct prepares the 2-cilia-structure to find the best offset to
%compare cilia and eventually calculate the force
%   in this version of the funciton, I assume that I already only have 2
%   frames passed as an input

%% input check

if numel(clicked_frames) ~= 2 || numel(points) ~= 2
    error('frames_numbers_to_analyse has to be a 2-elements vector');
end

% sort if necessary the struct, in case of user error
if diff([clicked_frames.frame_number]) < 0
    [~, idx] = sort([clicked_frames.frame_number]);
    clicked_frames = clicked_frames(idx);
    points = points(ind);
    warning('switching frames');
end

%% flags

flag_plot = 0;


%% build structure

for i = 2:-1:1
    
    % x and y interpolated coordinates
    cil(i).xx = points(i).dw_cilium_xx;
    cil(i).yy = points(i).dw_cilium_yy;
    
    % arclength
    cil(i).tt = cumsum(sqrt([0;diff(cil(i).xx)].^2 + [0;diff(cil(i).yy)].^2));
    cil(i).tt_um = cil(i).tt.*px2mum;
    
    
    % use knnsearch to find points closer to integer
    um_points_to_find = vertcat((0:1:cil(i).tt_um(end))', cil(i).tt_um(end));
    cil(i).idx_umm = knnsearch(cil(i).tt_um, um_points_to_find);
    
    cil(i).tr_um = cil(i).tt_um(cil(i).idx_umm);
    cil(i).tr = cil(i).tt(cil(i).idx_umm);
    cil(i).xr = cil(i).xx(cil(i).idx_umm);
    cil(i).yr = cil(i).yy(cil(i).idx_umm);
    
    % timestamp
    
    cil(i).timestamp = clicked_frames(i).timestamp;
    
end




if flag_plot
    figure;
    hold on
    box on
    axis image
    plot(cil(1).xx,cil(1).yy,'b');
    plot(cil(1).xr,cil(1).yr,'c.-');
    plot(cil(2).xx,cil(2).yy,'r');
    plot(cil(2).xr,cil(2).yr,'m.-');
end

end



