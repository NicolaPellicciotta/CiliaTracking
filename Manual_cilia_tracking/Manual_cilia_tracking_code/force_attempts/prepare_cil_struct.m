function [ cil ] = prepare_cil_struct( clicked_frames, points, frames_numbers_to_analyse, px2mum )
%prepare_cil_struct prepares the 2-cilia-structure to find the best offset to
%compare cilia and eventually calculate the force
%   Detailed explanation goes here

%% input check

if numel(frames_numbers_to_analyse) ~= 2
    error('frames_numbers_to_analyse has to be a 2-elements vector');
end

frames_numbers_to_analyse = sort(frames_numbers_to_analyse);


%% flags

flag_plot = 0;


%% build structure

for i = 2:-1:1
    
    % x and y interpolated coordinates
    cil(i).xx = points(frames_numbers_to_analyse(i)).dw_cilium_xx;
    cil(i).yy = points(frames_numbers_to_analyse(i)).dw_cilium_yy;
    
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
    
    cil(i).timestamp = clicked_frames(frames_numbers_to_analyse(i)).timestamp;
    
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

