%%%----- this is the workflow developed by Dr Luigi Feriani to ---------
%%%------- track the cilium waveform and find the force with RFT ------
%%%------ for details read the paper by E.Hamilton 2019


cilia_clicky_thingy_GUI   %%%% 
% selcet a .movie file . If you have a differnt format, first save it as
% series of tiffs and the use the movie2tiff.c file to convert it.
% parameters:
% clickrate: set the fps at which you want to click (if it is lower than the original one you are skipping some frames )
% with 3 and 4 you move to the previous or next frame.
% with alternative view you apply a background subtractiojn to the frame 
% press click point and click along the cilium (also the base). To finish
% click the last point with the right bottom of the mouse, and then you can
% move to the next frame. When you are done and you did at least one beating cycle
% just press save and exit. 
%this one at the end save a file .clclk.
%%% to explore the content of the file is  clclk = load('path/to/file.clclk', '-mat')

%% step 2
%%------------- step 2 ------------------------------------
%%----------- find surface and other staff -----------------

filename= 'chlamy_3000fps';
cilia_out=clclk_find_intersections(strcat(filename,'.clclk'));
%%%% filename is the name of the video without the extension .movie
%%%% close the first window and the second that ask you to open a file
%%%% then follow instruction.
%%%% there is some problems with the window figure 401, this one you need
%%%% to select the surface of the epithelium/organism and it must be a
%%%% closed polygon.
%%%% always insert the right pixel to micron ratio when asked.
%%%% if you don t need dewobbling just close the window that ask so.
%%%% if everything is good you gonna end up with a plot that ask you abou
%%%% the direction of the power and recovery stroke. then you must save
%%%% your varialble in this way:

save(strcat(filename,'.clclk_int'),'-struct','cilia_out');

%% step 3
%%%---------- Force calcultaion ---------------------
close all
make_clclk_force_calculations_and_figures_func(strcat(filename,'.clclk_int'));
% this first ask you if you want to remove some bad tracking, if not just
% close the window
% then you get a lot of plot and  should save a .clclk_force

% to load results from force calculations 

force = load(strcat(filename,'.clclk_force'), '-mat');

%% to do a plot of the trajectories and how manipulate the data

figure();
colors= jet(40);
for jj=1:40;
    color=colors(jj,:);
    plot(force.px2mum*force.points(jj).cilium_x,-force.px2mum*force.points(jj).cilium_y,'color',color);hold on;
end
colorbar()
axis equal

%%% useful data
force.px2mum; %%%% pixel to micron ratio
force.points(jj).cilium_x; %x-coordinate of the timepoint jj
force.points(jj).cilium_y; %y-coordinate of the timepoint jj (this one are inverted, indeed I put a - in the plot figure)
force.clicked_frames(jj).timestamp;   %%% timepoint jj in seconds
force.approximate_clickrate; %%%% frame rate of the clicking

%%% there are some empty timepoints so maybe just check for them
%%% other staff 
%%% force.Stroke.mean_cilium_length
%%%   force.Force.displ.cod_cy                centre of drag x
%%% force.Force.displ.cod_cy                 centre of drag y
%%% force.Force.displ.commoncyls_F_y_pN      force y
%%% force.Force.displ.commoncyls_F_x_pN      force x
%%%  force.Force.displ.commoncyls_F_para_pN  force para to the cell, to the
%%%  surface of the cell
%%%  force.Force.displ.commoncyls_F_perp_pN  force perp to the cell  
