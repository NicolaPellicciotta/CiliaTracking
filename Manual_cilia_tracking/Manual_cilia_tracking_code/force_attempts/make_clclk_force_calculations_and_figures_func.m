function [] = make_clclk_force_calculations_and_figures_func(clclk_int_fullfilename)

%% input check

if nargin < 1 || isempty(clclk_int_fullfilename)
    [clclk_int_filename, clclk_int_pathname] = uigetfile('*.clclk_int');
else
    [clclk_int_filename, clclk_int_pathname] = parse_filename(clclk_int_fullfilename); % this adds the path if a local file was given as an input
end

clclk_int_fullfilename = fullfile(clclk_int_pathname, clclk_int_filename);


%% prepare savenames

[savepath,clclk_int_basename] = fileparts(clclk_int_fullfilename);
displ_plot_fullsavename = fullfile(savepath,clclk_int_basename);
forcealongcilium_fig_fullsavename = [displ_plot_fullsavename,'_force_along_cilium'];

%% analyse and calculate force

clclk = clclk_calculate_force( clclk_int_fullfilename );

% plot results and save
hf = displ_plot(clclk,false);
hf1 = displ_plot(clclk,true);
% print2svg(hf,displ_plot_fullsavename,1,1);

% new plot of force along cilium and save
hf2 = plot_force_along_cilium(clclk, true);
hf3 = plot_force_along_cilium(clclk, false);
% print2svg(hf2,forcealongcilium_fig_fullsavename,1,1);


% save also file
clclk.clclk_force_fullfilename = fullfile(savepath, [clclk_int_basename '.clclk_force']);
save(clclk.clclk_force_fullfilename,'-struct','clclk')



end %function




