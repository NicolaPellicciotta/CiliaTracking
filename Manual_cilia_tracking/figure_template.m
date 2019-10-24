function [hf, ha] = figure_template()
%figure_template Silly funciton to create 8 by 6 cm figure
%   Detailed explanation goes here


markarr = '^odsv';
cmap = [27,158,119 ;
    217,95,2 ;
    117,112,179 ] ./ 255;

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;   % line width
mkrsz = 6;
hapos = [0.13 0.155 0.82 0.79];


hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes('Position',hapos);
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';
ha.FontSize = tfs;

ha.XLabel.FontSize = lfs;
ha.YLabel.FontSize = lfs;

ha.XLabel.Units = 'normalized';
ha.XLabel.VerticalAlignment = 'baseline';
ha.XLabel.Position(2) = -0.18;
ha.YLabel.Units = 'normalized';
ha.YLabel.VerticalAlignment = 'cap';
ha.YLabel.Position(1) = -0.15;