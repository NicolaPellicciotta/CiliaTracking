%%%%%% work flow for cilia manual tracking from luigi
data_dir='/home/np451/Desktop/sync_paper/Karl Ferdinand mouse brain cilia/atAno5_16000/formovie';
cd(data_dir);
filename='wt_';
filename='atAno_';
px2mu= 0.107 %%%% for wt

px2mu= 0.0543;
px2mu = 0.0452
cilia_clicky_thingy_GUI   %%%% this one at the end save a file .clclk
%%% to explore the content of the file is  clclk = load('path/to/file.clclk', '-mat')

cilia_out=clclk_find_intersections(strcat(filename,'.clclk'));

%%%%% press ok and the cancel when ask for a mask and for a file, after it
%%%%% plot and ask for rotation, wait 10 seconds otherwise the code crash
%%%%% not sure why px2mu for 60X 0.0933

save(strcat(filename,'clclk_int'),'-struct','cilia_out');


ciao=make_clclk_force_calculations_and_figures_func(strcat(filename,'clclk_int'));
%%%% this one should save a .clclk_force
%%% to load results from force calculations 
%60X_2.14Dec2018_15.20.clclk_force
%force = load(strcat(filename,'.clclk_force'),'-mat');
%force = load('60X_2.14Dec2018_15.20.clclk_force','-mat');
force = load('atAno_clclk_int.clclk_force','-mat');
%%% force.Stroke.mean_cilium_length
%%%   force.Force.displ.cod_cy centre of drag x
%%% force.Force.displ.cod_cy centre of drag y
%%% force.Force.displ.commoncyls_F_y_pN   force y
%%% force.Force.displ.commoncyls_F_x_pN   force x
%%%  force.Force.displ.commoncyls_F_para_pN  force para to the cell, to the
%%%  surface of the cell
%%%  force.Force.displ.commoncyls_F_perp_pN  force perp to the cell  



%% boring things to create a movie ifle from tiffs

%%%%% write the name in the right format that I ate

d=dir('*.tif');
imagename= 'quadri_trot_3000fps';
for i=1:numel(d)
    filename= d(i).name;
    number = str2num(filename(end-7:end-4));
    newname=strcat(imagename,num2str(number),'.tif');
    movefile( filename, newname);
end