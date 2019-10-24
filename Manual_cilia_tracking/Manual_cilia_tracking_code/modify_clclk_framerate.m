function [] = modify_clclk_framerate(folder, framerate)

% input check
if nargin < 1 || isempty(folder)
    folder = uigetdir('','Select folder that contains the .clclk files');
end

if nargin < 2 || isempty(framerate)
    fr_string = inputdlg('Enter framerate (fps)','Enter framerate',1);
    framerate = str2double(fr_string);
end


% look for clclk files, delete ghost files
fl = dir(fullfile(folder,'*.clclk'));
idx_bad = arrayfun(@(f) strcmp(f.name(1:2),'._'),fl);
fl(idx_bad) = [];

% loop on files
for fc = 1:numel(fl)
    
    filename = fullfile(fl(fc).folder, fl(fc).name);
    modify_one_clclk_framerate(filename, framerate);
end

end


function [] = modify_one_clclk_framerate(filename, framerate)
%modify_clclk_framerate Fix problems if a .movie was created and clicked
%with the wrong framerate

% quick input check
[~,~,ext] = fileparts(filename);
if ~strcmp(ext,'.clclk')
    error('wrong filetype')
end

disp(filename)

% save old file just in case
disp('Doing backup')
copyfile(filename,[filename,'_backup']);

% open file
clclk = load(filename, '-mat');

disp(['fixing framerate to ',num2str(framerate),'fps'])
% fix framerate in movie_object
clclk.movie_object.FrameRate = framerate;

% loop on clicked frames
for fc = [clclk.clicked_frames.frame_number]
    clclk.clicked_frames(fc).timestamp = clclk.clicked_frames(fc).frame_number / framerate;
end

% now save
disp('saving...')
save('-v7.3', filename, '-struct', 'clclk');
disp('Done!')
fprintf('\n');
end

