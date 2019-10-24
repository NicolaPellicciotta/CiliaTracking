function [unique_filelist] = compare_clicking(searchstring, targetdir)

% input check
if nargin < 2 || isempty(targetdir)
    targetdir = pwd;
end

if nargin < 1 || isempty(searchstring)
    
    % manual files selection
    successful_selection = false;
    while ~successful_selection
        [filenames, targetdir] = uigetfile({'*.*';'*.clclk';'*.clclk_int';'*.clclk_force';'*.clclk_curv'},...
            'Select files to compare',...
            'MultiSelect','on');
        if ~iscell(filenames)
            h = warndlg('Please select at least two files');
            uiwait(h);
        else
            successful_selection = true;
        end
    end %while
    
    % now use it to create filelist
    for fnc = numel(filenames):-1:1 %filename counter
        filelist(fnc) = dir(fullfile(targetdir, filenames{fnc}));
    end
    
else
    
    filelist = dir(fullfile(targetdir,searchstring));
    
end

disp('looking for duplicates within the following directory:')
disp(targetdir)


[duplicate_matrix, duplicate_array, unique_filelist] = compare_files(filelist);

if ~any(duplicate_array)
    disp('Congratulations, no duplicates')
else
    duplicate_array
end



end %function


function [duplicate_matrix, duplicate_array, unique_filelist] = compare_files(filelist)
%compare_clclks compares clclks by plotting the raw clicking points
%returns a duplicate_matrix in which entry (i,j) is true if clclk(i) ==
%clclk(j), false otherwise

figure;
hold on

cmap = lines(numel(filelist));
markarr = '^odsvp><h*+x.';

for cc = 1:numel(filelist)
    
    clclk = load(fullfile(filelist(cc).folder, filelist(cc).name),'-mat');
    
    % find clicked frames
    ind_clicked_frames = [clclk.clicked_frames.frame_number];
    
    % concatenate all x and y data
    x_plot = cat(1, clclk.points(ind_clicked_frames).cilium_x);
    y_plot = cat(1, clclk.points(ind_clicked_frames).cilium_y);
    
    % plot
    hp(cc) = plot(x_plot, y_plot, markarr(cc),...
        'color',cmap(cc,:),...
        'markersize',6);
    hp(cc).DisplayName = [num2str(cc),' --- ',clclk.savename] ;
    
    % save quantities for automatic recognition of duplicate
    S(cc).ind_clicked_frames = ind_clicked_frames;
    S(cc).x_plot = x_plot;
    S(cc).y_plot = y_plot;
    S(cc).savename = clclk.savename;
    
    clear clclk;
    
end % for


% for loop that identifies if any clclk is a duplicate
duplicate_matrix = false(numel(S));
duplicate_array = zeros(numel(S),1);

for c = 1:numel(S)
    for r = 1:numel(S)
        
        if numel(S(c).ind_clicked_frames) == numel(S(r).ind_clicked_frames)
            if numel(S(c).x_plot) == numel(S(r).x_plot) && numel(S(c).y_plot) == numel(S(r).y_plot)
                
                if all(S(c).x_plot == S(r).x_plot) && all(S(c).y_plot == S(r).y_plot)
                    duplicate_matrix(r,c) = true;
                    if r > c
                        disp(['File ', num2str(c) ,' == file ', num2str(r)])
                        duplicate_array(r) = c;
                    end
                end
            end
        end
        
    end % for r
end % for c

idx_unique_files = duplicate_array == false;
unique_filelist = filelist(idx_unique_files);

disp('savename of unique files');
for i = find(idx_unique_files(:)')
    fprintf('\n\tfile name:\t%s\t\tsavename:\t%s',filelist(i).name,S(i).savename);
end
fprintf('\n');

end