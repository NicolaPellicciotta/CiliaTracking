function [clclk] = clclk_find_intersections(clclk_filename, flag_debugging, flag_make_video)

%% input parsing

if nargin < 3 || isempty(flag_make_video)
    flag_make_video = false;
end %if

if nargin < 2 || isempty(flag_debugging)
    flag_debugging = false;
end %if


if nargin < 1 || isempty(clclk_filename)
    [clclk_filename, clclk_pathname] = uigetfile('*.clclk*');
else
    [clclk_pathname, clclk_filename, ext] = fileparts(clclk_filename);
    clclk_filename = [clclk_filename,ext];
    if isempty(clclk_pathname)
        clclk_pathname = pwd;
    end %if
end %if

clclk_fullfilename = fullfile(clclk_pathname, clclk_filename);
flag_make_video

%% constants

% folder in which to look for the string
target_dir = {'E:\Data\Cilia_Profile',...
    'E:\Data\Cilia_Amelia\2015_12_11\2015_12_11_profile\CC15_48B_88B_PIVlab_sessions'};

% folder in which to look for excel files with magnifications
magtable_dir = 'E:\Data\Cilia_Profile\Magnification_Tables';

% regular expression looking for 2digits, 1Ucase, 2 lcase,
% 4digits_2digits.2digits.2digits, hopefully matching the automatic
% date_time of temika
date_time_expr = '\d{2}[A-Z][a-z]{2}\d{4}_\d\d\.\d\d\.\d\d';
time_expr = '\d\d\.\d\d\.\d\d';

polynomium_grade = 4; %grade of cilium-interpolant polynomium

MeasuringArcLengths_mum = 3:7; %to find points along cilium and look at angles

%% look for the right PIV_session_kymo

% clclk file we want to match with a PIV_session
% clclk_filename = 'A.08Jul2015_15.28.53.clclk';

% date and time string extracted from the filename
date_time_str = regexp(clclk_filename, date_time_expr, 'match', 'once');
time_str = regexp(date_time_str, time_expr, 'match', 'once');

if ~isempty(time_str)
    expr_to_match = [regexprep(time_str,'\.','\\.'),'.*session_kymo\.mat'];
    
    cprintf('*[0 .5 0]','Looking for PIV file with matching date/time... ');
    matches = lookforfile(target_dir,expr_to_match);
else
    matches = [];
end


%% interactive bit to select one file if multiple results

if isempty(matches)
    
    hdlg = warndlg({'No matches found.'; 'Please manually select the right session_kymo.mat file'},'modal');
    uiwait(hdlg);
    [fn, pn] = uigetfile;
    if ~fn && ~pn %if user pressed cancel
        clclk.Intersections = [];
        cprintf('*[1 0 0]',...
            'PIV file not found! \nThe mask will have to be selected manually.\n');
        flag_debugging = true; %force to show plots so we know if stuff went wrong in manual selection of rotation etc
        kymo_filename = [];
    else
        kymo_filename = fullfile(pn, fn);
    end %if
    
end %if

if numel(matches) > 1 && iscellstr(matches)
    
    ok = false;
    while ok == false
        [Selection, ok] = listdlg('ListString',matches,...
            'SelectionMode','single',...
            'Name','Select one file',...
            'PromptString','Select the right match',...
            'ListSize',[6*length(matches{1}), 300]);
    end %while
    kymo_filename = matches{Selection};
    
elseif numel(matches) == 1 && iscellstr(matches)
    
    kymo_filename = matches{1};
    
end %if

cprintf('*[0 .5 0]','Done!\n');


%% load the session_kymo file (if it exists) and clclk file

cprintf('*[0 .5 0]','Loading files...\n')

disp(clclk_fullfilename);
disp(kymo_filename);

% clclk = load(clclk_fullfilename,'-mat','clicked_frames','points');
[clclk] = clclk_reader(clclk_fullfilename, 0, 'blue_gradient', true);

if ~isempty(kymo_filename)
    load(kymo_filename,'CommonProperties');
end %if

cprintf('*[0 .5 0]','Done!\n')


%% just in case the kymo are still with the wrong px2mum find magnification by looking at the Magnification Tables

mag_str = find_magnification(magtable_dir, regexprep(date_time_str,'\.','\\.'));
if ~isempty(mag_str)
    switch mag_str
        case '60X'
            CommonProperties.px2mum = 0.097;
        case '60X_1.5X'
            CommonProperties.px2mum = 0.0647;
        otherwise
            warning('Magnification not found')
    end %switch
end %if

%% if kymo file doesn't exist create all possible properties

if isempty(kymo_filename)
    
    % this will make kymo_frame identical to clclk_frame
    CommonProperties.firstframe = clclk.clicked_frames(min([clclk.clicked_frames.frame_number])).IM;
    %{
    % then ask user if to rotate it
    hf2 = figure(4001);
    imshow(CommonProperties.firstframe,[]);
    ldlg = [];
    
    while isempty(ldlg)
        ldlg = listdlg('ListString', {'none','90 clockwise','90 counterclockwise','180'},...
            'SelectionMode','single','PromptString','Select wanted rotation:');
    end %while
    
    close(hf2);
    switch ldlg
        case 1 %none
        case 2 %90cw
            CommonProperties.firstframe = rot90(CommonProperties.firstframe,-1);
        case 3 %90ccw
            CommonProperties.firstframe = rot90(CommonProperties.firstframe);
        case 4 %180
            CommonProperties.firstframe = rot90(CommonProperties.firstframe,2);
    end %switch
    %}
    % mask selection
    hf2 = figure(4001);
    imshow(CommonProperties.firstframe,[]);
    CommonProperties.mask = roipoly(mat2gray(CommonProperties.firstframe));
    close(hf2)
    
    % measuring lines
    dlgopt.Interpreter = 'Tex';
    answers = inputdlg({'Insert heights of measuring lines, [\mum]'},...
        'Insert parameters',1,...
        {num2str([2.5 3.5 4.5 5.5]);},...
        dlgopt);
    CommonProperties.MeasuringLineHeight_mum = str2num(answers{1});
    
    if ~isfield(CommonProperties,'px2mum')
        answers = inputdlg({'\mum/px ratio (0.0647 for Grasshopper and 60X + 1.5X magnification)'},...
            'Insert parameters',1,...
            {num2str(0.0647)},...
            dlgopt);
        CommonProperties.px2mum = str2double(answers{1});
    end %if
    
end %if


%% find a frame for each case

% kymo
kymo_frame = CommonProperties.firstframe;

% clclk
fcf = min([clclk.clicked_frames.frame_number]);% first clicked frame
clclk_frame = double(clclk.clicked_frames(fcf).IM);
% kymo_frame = clclk.clicked_frames(fcf).IM; %this was for emergency when
% video woulnd't register with kymo_frame 

%% decide rotation/transposition/flipud/fliplr to match kymo conditions
% kymo is the reference one, because had been rotated as the user wanted

cprintf('*[0 .5 0]','Finding affine transformation... ')

if size(kymo_frame) == size(clclk_frame) %either flipud or fliprl or rot180
    
    [hh, ww] = size(kymo_frame);
    
    norm_coeff   = max( max( normxcorr2(clclk_frame,          double(kymo_frame) ) ));
    flipud_coeff = max( max( normxcorr2(flipud(clclk_frame),  double(kymo_frame) ) ));
    fliplr_coeff = max( max( normxcorr2(fliplr(clclk_frame),  double(kymo_frame) ) ));
    rot180_coeff = max( max( normxcorr2(rot90(clclk_frame,2), double(kymo_frame) ) ));
    
    [~,whichtransf] = max([norm_coeff, flipud_coeff, fliplr_coeff, rot180_coeff]);
    
    switch whichtransf
        
        case 1 % identical
            cprintf('*[0 .5 0]','No transformation needed\n');
            TM = eye(3);  % affine transformation matrix
            tfun = @(x)x; % function to rotate the frame itself
            
        case 2 % flipud
            cprintf('*[0 .5 0]','Flip Upside-Down needed\n');
            TM = [1 0 0; 0 -1 hh; 0 0 1];% also shift y
            tfun = @(x)flipud(x);
            
        case 3 % fliplr
            cprintf('*[0 .5 0]','Flip Left-Right needed\n');
            TM = [-1 0 ww; 0 1 0; 0 0 1];% also shift x
            tfun = @(x)fliplr(x);
            
        case 4 % rot 180
            cprintf('*[0 .5 0]','180deg rotation needed\n');
            TM = [-1 0 ww; 0 -1 hh; 0 0 1];% shift both x and y
            tfun = @(x)rot90(x,2);
            
        otherwise
            error('Sth went wrong in detection of transformation');
            
    end %switch
    
elseif size(kymo_frame) == flip(size(clclk_frame)) %either rot90 or rot-90 or transpose
    
    [hh, ww] = size(kymo_frame);
    
    trans_coeff    = max( max( normxcorr2(clclk_frame',          double(kymo_frame) ) ));
    rot90ccw_coeff = max( max( normxcorr2(rot90(clclk_frame),    double(kymo_frame) ) )); %counterclockwise
    rot90cw_coeff  = max( max( normxcorr2(rot90(clclk_frame,-1), double(kymo_frame) ) )); %clockwise
    
    [~,whichtransf] = max([trans_coeff, rot90ccw_coeff, rot90cw_coeff]);
    
    switch whichtransf
        
        case 1 % transpose
            cprintf('*[0 .5 0]','Transposition needed\n');
            TM = [0 -1 hh; -1 0 hh; 0 0 1];
            tfun = @(x)x';
            
        case 2 % rot90 ccw
            cprintf('*[0 .5 0]','Ccw 90deg rotation needed\n');
            TM = [0 1 0; -1 0 hh; 0 0 1]; % a ccw rotation of the figure (YDir reverse) is a cw rotation of the axis of clicking BC fuck MatLab that's why
            tfun = @(x)rot90(x,1);
            
        case 3 % rot90 cw
            cprintf('*[0 .5 0]','Cw 90deg rotation needed\n');
            TM = [0 -1 ww; 1 0 0; 0 0 1]; % a cw rotation of the figure (YDir reverse) is a ccw rotation of the axis of clicking BC fuck MatLab that's why
            tfun = @(x)rot90(x,-1);
            
        otherwise
            error('Sth went wrong in detection of transformation');
            
    end %switch
    
else
    clclk.Intersections = [];
    cprintf('*[1 0 0]','The frames in the clclk and kymo files are of different sizes.\n');
    cprintf('*[1 0 0]','The software, as is, can''t cope with that.\n');
    cprintf('*[1 0 0]','Empty structure returned.\n');
    return
    
end %if


%% now apply transformation to all clicked frames, and use the for loop to sort and interpolate points as well

cprintf('*[0 .5 0]','Appliying transformation to clicked points... ');

% apply transformation to first clicked frame as well
clclk_frame = tfun(clclk_frame);

% line integral function
lineintegralfun = @(x,y)sum( sqrt( diff(x).^2 + diff(y).^2 ) );

idx_clicked_frames = [clclk.clicked_frames.frame_number];

for i = idx_clicked_frames
    
    % apply transformations:
    
    clclk.clicked_frames(i).IM = tfun(clclk.clicked_frames(i).IM);
    
    [clclk.points(i).cilium_x, clclk.points(i).cilium_y] = transf2D(TM,...
        clclk.points(i).cilium_x, clclk.points(i).cilium_y);
    
    [clclk.points(i).cilium_xx, clclk.points(i).cilium_yy] = transf2D(TM,...
        clclk.points(i).cilium_xx, clclk.points(i).cilium_yy);
    
    [clclk.points(i).cilium_sorted_x, clclk.points(i).cilium_sorted_y] = transf2D(TM,...
        clclk.points(i).cilium_sorted_x, clclk.points(i).cilium_sorted_y);
        
    clclk.points(i).cilium_length = lineintegralfun(clclk.points(i).cilium_xx,clclk.points(i).cilium_yy);
    
    %{
    % plots for debugging
    if flag_debugging
        cla
        hold on
        plot(clclk.points(i).cilium_xx,clclk.points(i).cilium_yy);
        colorset = cool(numel(x));
        arrayfun(@(j)plot(x(j),y(j),...
            'o','Color',colorset(j,:),'MarkerFaceColor','auto'),...
            1:numel(x));
        title(i)
        xlabel(['Sum of residuals = ',num2str(Sx.normr+Sy.normr)]);
        % pause
    end %if
    %}
    
end %for



%%%%%%%%%%%%%%%%% end code taken from clclk_reader %%%%%%%%%%%%%%

cprintf('*[0 .5 0]','Done!\n');


%% dewobbling videos that need so

clclk = clclk_dewobble(clclk, kymo_frame);
% from now on I can use the dw_* variables


% debug plot, to check that dewobble didn't break the sorting
figure
hold on
% plot(extr(:,1,1),extr(:,2,1),'or');
arrayfun( @(i)plot(clclk.points(i).dw_cilium_xx(1), clclk.points(i).dw_cilium_yy(1),'go'),idx_clicked_frames )
arrayfun( @(i)plot(clclk.points(i).dw_cilium_xx(end), clclk.points(i).dw_cilium_yy(end),'g^'),idx_clicked_frames )
% pause

%% create the measuring lines, starting from the mask

Line = create_measuring_lines(CommonProperties.mask, CommonProperties.MeasuringLineHeight_mum, CommonProperties.px2mum);


%% now find intersections of each cilium with all lines a first time
% this will give problems for 2 reasons:
% crossing with no pixels in common
% cilium gets out of region of green lines

cprintf('*[0 .5 0]','Looking for intersections... ');


% if flag_debugging
hf1 = figure(401);
clf
ha1 = axes(hf1);
imshow(CommonProperties.firstframe,[],'parent',ha1)
superimpose_a_mask(ha1,CommonProperties.mask,[1 0 0],0.4);
hold on
arrayfun(@(i)plot( Line(i).edge_coord_x, Line(i).edge_coord_y, 'g'  ),1:numel(Line));
arrayfun(@(i)plot( clclk.points(i).dw_cilium_xx(1), clclk.points(i).dw_cilium_yy(1), 'mo'  ),idx_clicked_frames);
arrayfun(@(i)plot( clclk.points(i).dw_cilium_xx(end), clclk.points(i).dw_cilium_yy(end), 'm^'  ),idx_clicked_frames);
% end %if
%%

for i = flip(idx_clicked_frames)
    
    %     if flag_debugging,
    figure(hf1);
    plot(ha1,clclk.points(i).dw_cilium_xx, clclk.points(i).dw_cilium_yy,'b');

    %     end
    
    % average separation between consecutive clicks
    avg_click_sep = mean(hypot( diff(clclk.points(i).dw_cilium_x), diff(clclk.points(i).dw_cilium_y) ));
    
    for lc = 1:numel(Line)
        
        XX = vertcat(ceil(clclk.points(i).dw_cilium_xx),...
            ceil(clclk.points(i).dw_cilium_xx),...
            floor(clclk.points(i).dw_cilium_xx),...
            floor(clclk.points(i).dw_cilium_xx));
        YY = vertcat(ceil(clclk.points(i).dw_cilium_yy),...
            floor(clclk.points(i).dw_cilium_yy),...
            ceil(clclk.points(i).dw_cilium_yy),...
            floor(clclk.points(i).dw_cilium_yy));
        
        [int_xy, int_ind, int_indline] = intersect( [XX,YY],...
            [Line(lc).edge_coord_x, Line(lc).edge_coord_y],...
            'rows','stable');
        
        % bring indices back to the right range after stacking 4 vectors
        int_ind = mod(int_ind-1,length(clclk.points(i).dw_cilium_xx))+1;
        int_indline = mod(int_indline-1,length(Line(lc).edge_coord_x))+1;
        
        % remove multiple occurrences...
        int_ind = unique(int_ind,'stable');
        int_indline = unique(int_indline,'stable');
        
        % check if double intersection:
        if numel(int_xy >= 2)
            
            % find intersections closest and furthest from basal
            min_int_ind = min(int_ind);
            max_int_ind = max(int_ind);
            
            % then measure their distance (euclidean)
            int_max_sep = hypot( diff(clclk.points(i).dw_cilium_xx([min_int_ind, max_int_ind])),...
                diff(clclk.points(i).dw_cilium_yy([min_int_ind, max_int_ind])) );
            
            % if it is more than twice the average clicking separation keep
            % only the lower group of clicks
            if int_max_sep > 2*avg_click_sep
                
                % take only the ind_int that corresponds to points within
                % avg_click_sep from the intersection the closest to the
                % basal body:
                
                % calculate all distances between first intersection and
                % all the others
                int_seps = hypot( (clclk.points(i).dw_cilium_xx(min_int_ind) - clclk.points(i).dw_cilium_xx(int_ind) ) ,...
                    (clclk.points(i).dw_cilium_yy(min_int_ind) - clclk.points(i).dw_cilium_yy(int_ind) ) );
                
                % and take ony the ones close enough
                int_ind = int_ind(int_seps < avg_click_sep);
                
            end %if
            
        end %if
        
        if ~isempty(int_xy) %found intersection 
            
            % take the one the closest to the basal body, whose distance to
            % the measuring line is at most 0.2px more than the minimum 
            % among all cilium-line distances
            [idx_along_cilium, D] = knnsearch([clclk.points(i).dw_cilium_xx(int_ind), clclk.points(i).dw_cilium_yy(int_ind)],...
                [Line(lc).edge_coord_x(int_indline), Line(lc).edge_coord_y(int_indline)]);
            
            clclk.points(i).inter_xx(lc) = clclk.points(i).dw_cilium_xx(int_ind(min(idx_along_cilium(D<(min(D)+0.2)))));
            clclk.points(i).inter_yy(lc) = clclk.points(i).dw_cilium_yy(int_ind(min(idx_along_cilium(D<(min(D)+0.2)))));
            %             clclk.points(i).inter_xx(lc) = clclk.points(i).dw_cilium_xx(int_ind(min(idx_along_cilium(D==min(D)))));
            %             clclk.points(i).inter_yy(lc) = clclk.points(i).dw_cilium_yy(int_ind(min(idx_along_cilium(D==min(D)))));
            
            %         if flag_debugging
            figure(hf1);
%             plot(clclk.points(i).dw_cilium_xx(int_ind),clclk.points(i).dw_cilium_yy(int_ind),'c.');
%             plot(clclk.points(i).inter_xx(lc),clclk.points(i).inter_yy(lc),'r.','MarkerSize',8);
%             pause
            %         end %if
            
        else            
            clclk.points(i).inter_xx(lc) = nan;
            clclk.points(i).inter_yy(lc) = nan;
        end %if
        
    end %for
    
end %for

% if flag_debugging
figure(hf1);
arrayfun(@(i)plot( clclk.points(i).inter_xx, clclk.points(i).inter_yy, 'r.'  ),idx_clicked_frames);


shg
% end %if



%% find intersections lines
% redefine the ranges of the edge_coord_x  and y lines, I'll need them
% later on for defining an axis, also, figures look nicer

minx = min(vertcat(clclk.points(idx_clicked_frames).dw_cilium_xx));
maxx = max(vertcat(clclk.points(idx_clicked_frames).dw_cilium_xx));

for lc = numel(Line):-1:1 % sneaky allocation
    
    clclk.Intersections(lc).height_mum = Line(lc).height_mum;
    
    % create polynomial line
    range = Line(lc).edge_coord_x >= minx & Line(lc).edge_coord_x <= maxx;
    clclk.Intersections(lc).line_x = Line(lc).edge_coord_x(range);
    clclk.Intersections(lc).line_y = Line(lc).edge_coord_y(range);
    
end %for

% sometimes these lines can have orientations different from one another. Trying to fix this now

% implicitly use the start of the outmost line as reference.
% flip lines so that start-to-start distances between neighbouring lines
% are smaller than start-to-end distances

for lc = numel(Line)-1 : -1 :1
    
    start_to_start_distance = pdist([clclk.Intersections(lc+1).line_x(1), clclk.Intersections(lc+1).line_y(1);...
        clclk.Intersections(lc).line_x(1), clclk.Intersections(lc).line_y(1)])
    
    start_to_end_distance = pdist([clclk.Intersections(lc+1).line_x(1), clclk.Intersections(lc+1).line_y(1);...
        clclk.Intersections(lc).line_x(end), clclk.Intersections(lc).line_y(end)])
    
    end_to_start_distance = pdist([clclk.Intersections(lc+1).line_x(end), clclk.Intersections(lc+1).line_y(end);...
        clclk.Intersections(lc).line_x(1), clclk.Intersections(lc).line_y(1)])
    
    end_to_end_distance = pdist([clclk.Intersections(lc+1).line_x(end), clclk.Intersections(lc+1).line_y(end);...
        clclk.Intersections(lc).line_x(end), clclk.Intersections(lc).line_y(end)])
    
    if start_to_start_distance + end_to_end_distance > start_to_end_distance + end_to_start_distance
        clclk.Intersections(lc).line_x = flip(clclk.Intersections(lc).line_x);
        clclk.Intersections(lc).line_y = flip(clclk.Intersections(lc).line_y);
        
        new_start_to_start_distance = pdist([clclk.Intersections(lc+1).line_x(1), clclk.Intersections(lc+1).line_y(1);...
            clclk.Intersections(lc).line_x(1), clclk.Intersections(lc).line_y(1)])
        new_start_to_end_distance = pdist([clclk.Intersections(lc+1).line_x(1), clclk.Intersections(lc+1).line_y(1);...
            clclk.Intersections(lc).line_x(end), clclk.Intersections(lc).line_y(end)])
        new_end_to_start_distance = pdist([clclk.Intersections(lc+1).line_x(end), clclk.Intersections(lc+1).line_y(end);...
            clclk.Intersections(lc).line_x(1), clclk.Intersections(lc).line_y(1)])
        new_end_to_end_distance = pdist([clclk.Intersections(lc+1).line_x(end), clclk.Intersections(lc+1).line_y(end);...
            clclk.Intersections(lc).line_x(end), clclk.Intersections(lc).line_y(end)])
    end %if
    
end %for

% debugging plot
for lc = numel(Line):-1:1 % sneaky allocation

    %     if flag_debugging,
    figure(hf1);
    hold on
    hpi(lc) = plot(ha1,clclk.Intersections(lc).line_x, clclk.Intersections(lc).line_y, 'm'); %handle plot intersection
    hpism(lc) = plot(ha1,clclk.Intersections(lc).line_x(1), clclk.Intersections(lc).line_y(1),...
        'mx', 'MarkerFaceColor','m'); %start marker
    hpiem(lc) = plot(ha1,clclk.Intersections(lc).line_x(end), clclk.Intersections(lc).line_y(end),...
        'm^', 'MarkerFaceColor','m'); %end marker
    %     end %if
    
end %for

cprintf('*[0 .5 0]','Done!\n'); % looking for intersections


%% dialog to ask whether sign of line is correct (e.g. if arrow matches power stroke direction)

button = questdlg({'Is the arrow pointing';'in the same direction as';'the power stroke?'},'Choose Line Orientation',...
    'Yes','No','Yes');
switch button
    case 'No'
        
        cprintf('*[1 .5 0]','Changing direction of measuring line to match power stroke direction... ');
        for lc = 1:numel(Line)
            
            %flip direction (arc coordinate) of lines
            clclk.Intersections(lc).line_x = flip(clclk.Intersections(lc).line_x);
            clclk.Intersections(lc).line_y = flip(clclk.Intersections(lc).line_y);
            
            %update plot
            figure(hf1)
            delete(hpism(lc))
            delete(hpiem(lc))
            hpism(lc) = plot(clclk.Intersections(lc).line_x(1), clclk.Intersections(lc).line_y(1),...
                'mx', 'MarkerFaceColor','m'); %start marker
            hpiem(lc) = plot(clclk.Intersections(lc).line_x(end), clclk.Intersections(lc).line_y(end),...
                'm^', 'MarkerFaceColor','m'); %end marker
            
        end %for
        cprintf('*[1 .5 0]','Done!\n');

    otherwise
    cprintf('*[0 .5 0]','Measuring line direction already matching power stroke.\n');
        
end

%% debugging plot

% this can also save frames for a video
if flag_debugging
    
    if flag_make_video
        [path, name, ext] = fileparts(clclk_fullfilename);
        mkdir(fullfile(path,name));
        cprintf('*[1 .7 0]','Saving frames... ');
    end %if
    
    hf2 = figure(402);
    clf
    cc = 0;
    for i = idx_clicked_frames
        
        cc = cc+1;
        
        if exist('hp','var'), delete(hp), end
        if exist('hi','var'), delete(hi), end
        
        frame = my_scalebar_burn(clclk.clicked_frames(i).dw_IM,5,CommonProperties.px2mum,[.95 .05]);
        ts = mat2gray(create_timestamp(clclk.clicked_frames(i).timestamp));
        ts = ts .*double((max(frame(:))- min(frame(:))));
        ts = ts + double(min(frame(:)));
        ts = uint16(ts);
        frame(1:size(ts,1),1:size(ts,2)) = ts;
        hi = imshow(frame,[],'border','tight');
%         hi = imshow(frame,[prctile(frame(:),5), prctile(frame(:),99)],'border','tight');
%         superimpose_a_mask(gcf,CommonProperties.mask,[1 0 0],0.4);
        hold on
        %arrayfun(@(i)plot( clclk.Intersections(i).line_x, clclk.Intersections(i).line_y, 'g'  ),1:numel(Line));
        hp = plot(clclk.points(i).dw_cilium_xx, clclk.points(i).dw_cilium_yy,'b.');
        %plot(clclk.points(i).inter_xx, clclk.points(i).inter_yy,'r.');
        
        drawnow;
        shg
        
        if flag_make_video
            hf2.PaperPositionMode = 'auto';
            hf2.InvertHardcopy = 'off';
            hf2.Color = 'none';
            savename = fullfile(path,name,sprintf('%.4d.png',cc));
            print(gcf,savename,'-dpng')
        end %if
        
    end %for
    
    clear path name ext savename
    
    if flag_make_video
        cprintf('*[1 .7 0]','Done!\n');
    end %if
    
end %if

%% time to calculate abs(velocities) now

% figure;
% hold on

cprintf('*[0 .5 0]','Calculating velocities... ');

for lc = 1:numel(Line)
    
    % put together intersections from single frames' structures
    xx = arrayfun(@(i)clclk.points(i).inter_xx(lc),idx_clicked_frames);
    yy = arrayfun(@(i)clclk.points(i).inter_yy(lc),idx_clicked_frames);
    
    % smooth to fix clicking manual errors
    smxx = smooth(xx)'; %transpose bc smooth puts stuff in column
    smyy = smooth(yy)';
    
    % extract times from structs
    ts = [clclk.clicked_frames(idx_clicked_frames).timestamp];
    
    % thow away nans
    idx_throw_away = isnan(xx) | isnan(yy) | isnan(ts);
    
    % crude calculation of velocity
    absv = hypot(diff(smxx(~idx_throw_away)), diff(smyy(~idx_throw_away)))...
        ./ diff(ts(~idx_throw_away));
    
    % definition of an axis
    DX = clclk.Intersections(lc).line_x(end) - clclk.Intersections(lc).line_x(1);
    DY = clclk.Intersections(lc).line_y(end) - clclk.Intersections(lc).line_y(1);
    
    % define sign of velocity wrt to axis
    signv = sign(sum(vertcat(diff(smxx(~idx_throw_away)),...
        diff(smyy(~idx_throw_away)))...
        .* repmat([DX;DY],1,sum(~idx_throw_away)-1)));
    
    % apply sign to velocity
    vel = absv .* signv;
    
    % write in clclk.Intersections structure
    clclk.Intersections(lc).inter_xx = xx(~idx_throw_away);
    clclk.Intersections(lc).inter_yy = yy(~idx_throw_away);
    clclk.Intersections(lc).inter_smxx = smxx(~idx_throw_away);
    clclk.Intersections(lc).inter_smyy = smyy(~idx_throw_away);
    clclk.Intersections(lc).inter_timestamp = ts(~idx_throw_away);
    clclk.Intersections(lc).inter_v = vel;
    clclk.Intersections(lc).inter_absv = absv;
    clclk.Intersections(lc).inter_signv = signv;
    
    % find average velocity in the 2 directions
    clclk.Intersections(lc).inter_avg_posv = mean(vel(signv > 0));
    clclk.Intersections(lc).inter_avg_negv = mean(vel(signv < 0));
    
    % %     plot(clclk.Intersections(lc).inter_timestamp(2:end),clclk.Intersections(lc).inter_v)
    % %
    % %     plot(clclk.Intersections(lc).inter_timestamp,clclk.Intersections(lc).inter_smxx)
    % %     plot(clclk.Intersections(lc).inter_timestamp,clclk.Intersections(lc).inter_smyy)
    % %     plot(clclk.Intersections(lc).inter_timestamp,clclk.Intersections(lc).inter_xx)
    % %     plot(clclk.Intersections(lc).inter_timestamp,clclk.Intersections(lc).inter_yy)
    %     % pause
    %     subplot(121)
    % hold on
    % plot(intersections_v)
    %     subplot(122)
    %     hold on
    %
    %     plot(hypot(diff(intersections_xx(~idx_throw_away)), diff(intersections_yy(~idx_throw_away)))...
    %         ./ diff(intersections_timestamp(~idx_throw_away)) .* sign(sum(vertcat(diff(intersections_smxx(~idx_throw_away)),...
    %         diff(intersections_smyy(~idx_throw_away)))...
    %         .* repmat([DX;DY],1,sum(~idx_throw_away)-1))))
    
    % find Amplitude as maximum pairwise distance between points on same line
    
    dummy_distances = pdist( [clclk.Intersections(lc).inter_smxx(:), clclk.Intersections(lc).inter_smyy(:)],'Euclidean'  );
    % if we want to know which points are the furthest apart:
    % dummy_distances = squareform(dummy_distances));
    %     [i, j] = ind2sub(size(dummy_distances), find(dummy_distances(:) == max(dummy_distances(:)))); 
    
    clclk.Intersections(lc).inter_amplitude = max(dummy_distances(:));
    
end %for


%% calculate ratio velocities.
%As the sign changes with orientation, we'll assume that power stroke is actually faster than recovery
% I'll use the number of loi

%this would be if power from left to right
try_ratio = [clclk.Intersections.inter_avg_posv] ./ abs([clclk.Intersections.inter_avg_negv]);
% the abs is because one is defined positive

if sum( try_ratio > 1 ) > 2 % then we guessed right
    
    weight = sum( try_ratio > 1 )/4; %how much to trust power/recovery categorisation
    for lc=1:4,
        clclk.Intersections(lc).v_ratio = ...
            clclk.Intersections(lc).inter_avg_posv ./ abs(clclk.Intersections(lc).inter_avg_negv);
        clclk.Intersections(lc).inter_weight = weight;
    end %for
    
elseif sum( try_ratio > 1) < 2 %then we're wrong and power stroke is from right to left
    
    weight = sum( try_ratio < 1 )/4;
    for lc=1:4,
        clclk.Intersections(lc).v_ratio = ...
            abs(clclk.Intersections(lc).inter_avg_negv) ./ clclk.Intersections(lc).inter_avg_posv;
        clclk.Intersections(lc).inter_weight = weight;
    end %for
    
else %we don't really know so we use the third line to decide
    
    if clclk.Intersections(3).inter_avg_posv > abs(clclk.Intersections(3).inter_avg_negv)
        for lc=1:4,
            clclk.Intersections(lc).v_ratio = ...
                clclk.Intersections(lc).inter_avg_posv ./ abs(clclk.Intersections(lc).inter_avg_negv);
        end
    else
        for lc=1:4,
            clclk.Intersections(lc).v_ratio = ...
                abs(clclk.Intersections(lc).inter_avg_negv) ./ clclk.Intersections(lc).inter_avg_posv;
        end %for
    end %if
end %if

cprintf('*[0 .5 0]','Done!\n'); %calculating velocities


%% generate a Stroke structure with average/general properties

% calculate average basal coordinates
clclk.Stroke.basal_avg_xx = mean(arrayfun(@(i)clclk.points(i).dw_cilium_xx(1),idx_clicked_frames));
clclk.Stroke.basal_avg_yy = mean(arrayfun(@(i)clclk.points(i).dw_cilium_yy(1),idx_clicked_frames));

% calculate convex hull (and area) of entire stroke
XX = vertcat(clclk.points(idx_clicked_frames).dw_cilium_xx); %cat all vectors together
YY = vertcat(clclk.points(idx_clicked_frames).dw_cilium_yy);
[K,V] = convhull(XX,YY);
clclk.Stroke.ConvHull_xx = XX(K);
clclk.Stroke.ConvHull_yy = YY(K);
clclk.Stroke.ConvHull_area = V;

% calculate average cilium length
clclk.Stroke.mean_cilium_length = mean([clclk.points(idx_clicked_frames).cilium_length]);
clclk.Stroke.std_cilium_length  = std([clclk.points(idx_clicked_frames).cilium_length]);


%% debugging plot

if flag_debugging
    hf3 = figure(403);
    clf
    
    convhull_mask = poly2mask(clclk.Stroke.ConvHull_xx,clclk.Stroke.ConvHull_yy,size(clclk_frame,1), size(clclk_frame,2));
    imshow(clclk_frame,[]);
    superimpose_a_mask(hf3,convhull_mask,[0 1 1],0.3);
    hold on
    plot(clclk.Stroke.ConvHull_xx, clclk.Stroke.ConvHull_yy,'c');
    hp = plot(clclk.Stroke.basal_avg_xx, clclk.Stroke.basal_avg_yy,'om');
    hp.MarkerFaceColor = hp.Color;
    arrayfun(@(i)plot(clclk.points(i).dw_cilium_xx,clclk.points(i).dw_cilium_yy,'b','LineWidth',1.2),idx_clicked_frames);
    
end %if


%% measure angles from average basal position to points defined at fixed arclength of cilium

%cumulative arclength (coordinata curvilinea)
arclengthfun = @(x,y)cumsum(sqrt( [0;diff(x)].^2 + [0;diff(y)].^2 ));

MeasuringArcLengths_px = MeasuringArcLengths_mum / CommonProperties.px2mum; %should be row

% if flag_debugging
%     hf = figure(404);
%     clf;
%     imshow(clclk_frame,[]);
%     hold on
%
% end %if

ColorSet = spring(numel(MeasuringArcLengths_px));

for i = idx_clicked_frames
    
    %define curvilinear coordinate along cilium
    tt = arclengthfun(clclk.points(i).dw_cilium_xx, clclk.points(i).dw_cilium_yy);
    
    % find where such coordinate is closer to the measuring values
    [~,idx_arclength] = min(abs( repmat(tt,1,numel(MeasuringArcLengths_px)) - repmat(MeasuringArcLengths_px, numel(tt), 1) ));
    
    if any(idx_arclength == length(tt));
        cprintf('*[1 0 0]','Cilium too short for this.\n');
    end
    
    % find angles of lines joining average basal position with points along
    % cilium
    dx = clclk.points(i).dw_cilium_xx(idx_arclength) - clclk.Stroke.basal_avg_xx;
    dy = clclk.points(i).dw_cilium_yy(idx_arclength) - clclk.Stroke.basal_avg_yy;
    clclk.points(i).angles = atan2d(dy,dx);
    
    % enforce the warning given before: if throw away measurements at
    % the tip of the cilium
    clclk.points(i).angles(idx_arclength == length(tt)) = nan;
    
    
    if flag_debugging,
        figure(hf3);
        scatter(clclk.points(i).dw_cilium_xx(idx_arclength),clclk.points(i).dw_cilium_yy(idx_arclength),[],ColorSet,'.');
        
    end %if
end %for

drawnow
shg
pause(1)
        
% find angular amplitude at the measuring lengths
angles_mat = [clclk.points(idx_clicked_frames).angles]; %numel(MeasuringArcLengths) by number of clicked frames
clclk.Stroke.angular_amplitude_deg = max(angles_mat,[],2) - min(angles_mat,[],2); %col vector, numel(MeasuringArcLengths) long
clclk.Stroke.MeasuringArcLengths_mum = MeasuringArcLengths_mum'; %col vector, as the angular amplitude


%% transform in microns

px2mum = CommonProperties.px2mum;

for lc = 1:4
    
    clclk.Intersections(lc).inter_v           = clclk.Intersections(lc).inter_v         .* px2mum;
    clclk.Intersections(lc).inter_absv        = clclk.Intersections(lc).inter_absv      .* px2mum;
    clclk.Intersections(lc).inter_avg_posv    = clclk.Intersections(lc).inter_avg_posv  .* px2mum;
    clclk.Intersections(lc).inter_avg_negv    = clclk.Intersections(lc).inter_avg_negv  .* px2mum;
    clclk.Intersections(lc).inter_amplitude   = clclk.Intersections(lc).inter_amplitude .* px2mum;
    
end %for

clclk.Stroke.mean_cilium_length = clclk.Stroke.mean_cilium_length * px2mum;
clclk.Stroke.std_cilium_length  = clclk.Stroke.std_cilium_length  * px2mum;
clclk.Stroke.ConvHull_area      = clclk.Stroke.ConvHull_area      * px2mum * px2mum;
clclk.Stroke.px2mum             = CommonProperties.px2mum;

end %function




%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%



function Line = create_measuring_lines(mask, heights_mum, px2mum)

%% initialisation of Line structure

N_lines = length(heights_mum);

Line.height_mum   = [];
Line.height_px    = [];
Line.edge_coord_x = [];
Line.edge_coord_y = [];

Line = Line(ones(N_lines,1)); %this makes it a struct array


for lc = 1:N_lines
    
    %dilate mask to find a line at a certain distance from initial mask
    Line(lc).height_mum = heights_mum(lc);
    Line(lc).height_px = round(heights_mum(lc) / px2mum);
    dilated_mask = imdilate(mask,strel('disk',Line(lc).height_px,0));
    
    %find the edge of the dilated mask
    edge_dilated_mask = bwmorph(edge(dilated_mask),'skel');                 %1 at the coordinates where we will want the values of the velocities
    
    %     %plot stuff
    %     superimpose_a_mask( hf, edge_dilated_mask, [0 1 0], 1);  %add edge of dilated mask
    %     drawnow;
    
    %find coordinates of the line on which we'll want the velocities
    edge_coord = bwboundaries(edge_dilated_mask); %output is a cell array
    Line(lc).edge_coord_x = edge_coord{1}(:,2);  % columns
    Line(lc).edge_coord_y = edge_coord{1}(:,1);  % rows
    
    %using bwboundaries works well but they are all duplicated:
    %trying to get rid of duplication issue
    ind = sub2ind(size(mask), Line(lc).edge_coord_y, Line(lc).edge_coord_x); %easier to spot duplicates with index
    ind = unique(ind,'stable');
    N_unique_points = numel(ind);   %will be useful later
    [Line(lc).edge_coord_y, Line(lc).edge_coord_x] = ind2sub(size(mask), ind);  %going back to row, columns
    
    %now points are unique but not in the order a human would draw a continuous line:
    consecutive_points_distance = sqrt(diff(Line(lc).edge_coord_x).^2+diff(Line(lc).edge_coord_y).^2);    %calculating the distance between points consecutive in the edge_coord vectors, to see where there are discontinuities
    breaking_points_ind = find(consecutive_points_distance > 2)+1; %actually should be enough >sqrt(2) to find discontinuities
    
    %label segments that are continuous
    segment_label = zeros(N_unique_points,1);
    segment_label(breaking_points_ind) = true;
    segment_label = cumsum(segment_label)+1;
    order_of_points = (1:N_unique_points)'; %store initial order of points
    
    %now reshuffle segments to minimise "line integral" with matlab ordering of
    %points: try all possible permutations and see which one's the one that
    %yields inimum "line integral"
    possible_permutations = perms(1:max(segment_label));    %this is a matrix, on each line a possile permutation of the segment label
    minvalue = Inf;
    minpermindex = 0;   %index of permutation with minimum consecutive_points_distance
    for pc = 1:size(possible_permutations,1)    %for each permutation
        dummy_order_of_points = [];
        for sc = possible_permutations(pc,:)           %create a dummy vector with order of segments specified by permutation
            dummy_order_of_points = vertcat(dummy_order_of_points,order_of_points(segment_label == sc));
        end
        %standard code for finding min
        total_consecutive_points_distance = sum(sqrt(diff(Line(lc).edge_coord_x(dummy_order_of_points)).^2+diff(Line(lc).edge_coord_y(dummy_order_of_points)).^2));
        if total_consecutive_points_distance < minvalue
            minvalue = total_consecutive_points_distance;
            minpermindex = pc;
        end
    end
    
    %here we should know the right permutation of segment labels that yields
    %minimum consecutive_points_distance, so we write it in "definitve"
    %variables
    dummy_order_of_points = [];
    dummy_segment_label = [];
    for sc = possible_permutations(minpermindex,:)           %create a dummy vector with order of segments specified by permutation index found to be the best
        dummy_order_of_points = vertcat(dummy_order_of_points,order_of_points(segment_label == sc));
        dummy_segment_label = vertcat(dummy_segment_label, segment_label(segment_label==sc));
    end
    order_of_points = dummy_order_of_points;    %store order in definitive variable
    segment_label = dummy_segment_label;        %store labels in definitive variable
    
    %now flip each segment that started with an odd label (empirically works)
    order_of_segment_label = unique(segment_label,'stable');
    for i=1:numel(order_of_segment_label)
        if mod(i,2)
            subset_of_order = order_of_points(segment_label==order_of_segment_label(i));
            order_of_points(segment_label==order_of_segment_label(i)) = flipud(subset_of_order);
        end
    end
    
    %finally, rearrange the coordinates according to the new order
    Line(lc).edge_coord_x = Line(lc).edge_coord_x(order_of_points);
    Line(lc).edge_coord_y = Line(lc).edge_coord_y(order_of_points);
    
end


end %function


function [mag_str] = find_magnification(magtable_dir,datetime_expr_to_match)

% list excel files in target folder
filelist = dir(fullfile(magtable_dir,'*.xlsx'));

% initialise cells
movienames = cell('');
magnifications = cell('');

for fc = 1:numel(filelist)
    
    %read xls file
    [~,xlsdata] = xlsread(fullfile(magtable_dir,filelist(fc).name));
    
    %first column is movie names, third proper magnifications
    movienames = vertcat(movienames,xlsdata{:,1});
    magnifications = vertcat(magnifications,xlsdata{:,3});
    
end %for

%tries and match movienames with string
matches = regexp(movienames,datetime_expr_to_match,'match');


if sum(~cellfun(@isempty,matches)) > 1 %hopefully will never happen because there'll be an error
    
    %find nonempty match
    idx_match = find(~cellfun(@isempty,matches));
    
    ok = false;
    while ok == false
        [Selection, ok] = listdlg('ListString',vertcat(matches{idx_match}),...
            'SelectionMode','single',...
            'Name','Select one file',...
            'PromptString','Select the right match',...
            'ListSize',[20*length(matches{idx_match(1)}{:}), 30*numel(idx_match)]);
    end %while
    
    idx_match = idx_match(Selection);
    
elseif sum(~cellfun(@isempty,matches)) == 1 %one match only
    
    %find nonempty match
    idx_match = find(~cellfun(@isempty,matches));
    
elseif sum(~cellfun(@isempty,matches)) == 0
    mag_str = [];
    return
end %if


%plot name of the video found
[~,matchingmoviename,ext] = fileparts(movienames{idx_match});
matchingmoviename = [matchingmoviename,ext];
cprintf('*[0 .5 0]',['Video found: ',matchingmoviename,'\n']);

%find magnification string
mag_str = magnifications{idx_match};
cprintf('*[0 .5 0]',['Magnification found: ',mag_str,'\n']);


end %function










