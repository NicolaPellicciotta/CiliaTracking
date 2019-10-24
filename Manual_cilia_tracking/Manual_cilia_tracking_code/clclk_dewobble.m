function [ clclk, dewobbled_stack, timestamp ] = clclk_dewobble( clclk_or_filename, first_frame )
%dewobble_movie uses image registration to return a frame_stack of
%dewobbled frames. It should get rid of rigid transformations
%(translations + rotations)

%% input check

if (nargin < 2 || isempty(first_frame))
    first_frame = []; %if it isn;t empty is because I'm passing it as an input from the kymo_mat session
end

if (nargin < 1 || isempty(clclk_or_filename))
    [filename, filepath] = uigetfile;
    fullfilename = fullfile(filepath, filename);
    
    cprintf('*[0 .5 0]','Loading clclk file... ');
    clclk = load(fullfilename,'-mat');
    cprintf('*[0 .5 0]','Done!\n')
end %if

if ~exist('clclk','var')
    
    if ischar(clclk_or_filename)
        
        [ filename, filepath ] = parse_filename( clclk_or_filename );
        
        
        fullfilename = fullfile(filepath, filename);
        
        cprintf('*[0 .5 0]','Loading clclk file... ');
        clclk = load(fullfilename,'-mat');
        cprintf('*[0 .5 0]','Done!\n')
        
    elseif isstruct(clclk_or_filename)
        
        cprintf('*[0 .5 0]','Input was a clclk structure, no need to load a file\n');
        clclk = clclk_or_filename;
        
    end
end %if


%% get movie parameters

idx_clicked_frames = [clclk.clicked_frames.frame_number]; %clicked frames

fs = cat(3,clclk.clicked_frames(idx_clicked_frames).IM); %3d stack
timestamp = cat(3,clclk.clicked_frames(idx_clicked_frames).timestamp); %times of frame acquisition

if isempty(first_frame)
    first_frame = fs(:,:,1);
end %if
[height, width, N_frames] = size(fs);

flag_read_all = true;


%% ask whether to dewobble or not

hsf = scroll_stack(fs);
hsf.Name = 'Dewobbling needed?';
uiwait(hsf);

button = questdlg('Does this video need dewobbling?','Dewobble?','Yes','No','Yes');


%% if dewobbling isn't needed mimic the creation of the dewobbled fields in the structure, and return the outputs

switch button
    case 'No'
        
        cprintf('*[0 .5 0]','Dewobbling not needed, creating structure fields... ');
        
        for fc = 1:N_frames
            
            cc = idx_clicked_frames(fc);
            clclk.clicked_frames(cc).dw_IM = clclk.clicked_frames(cc).IM;
            clclk.clicked_frames(cc).tform = affine2d(eye(3));
            clclk.points(cc).dw_cilium_x = clclk.points(cc).cilium_x;
            clclk.points(cc).dw_cilium_y = clclk.points(cc).cilium_y;
            if isfield(clclk.points,'cilium_xx')
                clclk.points(cc).dw_cilium_xx = clclk.points(cc).cilium_xx;
                clclk.points(cc).dw_cilium_yy = clclk.points(cc).cilium_yy;
            end %if
            clclk.trimming_coordinates.rmin = 1;
            clclk.trimming_coordinates.rmax = size(fs,1);
            clclk.trimming_coordinates.cmin = 1;
            clclk.trimming_coordinates.cmax = size(fs,2);
            if nargout > 1
                dewobbled_stack = fs;
            end %if
            
        end %for
        
        cprintf('*[0 .5 0]','Done!\n')
        return
        
    otherwise
        cprintf('*[0 .5 0]','Dewobbling needed, doing it now.\n');
end %switch


%% setting up a mask

cprintf('*[0 .5 0]','Select ROI for registration... ');

hf = figure(901);
hf.Name = 'Select ROI for registration:';

mask = roipoly(mat2gray(first_frame));
masked_first_frame = frame_preprocessing(first_frame, mask);

imshow(masked_first_frame,[]);
colorbar
shg
Rff = imref2d(size(first_frame));

cprintf('*[0 .5 0]','Done!\n')


%% registration

cprintf('*[0 .5 0]','Registering frames... ');

global dwfs;
dwfs = zeros(height,width,N_frames,'like',first_frame);
dwfs(:,:,1) = first_frame;

[optimizer,metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 150;


for fc = 1:N_frames
    
    cc = idx_clicked_frames(fc); %counter in the clclk struct array
    
    % load next frame
    frame_to_register = fs(:,:,fc);
    
    %apply mask
    processed_frame_to_register = frame_preprocessing(frame_to_register, mask);
    
    tform = imregtform( processed_frame_to_register, masked_first_frame,...
        'rigid',optimizer,metric, 'InitialTransformation', affine2d(eye(3)));
    
    %{
    cla
    imshowpair(imwarp(frame_to_register,tform,'OutputView',Rff), masked_first_frame);
    title(num2str(fc));
    drawnow;
    pause
        %}
        
        dwfs(:,:,fc) = imwarp(frame_to_register,tform,'OutputView',Rff);
        
        clclk.clicked_frames(cc).tform = tform;
        
end

close(hf)

cprintf('*[0 .5 0]','Done!\n')


%% check for misregistrations
% 
% cprintf('*[0 .5 0]','Checking for misregistrations...\n');
% 
% misregistrations = find_misregistrations;
% old_misregistrations = nan(size(misregistrations));
% 
% % while sum(misregistrations) > 0 &&  sum(misregistrations) < old_sum_misregistration %let's register again on the previous (registered) frame
% while sum(misregistrations) > 0 &&  any(misregistrations ~= old_misregistrations) %let's register again on the previous (registered) frame
%     
%     cprintf('*[0 .5 0]',sprintf('Found %d misregistration. Fixing them now...\n',sum(misregistrations)) );
% 
%     for fc = find(misregistrations)'
%         
%         cc = idx_clicked_frames(fc); %counter in the clclk struct array
%         
%         frame_to_register = fs(:,:,fc);
%         
%         %register the bad frames starting off the transformation matrix of the
%         %previous frame
%         processed_frame_to_register = frame_preprocessing(frame_to_register, mask);
%         
%         if fc > 1 % can't read a tform that isn't there
%             previous_frame_tform = clclk.clicked_frames(idx_clicked_frames(fc-1)).tform;
%         else
%             previous_frame_tform = affine2d(eye(3));    
%         end
%         
%         tform = imregtform( processed_frame_to_register, masked_first_frame,'rigid',...
%             optimizer, metric, 'InitialTransformation', previous_frame_tform);
%         
%         dwfs(:,:,fc) = imwarp(frame_to_register,tform,'OutputView',Rff);
%         clclk.clicked_frames(cc).tform = tform;
%         
%     end
%     
%     old_misregistrations = misregistrations;
%     misregistrations = find_misregistrations;
%     
% end
% 
% cprintf('*[0 .5 0]','Done!\n')

%% check for misregistrations (alternative)

cprintf('*[0 .5 0]','Checking for misregistrations...\n');

misregistrations = find_misregistrations;

cprintf('*[0 .5 0]',sprintf('Found %d misregistration. Fixing them now...\n',sum(misregistrations)) );

%go to each misregistration and deal with it separately
for fc = find(misregistrations)'
    
    if fc == 1
        continue %first frame is never a misregistration (hopefully)
    end
    
    cc = idx_clicked_frames(fc); %counter in the clclk struct array
    
    % prepare frame to register it
    frame_to_register = fs(:,:,fc);
    processed_frame_to_register = frame_preprocessing(frame_to_register, mask);
    
    % first try: register the bad frames starting off the transformation
    %matrix of the previous frame
    previous_frame_tform = clclk.clicked_frames(idx_clicked_frames(fc-1)).tform;
    tform = imregtform( processed_frame_to_register, masked_first_frame,'rigid',...
        optimizer, metric, 'InitialTransformation', previous_frame_tform);
    % apply the attempt
    dwfs(:,:,fc) = imwarp(frame_to_register,tform,'OutputView',Rff);
    clclk.clicked_frames(cc).tform = tform;
    
    % check whether it worked
    new_misregistrations = find_misregistrations;
    if ~any(fc == find(new_misregistrations))
        continue %fc is not classified as a misregistration anymore, so go to the next one
    end
    
    % we get here only if it didn't work
    
    % second try: register the bad frames starting off the transformation
    % matrix of the following frame
    
    if fc < N_frames &&... %because we have to read frame fc+1,
            ~any(fc+1 == find(new_misregistrations)) % also we don;t want to try to start with a wrong tform
        
        
        following_frame_tform = clclk.clicked_frames(idx_clicked_frames(fc+1)).tform;
        tform = imregtform( processed_frame_to_register, masked_first_frame,'rigid',...
            optimizer, metric, 'InitialTransformation', following_frame_tform);
        % apply the attempt
        dwfs(:,:,fc) = imwarp(frame_to_register,tform,'OutputView',Rff);
        clclk.clicked_frames(cc).tform = tform;
        
    end %if
    
    % check whether it worked
    new_misregistrations = find_misregistrations;
    if ~any(fc == find(new_misregistrations))
        continue %fc is not classified as a misregistration anymore, so go to the next one
    end
    
    % third try: start off the identity matrix
    initial_tform = affine2d(eye(3));
    tform = imregtform( processed_frame_to_register, masked_first_frame,'rigid',...
        optimizer, metric, 'InitialTransformation', initial_tform);
    % apply the attempt
    dwfs(:,:,fc) = imwarp(frame_to_register,tform,'OutputView',Rff);
    clclk.clicked_frames(cc).tform = tform;
    
    % check whether it worked
    new_misregistrations = find_misregistrations;
    if ~any(fc == find(new_misregistrations))
        continue %fc is not classified as a misregistration anymore, so go to the next one
    end
    
    
    % last attempt: start off an interpolation of previous and following
    % frame
    if fc < N_frames &&... %because we have to read frame fc+1,
            ~any(fc+1 == find(new_misregistrations)) % also we don;t want to try to start with a wrong tform
        
        prev_T = previous_frame_tform.T;
        foll_T = following_frame_tform.T;
        
        % translation is an average
        try_translation = (prev_T(3,[1 2]) + foll_T(3,[1 2])) ./2;
        
        % angle is average as well
        prev_angle = atan2(prev_T(2,1), prev_T(1,1));
        foll_angle = atan2(foll_T(2,1), foll_T(1,1));
        try_angle  = (prev_angle + foll_angle) /2;
        
        % set the tform
        try_T = [cos(try_angle) -sin(try_angle), 0; sin(try_angle) cos(try_angle) 0; try_translation 1];
        initial_tform = affine2d(try_T);
        
        tform = imregtform( processed_frame_to_register, masked_first_frame,'rigid',...
            optimizer, metric, 'InitialTransformation', initial_tform);
        % apply the attempt
        dwfs(:,:,fc) = imwarp(frame_to_register,tform,'OutputView',Rff);
        clclk.clicked_frames(cc).tform = tform;
        
        % check whether it worked
        new_misregistrations = find_misregistrations;
        if ~any(fc == find(new_misregistrations))
            continue %fc is not classified as a misregistration anymore, so go to the next one
        end
        
        % giving up, using the interpolated tform as final one
        dwfs(:,:,fc) = imwarp(frame_to_register,initial_tform,'OutputView',Rff);
        clclk.clicked_frames(cc).tform = initial_tform;
        cprintf('*[.7 .4 0]','Giving up, using an interpolation of the neighbouring frames'' transformation matrices. ')
        
    end %if
    
end

cprintf('*[0 .5 0]','Done!\n')




%% clear memory from wobbly stack

if flag_read_all
    clear global fs;
end


%% apply registration to frames and clicks in clclk

cprintf('*[0 .5 0]','Applying registration to clclk structure... ');

for fc = 1:N_frames
    
    cc = idx_clicked_frames(fc);
    
    % write dewobbled frame in clclk structure
    clclk.clicked_frames(cc).dw_IM = dwfs(:,:,fc);
    
    % apply transform to cilia coordinates
    [clclk.points(cc).dw_cilium_x, clclk.points(cc).dw_cilium_y] = ...
        transformPointsForward(clclk.clicked_frames(cc).tform,...
        clclk.points(cc).cilium_x,...
        clclk.points(cc).cilium_y);
    
    % and to interpolated ones if they exist
    if isfield(clclk.points,'cilium_xx')
        [clclk.points(cc).dw_cilium_xx, clclk.points(cc).dw_cilium_yy] = ...
            transformPointsForward(clclk.clicked_frames(cc).tform,...
            clclk.points(cc).cilium_xx,...
            clclk.points(cc).cilium_yy);
    end %if
    
end %for

cprintf('*[0 .5 0]','Done!\n')


%% finding the black parts to be able to cut them away for showing if needed

cprintf('*[0 .5 0]','Finding trimming coordinates... ');

mask = any(dwfs == 0,3);

[rmin, rmax, cmin, cmax] = find_coordinates_stack_trimming(mask);

clclk.trimming_coordinates.rmin = rmin;
clclk.trimming_coordinates.rmax = rmax;
clclk.trimming_coordinates.cmin = cmin;
clclk.trimming_coordinates.cmax = cmax;

cprintf('*[0 .5 0]','Done!\n')

% then just use this to trim:
% dewobbled_stack = dwfs(rmin:rmax,cmin:cmax,:);


%% prepare output

if nargout > 1
    dewobbled_stack = dwfs;
end %if
clear global dwfs;


%% failsafe exit

if ~exist('timestamp','var') && nargout > 2
    timestamp = [];
end


end


%%
%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%


function [rmin, rmax, cmin, cmax] = find_coordinates_stack_trimming(mask)

rectangle = false(size(mask));
rectangle(2:end-1,2:end-1) = true;

intersection = true;

while intersection
    
    rectangle = imerode(rectangle,ones(3));
    intersection = any(mask(:) & rectangle(:));
    
end %while

found_rectangle = false;

leftnhood   = false(3); leftnhood(2,1:2)   = true;
rightnhood  = false(3); rightnhood(2,2:3)  = true;
topnhood    = false(3); topnhood(1:2,2)    = true;
bottomnhood = false(3); bottomnhood(2:3,2) = true;

while ~found_rectangle
    
    old_rectangle = rectangle;
    
    try_rect = imdilate(rectangle,leftnhood);
    if ~any(mask(:) & try_rect(:))
        rectangle = try_rect;
    end
    
    try_rect = imdilate(rectangle,topnhood);
    if ~any(mask(:) & try_rect(:))
        rectangle = try_rect;
    end
    
    try_rect = imdilate(rectangle,rightnhood);
    if ~any(mask(:) & try_rect(:))
        rectangle = try_rect;
    end
    
    try_rect = imdilate(rectangle,bottomnhood);
    if ~any(mask(:) & try_rect(:))
        rectangle = try_rect;
    end
    
    if old_rectangle == rectangle
        found_rectangle = true;
    end
    
end %while

cmin = find(any(rectangle,1),1,'first');
cmax = find(any(rectangle,1),1,'last' );
rmin = find(any(rectangle,2),1,'first');
rmax = find(any(rectangle,2),1,'last' );


end %function


function [processed_frame] = frame_preprocessing(IM, mask)

if nargin < 2 || isempty(mask)
    mask = true(size(IM));
end

processed_frame = IM;
processed_frame_mask = IM;

processed_frame = imadjust(imgaussfilt(processed_frame,8));
processed_frame_mask = imadjust(imgaussfilt(processed_frame_mask,1));

processed_frame(mask) = processed_frame_mask(mask);


end %function


function [misregistrations] = find_misregistrations()

global dwfs;

N_frames = size(dwfs,3);
re = round(size(dwfs,1)/10);
ce = round(size(dwfs,2)/10);

black_array_top = arrayfun(@(fc)sum( sum(dwfs(1:re,:,fc) == 0 )), (1:N_frames)' );
sm_black_array_top = smooth(black_array_top,'rlowess');
fluct_top = abs(black_array_top - sm_black_array_top);

black_array_bottom = arrayfun(@(fc)sum( sum(dwfs(end-re:end,:,fc) == 0 )), (1:N_frames)' );
sm_black_array_bottom = smooth(black_array_bottom,'rlowess');
fluct_bottom = abs(black_array_bottom - sm_black_array_bottom);

black_array_left = arrayfun(@(fc)sum( sum(dwfs(:,1:ce,fc) == 0 )), (1:N_frames)' );
sm_black_array_left = smooth(black_array_left,'rlowess');
fluct_left = abs(black_array_left - sm_black_array_left);

black_array_right = arrayfun(@(fc)sum( sum(dwfs(:,end-ce:ce,fc) == 0 )), (1:N_frames)' );
sm_black_array_right = smooth(black_array_right,'rlowess');
fluct_right = abs(black_array_right - sm_black_array_right);

flucts = fluct_top + fluct_bottom + fluct_left + fluct_right;
plot(flucts)

% [~, misreg_loc] = findpeaks(misreg,'MinPeakProminence',1.5*std(misreg));
% %this works pretty well except in a couple of cases

trim_flucts = flucts( flucts >= prctile(flucts,10) & flucts <= prctile(flucts,90));
trim_flucts_std = std(trim_flucts);
misreg_loc = find(flucts > 5*trim_flucts_std);

figure
hold on
plot(flucts)
plot(misreg_loc,flucts(misreg_loc),'ro')
drawnow
pause
close

misregistrations = false(N_frames,1);
misregistrations(misreg_loc) = true;

end %function


