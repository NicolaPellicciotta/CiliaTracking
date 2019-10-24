function [clclk] = clclk_reader(clclkfilename, flag_debugging, colormap_name, flag_force_userIDbasetip)
%clclk_reader is a custom reader for files created with
%cilia_clicky_thingy_GUI. Open the files and plots the interpolated
%position, then saves a good quality figure

%% input check

if nargin < 4 || isempty(flag_force_userIDbasetip)
    flag_force_userIDbasetip = 0;
end

if nargin < 3 || isempty(colormap_name)
    colormap_name = [];
end

if nargin < 2 || isempty(flag_debugging)
    flag_debugging = false;
end

if nargin < 1 || isempty(clclkfilename)
    % if isempty(clclkfilename)
    [clclkfilename, clclkfilepath] = uigetfile('*.clclk','Select a .clclk file.');
else
    [clclkfilepath,clclkfilename,ext] = fileparts(clclkfilename);
    clclkfilename = [clclkfilename,ext];
end

if isempty(clclkfilepath)
    clclkfilepath = pwd;
end %if

clclkfullfilename = fullfile(clclkfilepath, clclkfilename);


%% load clclk file in clclk structure

clclk = load(clclkfullfilename,'-mat');


%% set parameters

idx_clicked_frames = [clclk.clicked_frames(:).frame_number];

firstclickedframe = min(idx_clicked_frames);
lastclickedframe = max(idx_clicked_frames);

polynomium_grade = 4;


% numbers for scalebar

scalebar_length_mum = 2;
% clclk.px2mum = 0.0647;
if ~isfield(clclk,'px2mum') % line introduced for vitos samples in which everything is at px2mum = 0.093 and i put it manually into the clclk
    clclk.px2mum = parse_filename_for_px2mum(clclkfilename);
end
% if not found then look it up in the magnification table
if isempty(clclk.px2mum)
    date_time_expr = '\d{2}[A-Z][a-z]{2}\d{4}_\d\d\.\d\d\.\d\d';
    date_time_str = regexp(clclkfilename,date_time_expr,'match','once');
    date_time_str = regexprep(date_time_str,'\.','\\.');
    mag_str = find_magnification('E:\Data\Cilia_Profile\Magnification_Tables',date_time_str);
    clclk.px2mum = parse_filename_for_px2mum(mag_str);
end %if

% if still not found then enter manually
if isempty(clclk.px2mum)
    clclk.px2mum = input('Please enter pixel to micron ratio: ');
end % if

scalebar_length_px = scalebar_length_mum / clclk.px2mum;
scalebar_dx_pos = -Inf;
scalebar_y_pos = -Inf;


%% find basal points (dunno which order cilia were clicked

for i = idx_clicked_frames
    
    [x,y] = sortbypath_pdist(clclk.points(i).cilium_x,clclk.points(i).cilium_y); %sort points along cilium
    
    clclk.points(i).cilium_sorted_x = x;
    clclk.points(i).cilium_sorted_y = y;
    
end %for


extr = nan( numel(idx_clicked_frames), 2, 2 ); % 1st dim: frame, 2nd dim: x o y, 3rd dim: group 1 or 2

jc = 0;
for i = idx_clicked_frames
    
    jc = jc + 1;

    % points as of now calssified as bases
    extr(jc, 1, 1) = clclk.points(i).cilium_sorted_x(1);
    extr(jc, 2, 1) = clclk.points(i).cilium_sorted_y(1);
    
    % points as of now classified as tips
    extr(jc, 1, 2) = clclk.points(i).cilium_sorted_x(end);
    extr(jc, 2, 2) = clclk.points(i).cilium_sorted_y(end);
    
end %for

% clusterise
kind_1 = kmedoids( extr(:,:,1) ,2, 'Replicates',5); %clusterise the points so far classified as basal
kind_2 = kmedoids( extr(:,:,2) ,2, 'Replicates',5); %independently clusterise the points so far classified as distal

% now either kind_1,2 are identical, or the opposite of one another
if all(kind_1 == kind_2) && ~flag_force_userIDbasetip 
    
    % basal and distal are not consistently named, flip one
    kind_2 = mod(kind_2,2) + 1;
    
elseif all((kind_1-1) == ~(kind_2-1)) && ~flag_force_userIDbasetip 
    
    % this is ok as basal and distal are labelled consistently (each basal
    % is different than each distal)
    
else %unconsistent clustering!
    disp('unconsistent clustering, retrying')
    
    % now try to cluster all points at the same time, as it should be
    % easier. Also, start with one base and tip as medoids
    kind = kmedoids(vertcat(extr(:,:,1), extr(:,:,2)), 2, 'Start', squeeze(extr(1,:,:))' ); 
    
    % split output for consistency
    kind_1 = kind(1:end/2);
    kind_2 = kind(end/2+1:end);
    
    if any(kind_1-1 == kind_2-1) || flag_force_userIDbasetip 
        
        disp('clustering failed again');
        
        % user interface to ask for draggable rectangle
        all_x = vertcat(extr(:,1,1), extr(:,1,2)); %first half is classified as base, second as tip
        all_y = vertcat(extr(:,2,1), extr(:,2,2)); %first half is classified as base, second as tip
        
        figure;
        imshow(clclk.clicked_frames(firstclickedframe).IM,[]);
        title('Drag rectangle around basal points');
        hold on
        plot(all_x, all_y, '.');
        h = imrect;
        rect = wait(h);
        
        % find which points are in the rectangle
        idx_base = false(size(all_x));
        idx_base = isinrange(all_x,rect(1),sum(rect([1 3]))) & isinrange(all_y,rect(2),sum(rect([2 4])));
        
        % now make notation consistent
        kind = 2 - idx_base;
        kind_1 = kind(1:end/2);
        kind_2 = kind(end/2+1:end);
    
    end
    
end

% group points for easier plotting
XY1 = vertcat(extr(kind_1==1,:,1), extr(kind_2==1,:,2));
XY2 = vertcat(extr(kind_1==2,:,1), extr(kind_2==2,:,2));

% calculate convex hull and area, useful to classify top/bottom
CH1 = convhull(XY1(:,1),XY1(:,2));
CHA1 = polyarea(XY1(CH1,1),XY1(CH1,2));
CH2 = convhull(XY2(:,1),XY2(:,2));
CHA2 = polyarea(XY2(CH2,1),XY2(CH2,2));

% decide which one is the basal. In theory, the one with less eccentricity
% and a smaller mean distance to the centroids is the basal group
CXY1 = mean(XY1);
CXY2 = mean(XY2);

meanD1 = mean( sqrt(sum((XY1 - CXY1).^2,2)));
[~,~,latent1]=pca(XY1 - CXY1);

meanD2 = mean( sqrt(sum((XY2 - CXY2).^2,2)));
[~,~,latent2]=pca(XY2 - CXY2);


if (latent1(1)/latent1(2) >= latent2(1)/latent2(2)) && (CHA1 >= CHA2) && (meanD1 >= meanD2) && ~flag_force_userIDbasetip
    
    % first group is definitely distal, switch
    disp('switching');
    tmp = kind_2;
    kind_2 = kind_1;
    kind_1 = tmp;
    clear tmp
    
    % reassign extremes for easy plotting
    XY1 = vertcat(extr(kind_1==1,:,1), extr(kind_2==1,:,2));
    XY2 = vertcat(extr(kind_1==2,:,1), extr(kind_2==2,:,2));
    CXY1 = mean(XY1);
    CXY2 = mean(XY2);
    
elseif (latent1(1)/latent1(2) < latent2(1)/latent2(2)) && (CHA1 < CHA2) && (meanD1 < meanD2) && ~flag_force_userIDbasetip
    % second group is clearly distal, nothing to do i think
    
else %ambiguous ask the user

    % user interface
    figure('Name','click on base of cilium');
    imshow(clclk.clicked_frames(firstclickedframe).IM,[])
    hold on
    plot(XY1(:,1),XY1(:,2),'r.')
    plot(XY2(:,1),XY2(:,2),'g.')
    [xin,yin] = ginputc(1);
    
    % select basal cluster as the one with the centroid closest to the input point
    xyin = [xin,yin];
    [~,ind_base] = min( sqrt(sum((xyin - vertcat(CXY1,CXY2)).^2,2))); 
    
    % and switch if necessary
    if ind_base == 2
        tmp = kind_2;
        kind_2 = kind_1;
        kind_1 = tmp;
        clear tmp
        XY1 = vertcat(extr(kind_1==1,:,1), extr(kind_2==1,:,2));
        XY2 = vertcat(extr(kind_1==2,:,1), extr(kind_2==2,:,2));
        CXY1 = mean(XY1);
        CXY2 = mean(XY2);
    end
end


% plot to check everything worked
if flag_debugging
    
    figure
    imshow(clclk.clicked_frames(firstclickedframe).IM,[])
    hold on
    
    
    plot(extr(:,1,1),extr(:,2,1),'ro')
    plot(extr(:,1,2),extr(:,2,2),'rs')
    
    plot(XY1(:,1),XY1(:,2),'b.')
    plot(XY2(:,1),XY2(:,2),'y.')
    
    plot(CXY1(1), CXY1(2),'bo');
    plot(CXY2(1), CXY2(2),'yo');
    shg
end


% now all the cilia that need to be flipped should have 2 written into kind_1
for i = idx_clicked_frames(kind_1==2)
    
    clclk.points(i).cilium_sorted_x = flip(clclk.points(i).cilium_sorted_x);
    clclk.points(i).cilium_sorted_y = flip(clclk.points(i).cilium_sorted_y);
    
end %for i

figure
imshow(clclk.clicked_frames(firstclickedframe).IM,[],'border','tight')
hold on
for i = idx_clicked_frames
    
    plot(clclk.points(i).cilium_sorted_x(1), clclk.points(i).cilium_sorted_y(1),'b.');
    plot(clclk.points(i).cilium_sorted_x(end), clclk.points(i).cilium_sorted_y(end),'y.');
    
end %for i




%% fit clicked points with polynomials

% first define standard arclength everything is calculated at (common
% across timepoints)
max_t = max( arrayfun(...
    @(i) sum(hypot(diff(clclk.points(i).cilium_sorted_x), diff(clclk.points(i).cilium_sorted_y))), idx_clicked_frames) );
tt = linspace(0, max_t, 500);

for i = idx_clicked_frames

    x = clclk.points(i).cilium_sorted_x;
    y = clclk.points(i).cilium_sorted_y;
    t = cumsum(sqrt([0;diff(x)].^2 + [0;diff(y)].^2));  % define arclength-coordinate
%         tt = cell2mat(arrayfun(@(i)linspace(t(i),t(i+1),30),1:numel(t)-1,'UniformOutput',false)); % this is for plotting smothly
%     tt = linspace(t(1),t(end),30*numel(t)); % this is for plotting smothly
    
    % fit points with polynomial and evaluate it at plotting points
    if length(x) > polynomium_grade
        [Px,Sx] = polyfit(t,x,polynomium_grade);
        [Py,Sy] = polyfit(t,y,polynomium_grade);
    else
        error('Cilium not clicked enough times: we need at the very least 5 clicks per frame (but the more the better)');
    end
    
    clclk.points(i).cilium_Px = Px;
    clclk.points(i).cilium_Py = Py;
    xx = polyval(Px,tt);
    yy = polyval(Py,tt);
    
    % tt is obviously too long for all cilia. chop xx and yy at closest
    % point with sorted x,y
    [~,wh] = min(   hypot(x(end) - xx, y(end)-yy)  );
    
    clclk.points(i).cilium_xx = xx(1:wh)';
    clclk.points(i).cilium_yy = yy(1:wh)';
    clclk.points(i).tt = tt(1:wh)';
    clclk.points(i).curvature = calculate_abscurvature(clclk.points(i).cilium_xx, clclk.points(i).cilium_yy, ...
        clclk.points(i).tt);
    
    clclk.points(i).normr = Sx.normr + Sy.normr;
    
    scalebar_dx_pos = max(scalebar_dx_pos, max(clclk.points(i).cilium_xx));
    scalebar_y_pos = max(scalebar_y_pos, max(clclk.points(i).cilium_yy));
    
    
end % for


%% remove misclickings (same frame clicked twice) by looking at outliers in the residuals from the fit

mm = median([clclk.points(idx_clicked_frames).normr]);
ss = iqr([clclk.points(idx_clicked_frames).normr]);

% idx_clicked_frames([clclk.points(idx_clicked_frames).normr] < mm - 3*ss | [clclk.points(idx_clicked_frames).normr] > mm + 3*ss) = [];


misclicking = ([clclk.points(idx_clicked_frames).normr] < mm - 3*ss |...
    [clclk.points(idx_clicked_frames).normr] > mm + 3*ss);


%now set to empty all fields of the structures at that frame
clicks_to_remove = idx_clicked_frames(misclicking);

for i = clicks_to_remove
    
    fields = fieldnames(clclk.clicked_frames);
    for j = 1:numel(fields)
        clclk.clicked_frames(i).(fields{j}) = [];
    end
    
    fields = fieldnames(clclk.points);
    for j = 1:numel(fields)
        clclk.points(i).(fields{j}) = [];
    end
    
end %for


idx_clicked_frames(misclicking) = [];





%% plots for debugging
if flag_debugging
    for i = idx_clicked_frames
        
        
        x = clclk.points(i).cilium_sorted_x;
        y = clclk.points(i).cilium_sorted_y;
        
        cla
        hold on
        plot(clclk.points(i).cilium_xx,clclk.points(i).cilium_yy);
        colorset = cool(numel(x));
        
        arrayfun(@(j)plot(x(j),y(j),...
            'o','Color',colorset(j,:),'MarkerFaceColor','auto'),...
            1:numel(x));
        title(i)
        xlabel(['Sum of residuals = ',num2str(Sx.normr+Sy.normr)]);
        pause
    end
    
end


%% code for selecting colormap

if isempty(colormap_name)
    colormap_list = {'blue_gradient';...
        'parula';...
        'jet';...
        'hsv';...
        'hot';...
        'cool';...
        'spring';...
        'summer';...
        'autumn';...
        'winter';...
        'gray';...
        'bone';...
        'copper';...
        'pink';...
        'redblue';...
        };
    
    [colormap_selection, ok_status] = listdlg('PromptString','Select a colormap:',...
        'SelectionMode','single',...
        'ListString',colormap_list);
    
    colormap_chosen = colormap_list{colormap_selection};
else
    colormap_chosen = colormap_name;
end %if


%% GUI for selecting first and last clicked frame

hf = figure;
ha = axes('Parent',hf,...
    'Units','normalized',...
    'Position',[0.05 0.1 0.9 0.85],...
    'Box','on',...
    'NextPlot','add',...
    'DataAspectRatioMode','manual',...
    'PlotboxAspectRatioMode','manual',...
    'XTick',[],'YTick',[],...
    'YDir','Reverse');

% sliders
hsld1 = uicontrol(hf,'Style', 'slider',...
    'Min',firstclickedframe,'Max',lastclickedframe,'Value',firstclickedframe,...
    'Units','normalized',...
    'SliderStep',[1/(lastclickedframe-firstclickedframe), 1/(lastclickedframe-firstclickedframe)],...
    'Position', [0.15 0.05 0.7 0.04],...
    'Callback', @updatefigure);

hsld2 = uicontrol(hf,'Style', 'slider',...
    'Min',firstclickedframe,'Max',lastclickedframe,'Value',lastclickedframe,...
    'SliderStep',[1/(lastclickedframe-firstclickedframe), 1/(lastclickedframe-firstclickedframe)],...
    'Units','normalized',...
    'Position', [0.15 0.01 0.7 0.04],...
    'Callback', @updatefigure);

hl = addlistener([hsld1, hsld2], 'Value', 'PostSet',@updatefigure );

% savebutton
hsb1 = uicontrol(hf,'Style', 'PushButton',...
    'Units','normalized',...
    'Position', [0.85 0.01 0.1 0.09],...
    'String','Exp Fig',...
    'FontName','Helvetica',...
    'FontSize',11,...
    'FontWeight','b',...
    'BackgroundColor',[1 0.2 0.2],...
    'Callback', @savefigure);

hsb2 = uicontrol(hf,'Style', 'PushButton',...
    'Units','normalized',...
    'Position', [0.05 0.01 0.1 0.09],...
    'String','Exp Vid',...
    'FontName','Helvetica',...
    'FontSize',11,...
    'FontWeight','b',...
    'BackgroundColor',[0.2 1 0.2],...
    'Callback', @createvideo);

% arrayfun(@(i) plot(clclk.points(i).cilium_xx,clclk.points(i).cilium_yy,'LineWidth',3),idxclicked)
updatefigure

drawnow
ha.XLim = ha.XLim + [-10 10];
ha.YLim = ha.YLim + [-10 10];
ha.XLimMode = 'manual';
ha.YLimMode = 'manual';


%% update figure GUI

    function [] = updatefigure(~,~)
        %     function [] = updatefigure(source,callbackdata)
        
        minplotframe = floor(hsld1.Value);
        maxplotframe = ceil(hsld2.Value);
        numplotframe = maxplotframe - minplotframe + 1;
        
        cla(gca)
        delete(findobj(ha.Parent.Children,'Type','colorbar'))
        
        
        plotidx = idx_clicked_frames(idx_clicked_frames >= minplotframe & idx_clicked_frames <= maxplotframe);
        
        if strcmp(colormap_chosen,'blue_gradient')
            colorset = [ zeros(numplotframe,2), ones(numplotframe,1), linspace(0.2,1,numplotframe)' ]; %alphablue
        else
            colorset = eval([colormap_chosen,'(',num2str(numplotframe),')']);
        end %if
        
        
%         colorset = [ repmat(linspace(0.9,0,numplotframe)',1,2), ones(numplotframe,1) ];
        
        
                arrayfun(@(i) plot(gca,clclk.points(i).cilium_xx,clclk.points(i).cilium_yy,...
                    'LineWidth',3,'Color',colorset(i - minplotframe + 1,:)),plotidx);
                
                
        % curvature-coded
%         arrayfun(@(i) scatter(gca,clclk.points(i).cilium_xx,clclk.points(i).cilium_yy,...
%             10,clclk.points(i).curvature,'filled'),...
%             plotidx);
        
        
        %         arrayfun(@(i) plot(gca,clclk.points(i).cilium_x,clclk.points(i).cilium_y,...
        %             'ok'),plotidx);
        
        % scalebar
        plot(gca,[scalebar_dx_pos - scalebar_length_px, scalebar_dx_pos], scalebar_y_pos.*[1 1],...
            'LineWidth',5,'Color','k');
        
        % colorbar/colormap
        if ~strcmp(colormap_chosen,'blue_gradient') % doesn't work with alpha values
            colormap(colorset);
            ha.CLim = [0 numplotframe] ./ clclk.movie_object.FrameRate;
            
            hc = colorbar(gca,'south');
            hc.Position(2) = ha.Position(2);
            hc.Position(4) = hc.Position(4)/2;
            hc.Label.String = 'Time, [s]';
            hc.Label.FontSize = 10;
            
            
%             hc = colorbar(gca,'eastoutside');
%             hc.Label.String = 'Time, [s]';
%             hc.Label.FontSize = 12;
            
        end %if
        
    end


%% save figure

    function [] = savefigure(~,~)
        
        minplotframe = floor(hsld1.Value);
        maxplotframe = ceil(hsld2.Value);
        
        hf2 = figure;
        h_dummy = copyobj([ha, findobj(ha.Parent.Children,'Type','colorbar')], hf);
        ha2 = h_dummy(1);
        hc2 = h_dummy(2);
        ha2.Parent = hf2;
        ha2.XTick = [];
        ha2.YTick = [];
        
%         hc2.Parent = hf2;
        
        %         drawnow
        
        
        
        [savepath, savename, ~] = fileparts(clclkfullfilename);
        
        savename = [savename,'_[',num2str(minplotframe),',',num2str(maxplotframe),']'];
        
        saveas(hf2,fullfile(savepath,[savename,'.fig']))
        
        hf2.Color = 'none';
        hf2.PaperPositionMode = 'auto';
        hf2.PaperUnits = 'normalized';
        hf2.InvertHardcopy = 'off';
        hf2.Renderer = 'painter';
        
        print(hf2,fullfile(savepath,[savename,'.svg']),'-dsvg')
        winopen(savepath)
    end


%% create video

    function [] = createvideo(~,~)
        
        % find extremes for saving
        minplotframe = floor(hsld1.Value);
        maxplotframe = ceil(hsld2.Value);
        
        % create savename
        [savepath, savename, ~] = fileparts(clclkfullfilename);
        savename = [savename,'_[',num2str(minplotframe),',',num2str(maxplotframe),']'];
        fullsavename = fullfile(savepath,[savename,'.avi']);
        [savename, savepath] = uiputfile('*.*','Save clclk movie for display purposes',fullsavename);
        fullsavename = fullfile(savepath,savename);
        
        % create video object
        wo = VideoWriter(fullsavename,'Uncompressed AVI');
        
        % set video properties
        wo.FrameRate = 5;
        
        % open it
        wo.open;
        
        % set contrast on clicked frames and convert to uint8
        fs = cat(3, clclk.clicked_frames(idx_clicked_frames).IM); % 3d frame stack
        fs = mat2gray(fs); %double between 0 and 1
        fs = uint8(fs*255);
        
        [height, width, N_frames] = size(fs);
        
        % for loop on clicked frames
        for fc = 1:N_frames
            
            % uint8 frame
            frame = fs(:,:,fc);
            
            % find pixels that belong to cilium to be made red
            
            dummy_px_yy = round(clclk.points(idx_clicked_frames(fc)).cilium_yy);
            dummy_px_xx = round(clclk.points(idx_clicked_frames(fc)).cilium_xx);
            
            % find those out of frame (bc of rounding and interpolation)
            idx_px_oor = dummy_px_yy > height | dummy_px_yy <= 0 |...
                dummy_px_xx > width | dummy_px_xx <= 0;
            
            dummy_px_yy = dummy_px_yy(~idx_px_oor);
            dummy_px_xx = dummy_px_xx(~idx_px_oor);
            
            ciliumind = sub2ind([height, width],...
                dummy_px_yy, dummy_px_xx);
            
            % prepare colours
            red_frame = frame;
            green_frame = frame;
            blue_frame = frame;
            
            % write cilium in the red channel of the rgb frame
            red_frame(ciliumind)   = 255;
            green_frame(ciliumind) = 0;
            blue_frame(ciliumind)  = 0;
            
            % assemble rgb frame
            rgbframe = cat(3, red_frame, green_frame, blue_frame );
            
            % write rgb frame into video
            wo.writeVideo(rgbframe)
            
        end
        
        % close video object
        wo.close;
        
    end

end


function [px2mum, Magnification] = parse_filename_for_px2mum(filename)

if isempty(filename)
    cprintf('*[1 .3 0]','\n Input string was empty.\n');
    px2mum = [];
else
    if contains(filename,'4X','IgnoreCase',true)
        if contains(filename,'1.5X','IgnoreCase',true)
            Magnification = '6X';
            px2mum = 0.9727;
        else
            Magnification = '4X';
            px2mum = 1.459;
        end
    elseif contains(filename,'10X','IgnoreCase',true)
        if contains(filename,'1.5X','IgnoreCase',true)
            Magnification = '15X';
            px2mum = 0.3893;
        else
            Magnification = '10X';
            px2mum = 0.584;
        end
    elseif contains(filename,'20X','IgnoreCase',true)
        if contains(filename,'1.5X','IgnoreCase',true)
            Magnification = '30X';
            px2mum = 0.195;
        else
            Magnification = '20X';
            px2mum = 0.292;
        end
    elseif contains(filename,'30X','IgnoreCase',true)
        if contains(filename,'1.5X','IgnoreCase',true)
            Magnification = '45X';
            px2mum = 0.13;
        else
            Magnification = '30X';
            px2mum = 0.195;
        end
    elseif contains(filename,'40X','IgnoreCase',true)
        if contains(filename,'1.5X','IgnoreCase',true)
            Magnification = '60X';
            px2mum = 0.0973;
        else
            Magnification = '40X';
            px2mum = 0.146;
        end
    elseif contains(filename,'60X','IgnoreCase',true)
        if contains(filename,'1.5X','IgnoreCase',true)
            Magnification = '90X';
            px2mum = 0.0647;
        else
            Magnification = '60X';
            px2mum = 0.097;
        end
    else
        cprintf('*[1 .3 0]','\n Couldn''t find magnification string in namefile.\n');
        px2mum = [];
    end
end % if isempty

if exist('Magnification','var')
    disp(['Total Magnification: ',Magnification]);
end

end



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


if sum(~cellfun(@isempty,matches)) > 1
    
    ok = false;
    while ok == false
        [Selection, ok] = listdlg('ListString',matches,...
            'SelectionMode','single',...
            'Name','Select one file',...
            'PromptString','Select the right match',...
            'ListSize',[6*length(matches{1}), 300]);
    end %while
    
    idx_match = Selection;
    
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































