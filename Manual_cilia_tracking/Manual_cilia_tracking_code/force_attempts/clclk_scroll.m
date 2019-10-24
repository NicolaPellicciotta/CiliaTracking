function [] = clclk_scroll( clclk )
%clclk_scroll Allows to scroll via a slider/mousewheel along the clicked
%frames of a profile video
%   



%% parameters

idx_clicked_frames = [clclk.clicked_frames.frame_number];
N_frames = numel(idx_clicked_frames);

flag_int = isfield('clclk.clicked','dw_IM');




%% GUI

hf = figure;
hf.Color = 'w';

setappdata(hf,'idx_bad_frames',false(size(idx_clicked_frames)));

% callback function for mousewheel

hf.WindowScrollWheelFcn = @scroll_figure;

% position axis
ha = axes('Units','Normalized',...
    'Position',[0 0.1/1.1 1 .9/1.1]);
ha.Tag = 'ha';
hold on;

% show first clicked frame
if flag_int
    hi = imshow(clclk.clicked_frames(idx_clicked_frames(1)).dw_IM,[],'Parent',ha);
else
    hi = imshow(clclk.clicked_frames(idx_clicked_frames(1)).IM,[],'Parent',ha);
end

%show cilium in first frame
if flag_int
    hp = plot(ha,clclk.points(idx_clicked_frames(1)).dw_cilium_xx,...
        clclk.points(idx_clicked_frames(1)).dw_cilium_yy,'b','LineWidth',2,'clipping','off');
else
    hp = plot(ha,clclk.points(idx_clicked_frames(1)).cilium_x,...
        clclk.points(idx_clicked_frames(1)).cilium_y,'g*','MarkerSize',4,'clipping','off');
end

ibd = getappdata(hf,'idx_bad_frames');
if ibd(1) == true
    if flag_int
        hp.LineStyle = '--';
    else
        hp.MarkerStyle = '.';
    end
end

axis image

% ha.XLim = minmax_xx + diff(minmax_xx) * 0.1 * [-1 +1];
% ha.YLim = minmax_yy + diff(minmax_yy) * 0.1 * [-1 +1];

% ha.Visible = 'off';
ha.XTick = [];
ha.YTick = [];

ha.YDir = 'reverse';

ha.Title.String = ['clicked frame = ',num2str(1),', frame = ',num2str(idx_clicked_frames(1))];



% gui elements
hs = uicontrol('Style','Slider',...
    'Min',1,'Max',N_frames,...
    'Value',1,...
    'Units','Normalized',...
    'SliderStep',[1 1],...
    'Position',[0 0 1 0.1/1.1]);

hl = addlistener(hs,'Value','PostSet',@update_figure);


uiwait(hf);


%% callback functions

    function [] = update_figure(~,~)
        
        
        frame2show = round(hs.Value);
        figure(hf);
        axes(ha);
        
        
        % delete all plots
        delete(ha.Children);
        
               
        % show image
        if flag_int
            hi = imshow(clclk.clicked_frames(idx_clicked_frames(frame2show)).dw_IM,[],'Parent',ha);
        else
            hi = imshow(clclk.clicked_frames(idx_clicked_frames(frame2show)).IM,[],'Parent',ha);
        end
        
        % plot cilium
        if flag_int
            hp = plot(ha,clclk.points(idx_clicked_frames(frame2show)).dw_cilium_xx,...
                clclk.points(idx_clicked_frames(frame2show)).dw_cilium_yy,'b','LineWidth',2,'clipping','off');
        else
            hp = plot(ha,clclk.points(idx_clicked_frames(frame2show)).cilium_x,...
                clclk.points(idx_clicked_frames(frame2show)).cilium_y,'g*','MarkerSize',4,'clipping','off');
        end
        
       
        
        
        ha.Title.String = ['clicked frame = ',num2str(frame2show),', frame = ',num2str(idx_clicked_frames(frame2show))];
        
        
        
        
        setappdata(hf,'hp',hp)
        setappdata(hf,'frame2show',frame2show);
        
        
    end

    



end



function [] = scroll_figure(hObject, callbackdata)


frame2show = round(callbackdata.Source.Children(1).Value + callbackdata.VerticalScrollCount);
frame2show = max(1,frame2show);
frame2show = min(frame2show,callbackdata.Source.Children(1).Max);
callbackdata.Source.Children(1).Value = frame2show;

end

%
% function [] = set_button(hObject, callbackdata)
%
%
%
%
% end