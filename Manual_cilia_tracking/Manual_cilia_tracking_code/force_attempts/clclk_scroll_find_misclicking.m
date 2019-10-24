function [idx_bad_frames] = clclk_scroll_find_misclicking( clclk )
%clclk_scroll_find_misclicking Allows to scroll via a slider/mousewheel along the clicked
%frames of a profile video to find wrong clicks
%   Use this *before* calculating the displacements, to check on which
%   "frames" you want to ignore




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


hb = uicontrol('Style','radiobutton',...
    'Units','normalized',...
    'Position',[0 1/1.1 0.1 0.05],...
    'BackgroundColor','w',...
    'String','bad click');
hbl = addlistener(hb,'Value','PostSet',@update_bad_clicks);


hs = uicontrol('Style','Slider',...
    'Min',1,'Max',N_frames,...
    'Value',1,...
    'Units','Normalized',...
    'SliderStep',[1 1],...
    'Position',[0 0 1 0.1/1.1]);

hl = addlistener(hs,'Value','PostSet',@update_figure);

% idx_bad_frames = getappdata(hf,'idx_bad_frames');

uiwait(hf);

if nargout < 1
    clear idx_bad_frames
end

%% callback functions

    function [] = update_figure(~,~)
        
        idbafr = getappdata(hf,'idx_bad_frames');
        
        frame2show = round(hs.Value);
        figure(hf);
        axes(ha);
        
        
        % delete all plots
        delete(ha.Children);
        
        % get slider value
        hb.Value = idbafr(frame2show);
        
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
                clclk.points(idx_clicked_frames(frame2show)).cilium_y,'g*','MarkerSize',6,'clipping','off');
        end
        
        % and change style if bad
        if idbafr(frame2show) == true
            if flag_int
                hp.LineStyle = '--';
            else
                hp.Marker = '.';
            end
        end
        
        
        ha.Title.String = ['clicked frame = ',num2str(frame2show),', frame = ',num2str(idx_clicked_frames(frame2show))];
        
        
        if frame2show > 1
            
            %             hold on
            if flag_int
                hp(2) = plot(ha,clclk.points(idx_clicked_frames(frame2show-1)).dw_cilium_xx,...
                    clclk.points(idx_clicked_frames(frame2show-1)).dw_cilium_yy,'Color',[.7 .7 .7],'LineWidth',2,'clipping','off');
            else
                hp(2) = plot(ha,clclk.points(idx_clicked_frames(frame2show-1)).cilium_x,...
                    clclk.points(idx_clicked_frames(frame2show-1)).cilium_y,'*','Color',[.5 1 .6],'MarkerSize',2,'clipping','off');
            end
            
            if idbafr(frame2show-1) == true
                if flag_int
                    hp(2).LineStyle = '--';
                else
                    hp(2).Marker = '.';
                end
            end
            
        end
        
        if frame2show > 2
            if flag_int
            hp(3) = plot(ha,clclk.points(idx_clicked_frames(frame2show-2)).dw_cilium_xx,...
                clclk.points(idx_clicked_frames(frame2show-2)).dw_cilium_yy,'Color',[.85 .85 .85],'LineWidth',2,'clipping','off');
            else
                hp(3) = plot(ha,clclk.points(idx_clicked_frames(frame2show-2)).cilium_x,...
                    clclk.points(idx_clicked_frames(frame2show-2)).cilium_y,'*','Color',[.7 1 .8],'MarkerSize',2,'clipping','off');
           
            end
            
            if idbafr(frame2show-2) == true
                if flag_int
                    hp(3).LineStyle = '--';
                else
                    hp(3).Marker = '.';
                end
            end
        end
        
        setappdata(hf,'hp',hp)
        setappdata(hf,'frame2show',frame2show);
        setappdata(hf,'idx_bad_frames',idbafr);
        
        idx_bad_frames = getappdata(hf,'idx_bad_frames');
        
    end

    function idx_bad_frames = update_bad_clicks(~,~)
        
        frame2show = getappdata(hf,'frame2show');
        
        idx_bad_frames = getappdata(hf,'idx_bad_frames');
        idx_bad_frames(frame2show) = hb.Value;
        
        setappdata(hf,'idx_bad_frames',idx_bad_frames);
        
        set(findobj(gca,'Type','Line','Color','b'),'LineStyle','--');
        
        
        
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