function [ which_stroke, ind_clicked ] = clclk_select_whichstroke( clclk )
%clclk_select_powerstroke GUI for choosing which frames are power/recovery
%stroke. Returns 1 for power, 2 for recovery
%   Needs clclk after it has been read by clclk_reader



%% remove misclicking (because this is only done as a display thing in clclk_reader)

ind_clicked = [clclk.clicked_frames.frame_number];

% mm = median([clclk.points(ind_clicked).normr]);
% ss = iqr([clclk.points(ind_clicked).normr]);
% ind_clicked([clclk.points(ind_clicked).normr] < mm - 3*ss | [clclk.points(ind_clicked).normr] > mm + 3*ss) = [];


normr_array = [clclk.points(ind_clicked).normr];
fluct_normr_array = normr_array - reshape(smooth(normr_array,5),1,[]);
mm = median(fluct_normr_array);
ss = iqr(fluct_normr_array);
ind_clicked(fluct_normr_array > mm + 3*ss) = [];

N_frames = numel(ind_clicked);

%% GUI

hf = figure;
hf.Color = 'w';

setappdata(hf,'idx_power_stroke',false(size(ind_clicked)));


% callback function for mousewheel

hf.WindowScrollWheelFcn = @scroll_figure;

% position axis
ha = axes('Units','Normalized',...
    'Position',[0 0.1/1.1 1 .9/1.1]);
ha.Tag = 'ha';
hold on;


%show cilium in first frame
hp = plot(ha,clclk.points(ind_clicked(1)).cilium_xx,...
    clclk.points(ind_clicked(1)).cilium_yy,'g','LineWidth',2,'clipping','off');


ips = getappdata(hf,'idx_power_stroke');
if ips(1) == true    
    hp.Color = 'r';
else
    hp.Color = 'b';
end


axis image

minmax_xx = minmax(vertcat(clclk.points(ind_clicked).cilium_xx));
minmax_yy = minmax(vertcat(clclk.points(ind_clicked).cilium_yy));

ha.XLim = minmax_xx + diff(minmax_xx) * 0.1 * [-1 +1];
ha.YLim = minmax_yy + diff(minmax_yy) * 0.1 * [-1 +1];

% ha.Visible = 'off';
ha.XTick = [];
ha.YTick = [];

ha.YDir = 'reverse';

ha.Title.String = ['clicked frame = ',num2str(1),', frame = ',num2str(ind_clicked(1))];



% gui elements


hb = uicontrol('Style','radiobutton',...
    'Units','normalized',...
    'Position',[0 1/1.1 0.1 0.05],...
    'BackgroundColor','w',...
    'String','power stroke');
hbl = addlistener(hb,'Value','PostSet',@update_power_stroke);


hs = uicontrol('Style','Slider',...
    'Min',1,'Max',N_frames,...
    'Value',1,...
    'Units','Normalized',...
    'SliderStep',[1 1],...
    'Position',[0 0 1 0.1/1.1]);

hl = addlistener(hs,'Value','PostSet',@update_figure);

% idx_power_stroke = getappdata(hf,'idx_power_stroke');

uiwait(hf);

which_stroke = nan(size(ind_clicked));
which_stroke(idx_power_stroke) = 1;
which_stroke(~idx_power_stroke) = 2;

if nargout < 1
    clear idx_power_stroke
end

%% callback functions

    function [] = update_figure(~,~)
        
        idpowstr = getappdata(hf,'idx_power_stroke');
        
        frame2show = round(hs.Value);
        figure(hf);
        axes(ha);
        
        
        % delete all plots
        delete(ha.Children);
        
        % get slider value
        hb.Value = idpowstr(frame2show);
        
                
        % plot cilium
        
        hp = plot(ha,clclk.points(ind_clicked(frame2show)).cilium_xx,...
            clclk.points(ind_clicked(frame2show)).cilium_yy,'g','LineWidth',2,'clipping','off');
        
        % and change style if bad
        if idpowstr(frame2show) == true
            hp.Color = 'r';
        else
            hp.Color = 'b';
        end
        
        
        ha.Title.String = ['clicked frame = ',num2str(frame2show),', frame = ',num2str(ind_clicked(frame2show))];
        
        
        if frame2show > 1
            
            %             hold on
            hp(2) = plot(ha,clclk.points(ind_clicked(frame2show-1)).cilium_xx,...
                clclk.points(ind_clicked(frame2show-1)).cilium_yy,'Color','g',...
                'LineWidth',2,'LineStyle','--','clipping','off');
            
            if idpowstr(frame2show-1) == true
                hp(2).Color = 'r';
            else
                hp(2).Color = 'b';
            end
            
        end
        
        if frame2show > 2
            
            hp(3) = plot(ha,clclk.points(ind_clicked(frame2show-2)).cilium_xx,...
                clclk.points(ind_clicked(frame2show-2)).cilium_yy,'Color','g',...
                'LineWidth',2,'LineStyle',':','clipping','off');
            
            if idpowstr(frame2show-2) == true
                hp(3).Color = 'r';
            else
                hp(3).Color = 'b';
            end
            
        end
        
        setappdata(hf,'hp',hp)
        setappdata(hf,'frame2show',frame2show);
        setappdata(hf,'idx_power_stroke',idpowstr);
        
        idx_power_stroke = getappdata(hf,'idx_power_stroke');
        
    end

    function idx_power_stroke = update_power_stroke(~,~)
        
        frame2show = getappdata(hf,'frame2show');
        
        idx_power_stroke = getappdata(hf,'idx_power_stroke');
        idx_power_stroke(frame2show:end) = hb.Value;
        
        setappdata(hf,'idx_power_stroke',idx_power_stroke);
        
        if idx_power_stroke(frame2show)
            set(findobj(gca,'Type','Line','LineStyle','-'),'Color','r');
        else
            set(findobj(gca,'Type','Line','LineStyle','-'),'Color','b');
            
        end
        
        
    end



end



function [] = scroll_figure(hObject, callbackdata)


frame2show = round(callbackdata.Source.Children(1).Value + callbackdata.VerticalScrollCount);
frame2show = max(1,frame2show);
frame2show = min(frame2show,callbackdata.Source.Children(1).Max);
callbackdata.Source.Children(1).Value = frame2show;

end

