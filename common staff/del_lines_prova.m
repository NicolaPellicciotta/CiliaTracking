f = figure;
hold on;
x=1:100;
y=zeros(10,100);
tracks_to_delete = 0;

for i=1:10,
    y(i,:) = i.*x;
    h(i) = plot(x,y(i,:));
end

% stop button
STOP = 0;
%STOP_ = uicontrol('style','togglebutton','string','STOP','min',0,'max',1,'value',0);
STOP_ = uicontrol('style','pushbutton','string','STOP','Callback','STOP = 1; close(f)');

set(STOP_,'units','normalized','position',[.9 .0 .1 .05]);

UNDO = 0;
UNDO_ = uicontrol('style','togglebutton','string','UNDO','min',0,'max',1,'value',0);
set(UNDO_,'units','normalized','position',[.45 .0 .1 .05]);
while STOP == 0;
    p = gco;
    if ~isempty(p)
        i = find(h == p);
        if tracks_to_delete(end) ~= i
            tracks_to_delete(end+1) = i
        end
        set(h(i),'Visible','off');
    end
    
    UNDO = get(UNDO_,'value');
    
    if UNDO == 1 && length(tracks_to_delete) > 1
        set(h(tracks_to_delete(end)),'Visible','on');
        tracks_to_delete(end) = []
    end
    set(UNDO_,'value',0);
    
    pause(0.05)
    %STOP = get(STOP_,'value');
    
end

tracks_to_delete = tracks_to_delete(2:end); %tolgo lo 0 iniziale
