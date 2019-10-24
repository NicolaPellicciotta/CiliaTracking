function [hf] = scroll_stack( fs )
%scroll_stack Allows to scroll via a slider along the third dimension of a
%3D matrix

if ndims(fs) ~= 3
    error('This function only works with 3D matrices.');
end %if


%% parameters

[hh, ww, N_frames] = size(fs);
scrsz = get(0,'ScreenSize');

scrw = scrsz(3); %screen width
scrh = scrsz(4); %screen height


%% GUI

hf = figure;

% position figure in the middle of the screen (not robust against the
% frames being bigger than the screen)
hf.Position(3) = ww;
hf.Position(4) = hh;
hf.Position(1) = round( (scrw-ww)/2 );
hf.Position(2) = round( (scrh-hh*1.1)/2 );
hf.WindowScrollWheelFcn = @scroll_figure;

% position axis
ha = axes('Units','Normalized',...
    'Position',[0 0.1/1.1 1 1/1.1]);
ha.Tag = 'ha';

%show first frame
hi = imshow(fs(:,:,1),[],'Parent',ha);


% gui elements
hs = uicontrol('Style','Slider',...
    'Min',1,'Max',N_frames,...
    'Value',1,...
    'Units','Normalized',...
    'Position',[0 0 1 0.1/1.1]);

hl = addlistener(hs,'Value','PostSet',@update_figure);

%% callback functions

    function [] = update_figure(~,~)
        
        
        frame2show = round(hs.Value);
        figure(hf);
        dummy_frame = fs(:,:,frame2show);
        
        if isfloat(fs)
            dummy_frame = mat2gray(dummy_frame);
        end
        
        Low_High = stretchlim(dummy_frame(dummy_frame~=0), [0 1]);
        
        hi = imshow(imadjust(dummy_frame,Low_High),[]);
        
    end



end



function [] = scroll_figure(hObject, callbackdata)


frame2show = round(callbackdata.Source.Children(1).Value + callbackdata.VerticalScrollCount);
frame2show = max(1,frame2show);
frame2show = min(frame2show,callbackdata.Source.Children(1).Max);
callbackdata.Source.Children(1).Value = frame2show;

end