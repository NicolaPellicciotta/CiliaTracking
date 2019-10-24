function varargout = cilia_clicky_thingy_GUI_nomovie(varargin)
% CILIA_CLICKY_THINGY_GUI MATLAB code for cilia_clicky_thingy_GUI.fig
%      CILIA_CLICKY_THINGY_GUI, by itself, creates a new CILIA_CLICKY_THINGY_GUI or raises the existing
%      singleton*.
%
%      H = CILIA_CLICKY_THINGY_GUI returns the handle to a new CILIA_CLICKY_THINGY_GUI or the handle to
%      the existing singleton*.
%
%      CILIA_CLICKY_THINGY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CILIA_CLICKY_THINGY_GUI.M with the given input arguments.
%
%      CILIA_CLICKY_THINGY_GUI('Property','Value',...) creates a new CILIA_CLICKY_THINGY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cilia_clicky_thingy_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cilia_clicky_thingy_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cilia_clicky_thingy_GUI

% Last Modified by GUIDE v2.5 09-Feb-2016 19:00:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @cilia_clicky_thingy_GUI_nomovie_OpeningFcn, ...
    'gui_OutputFcn',  @cilia_clicky_thingy_GUI_nomovie_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cilia_clicky_thingy_GUI is made visible.
function cilia_clicky_thingy_GUI_nomovie_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cilia_clicky_thingy_GUI (see VARARGIN)

% save frames in a format .frames
[handles.filename, handles.filepath] = uigetfile('*.mat');

a = load(fullfile(handles.filepath,handles.filename));
handles.fs=a.fs; 

%handles.filename='test';handles.filepath = pwd;

handles.movie_object.FrameRate= 400; %%%% NICOLA  

handles.movie_object.NumberOfFrames= size(handles.fs,3);


% handles.filepath='E:\Data\Cilia_Amelia\2015_12_11\2015_12_11_profile';
% handles.filename='60X_1.5X_BF_profile_CC15_49B.11Dec2015_16.42.34.movie';
%end debug

%Nhandles.movie_object = moviereader(fullfile(handles.filepath,handles.filename));

handles.fc = 1; %frame counter
%N%%handles.total_frames.String = ['/',num2str(handles.movie_object.NumberOfFrames)];
handles.total_frames.String = ['/',num2str(size(handles.fs,3))];


handles.approximate_clickrate = 300; %fps: this will be used to determine automatically the separation
% between clickable frames
handles.edit_clickrate.String = num2str(handles.approximate_clickrate);

%Nhandles.frame_step = max(1,round(handles.movie_object.FrameRate / handles.approximate_clickrate));
handles.frame_step = max(1,round(handles.movie_object.FrameRate / handles.approximate_clickrate));
%max with 1 is if we have to click any movie with less than clickrate fps

% here we'll put all the points
handles.points = repmat(struct('cilium_x',[],'cilium_y',[]), handles.movie_object.NumberOfFrames ,1 );

% here we'll put the frames we have licked (so no need to find original
% video)
handles.clicked_frames = repmat(struct('frame_number',[],'timestamp',[],'IM',[]),...
    handles.movie_object.NumberOfFrames ,1 );

% variable for toggle switch value to show either frame or fluctuation from
% local temporal mean
handles.alternative_view = false;

% interval on which to take time average of frames for alternative view
handles.avint = 7;

% creates savename
[~,name,ext] = fileparts(handles.filename);
handles.savepath = handles.filepath;
handles.savename = [name, '.clclk'];

% Choose default command line output for cilia_clicky_thingy_GUI
handles.output = hObject;

% Update handles structure
update_updateable_nomovie(hObject, handles);

% UIWAIT makes cilia_clicky_thingy_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = cilia_clicky_thingy_GUI_nomovie_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%% Options Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function edit_clickrate_Callback_nomovie(hObject, eventdata, handles)
% hObject    handle to edit_clickrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_clickrate as text
%        str2double(get(hObject,'String')) returns contents of edit_clickrate as a double
temp_clickrate = str2double(get(hObject,'String'));

if temp_clickrate <= handles.movie_object.FrameRate && temp_clickrate > 0
    handles.frame_step = max(1,round(handles.movie_object.FrameRate / temp_clickrate));
    handles.approximate_clickrate = temp_clickrate;
    update_updateable_nomovie(hObject,handles)
elseif temp_clickrate <= 0
    warndlg('Please insert a number bigger than 0');
elseif temp_clickrate > handles.movie_object.FrameRate
    warndlg(['Exceeding Maximum FrameRate (',num2str(handles.movie_object.FrameRate),' fps)']);
end


% --- Executes during object creation, after setting all properties.
function edit_clickrate_CreateFcn_nomovie(hObject, eventdata, handles)
% hObject    handle to edit_clickrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%%%%%%%%%%%%%%%%%%%%% Movie playing commands %%%%%%%%%%%%%%%%%%%%%%



% --- Executes on button press in previous.
function previous_Callback_nomovie(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp_fc = handles.fc - handles.frame_step;
if temp_fc > 0
    handles.fc = temp_fc;
    update_updateable_nomovie(hObject,handles)
end




% --- Executes on button press in next.
function next_Callback_nomovie(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% increase number frame counter
temp_fc = handles.fc + handles.frame_step;
if temp_fc <= handles.movie_object.NumberOfFrames
    handles.fc = temp_fc;
    update_updateable_nomovie(hObject,handles)
else
    warndlg('Exceeding Maximum Number Of Frames');
end



% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn_nomovie(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key
    case 'rightarrow'
        next_Callback_nomovie(hObject, eventdata, handles)
    case 'leftarrow'
        previous_Callback_nomovie(hObject, eventdata, handles)
    otherwise
end


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn_nomovie(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.VerticalScrollCount
    case 1
        next_Callback_nomovie(hObject, eventdata, handles)
    case -1
        previous_Callback_nomovie(hObject, eventdata, handles)
    otherwise
end


function current_frame_Callback_nomovie(hObject, eventdata, handles)
% hObject    handle to current_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_frame as text
%        str2double(get(hObject,'String')) returns contents of current_frame as a double

temp_fc = str2double(get(hObject,'String'));
if temp_fc <= handles.movie_object.NumberOfFrames && temp_fc > 0
    handles.fc = round(temp_fc/handles.frame_step)*handles.frame_step+1;
    update_updateable_nomovie(hObject,handles)
elseif temp_fc <= 0
    warndlg('Please insert a number bigger than 0');
elseif temp_fc > handles.movie_object.NumberOfFrames
    warndlg(['Exceeding Maximum Number Of Frames (',num2str(handles.movie_object.NumberOfFrames),')']);
end



% --- Executes during object creation, after setting all properties.
function current_frame_CreateFcn_nomovie(hObject, eventdata, handles)
% hObject    handle to current_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%%%%%%%%%%%%%%%%% points acquisition functions %%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in delete_frame_points.
function delete_frame_points_Callback_nomovie(hObject, eventdata, handles)
% hObject    handle to delete_frame_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.clicked_frames(handles.fc).frame_number = [];
handles.clicked_frames(handles.fc).timestamp = [];
handles.clicked_frames(handles.fc).IM = [];
handles.points(handles.fc).cilium_x = [];
handles.points(handles.fc).cilium_y = [];
update_updateable_nomovie(hObject,handles)


% --- Executes on button press in click_points.
function click_points_Callback_nomovie(hObject, eventdata, handles)
% hObject    handle to click_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[temp_x,temp_y] = my_getpts(handles.axes1);

if isempty(handles.clicked_frames(handles.fc).IM)
%N    [im,ts] = handles.movie_object.read(handles.fc);
    [im,ts] = handles.fs(:,:,(handles.fc));
    handles.clicked_frames(handles.fc).frame_number = handles.fc;
    handles.clicked_frames(handles.fc).timestamp = ts;
    handles.clicked_frames(handles.fc).IM = im;
end

handles.points(handles.fc).cilium_x = [handles.points(handles.fc).cilium_x; temp_x];
handles.points(handles.fc).cilium_y = [handles.points(handles.fc).cilium_y; temp_y];

update_updateable_nomovie(hObject,handles)



% --- Executes on button press in save_and_exit.
function save_and_exit_Callback_nomovie(hObject, eventdata, handles)
% hObject    handle to save_and_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.savename, handles.savepath] = uiputfile('*.clclk','Save clicked cilia...',...
    fullfile(handles.savepath,handles.savename));

clickyres.filepath                  = handles.filepath;
clickyres.filename                  = handles.filename;
clickyres.movie_object              = handles.movie_object;

clickyres.approximate_clickrate     = handles.approximate_clickrate;
clickyres.frame_step                = handles.frame_step;
clickyres.clicked_frames            = handles.clicked_frames;
clickyres.points                    = handles.points;
clickyres.savepath                  = handles.savepath;
clickyres.savename                  = handles.savename;


save('-v7.3', fullfile(clickyres.savepath,clickyres.savename), '-struct', 'clickyres');
delete(handles.figure1)












%% Updates all that there is to update
function update_updateable_nomovie(hObject, handles)
update_frame_and_points_nomovie(hObject, handles)
handles.current_frame.String = num2str(handles.fc);
handles.static_frame_step.String = num2str(handles.frame_step);
% Update handles structure
guidata(hObject, handles);


%% Updates the image and points
function update_frame_and_points_nomovie(hObject, handles)

if handles.alternative_view == false
    imshow(imadjust(handles.fs(:,:,(handles.fc))),[]);
else %if toggle for alternative view is switched on
%N    cf = double(handles.movie_object.read(handles.fc)); % current frame
    cf = double(handles.fs(:,:,handles.fc)); % current frame
    ds = max(0, handles.fc - floor(handles.avint/2));
    de = min(handles.movie_object.NumberOfFrames, handles.fc + floor(handles.avint/2));
    dfs = double(handles.fs(:,:,ds:de)); %read short stack
    mdf = mean(dfs,3); %time average it
    imshow(cf - mdf,[]); %show difference between current frame and local time-average
end

hold(handles.axes1,'on');

frames_to_plot = handles.fc - [2, 1, 0] .* handles.frame_step;
colours = cool(2*numel(frames_to_plot));
MarkerSize = 2;


for i = 1:numel(frames_to_plot)
    if frames_to_plot(i) > 0
        plot(handles.points(frames_to_plot(i)).cilium_x, handles.points(frames_to_plot(i)).cilium_y,...
            'Marker','x',...
            'MarkerSize',MarkerSize,...
            'LineStyle','none',...
            'Color',colours(numel(frames_to_plot) + i,:))
    end
end

hold(handles.axes1,'off');


% --- Executes on button press in alternative_view.
function alternative_view_Callback_nomovie(hObject, eventdata, handles)
% hObject    handle to alternative_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alternative_view
handles.alternative_view = get(hObject,'Value');
update_updateable_nomovie(hObject,handles)



