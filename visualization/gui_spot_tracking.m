function varargout = gui_spot_tracking(varargin)
% GUI_SPOT_TRACKING MATLAB code for gui_spot_tracking.fig
%      GUI_SPOT_TRACKING, by itself, creates a new GUI_SPOT_TRACKING or raises the existing
%      singleton*.
%
%      H = GUI_SPOT_TRACKING returns the handle to a new GUI_SPOT_TRACKING or the handle to
%      the existing singleton*.
%
%      GUI_SPOT_TRACKING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SPOT_TRACKING.M with the given input arguments.
%
%      GUI_SPOT_TRACKING('Property','Value',...) creates a new GUI_SPOT_TRACKING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_spot_tracking_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_spot_tracking_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_spot_tracking

% Last Modified by GUIDE v2.5 11-Jul-2019 10:14:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_spot_tracking_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_spot_tracking_OutputFcn, ...
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


% --- Executes just before gui_spot_tracking is made visible.
function gui_spot_tracking_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_spot_tracking (see VARARGIN)

filenames = varargin{1}{1};
handles.imagedata = cellfun(@imread, filenames, 'uniformoutput', false);
handles.data_spots = varargin{1}{2};
handles.spot_ind = varargin{1}{3};
handles.boxes = varargin{1}{4};

% Display the first one 
image_fluor = handles.imagedata{1};

t = linspace(0,2*pi,50);

this_area = handles.boxes(handles.spot_ind,:); %handles.data_spots{1}(handles.spot_ind).box;
min_x = this_area(1);
max_x = this_area(2);
min_y = this_area(3);
max_y = this_area(4);
%zoomed_fluor = image_fluor(max(1,topy):min(size(image_fluor,1),topy+height),max(1,topx):min(size(image_fluor,2),topx+width));
zoomed_fluor = image_fluor(max(1,min_y):min(size(image_fluor,1),max_y),max(1,min_x):min(size(image_fluor,2),max_x));

%imshow(imadjust(zoomed_fluor), 'Border','tight');
% if(~isempty(handles.data_spots{1}(handles.spot_ind).threshold_spot))
%     max_im = handles.data_spots{1}(handles.spot_ind).threshold_spot;
% else
%     max_im = max(zoomed_fluor(:));
% end
imshow(imadjust(zoomed_fluor));%[min(zoomed_fluor(:)),max_im]);
xlabel(['frame ' num2str(1)]);
if(~isempty(handles.data_spots{1}(handles.spot_ind).x))
    hold on
    a = handles.data_spots{1}(handles.spot_ind).spot_major_axis/2;
    b = handles.data_spots{1}(handles.spot_ind).spot_minor_axis/2;
    Xc = handles.data_spots{1}(handles.spot_ind).x - max(1,min_x) + 1;
    Yc = handles.data_spots{1}(handles.spot_ind).y - max(1,min_y) + 1;
    phi = deg2rad(-handles.data_spots{1}(handles.spot_ind).spot_orient);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    % plot spot position
    plot(Xc, Yc, 'xr', 'MarkerSize',10, 'LineWidth', 1.5)
    hold on
    % plot ellipses fitted to spots
    plot(x,y,'magenta','Linewidth',1.5)
end

% Update the slider to accomodate all of the images
set(handles.frameslider, 'Min', 1, 'Max', numel(filenames), ...
    'SliderStep', [1 1]/(numel(filenames) - 1), 'Value', 1)

% Choose default command line output for gui_spot_tracking
handles.output = get(handles.frameslider, 'value');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_spot_tracking wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_spot_tracking_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% chose current frame value as command line output
varargout{1} = handles.output; %get(handles.frameslider, 'value');

% delete the figure
delete(handles.figure1);


% --- Executes on slider movement.
function frameslider_Callback(hObject, ~, handles)
% hObject    handle to frameslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

t = linspace(0,2*pi,50);

% get the value of the current frame
slideval = get(hObject, 'Value');
frame = round(slideval);

image_fluor = handles.imagedata{frame};

this_area = handles.boxes(handles.spot_ind,:); %handles.data_spots{1}(handles.spot_ind).box; % the box is the same for all frames!
min_x = this_area(1);
max_x = this_area(2);
min_y = this_area(3);
max_y = this_area(4);
zoomed_fluor = image_fluor(max(1,min_y):min(size(image_fluor,1),max_y),max(1,min_x):min(size(image_fluor,2),max_x));
% if(~isempty(handles.data_spots{frame}(handles.spot_ind).threshold_spot))
%     max_im = handles.data_spots{frame}(handles.spot_ind).threshold_spot;
% else
%     max_im = max(zoomed_fluor(:));
% end
imshow(imadjust(zoomed_fluor)) %[min(zoomed_fluor(:)),max_im]);
xlabel(['frame ' num2str(frame)]);
index = arrayfun(@(x) x.spot_ind == handles.spot_ind, handles.data_spots{frame});
if(~all(index == 0)) %~isempty(handles.data_spots{frame}(index).x)
    hold on
    a = handles.data_spots{frame}(index).spot_major_axis/2;
    b = handles.data_spots{frame}(index).spot_minor_axis/2;
    Xc = handles.data_spots{frame}(index).x - max(1,min_x) + 1;
    Yc = handles.data_spots{frame}(index).y - max(1,min_y) + 1;
    phi = deg2rad(-handles.data_spots{frame}(index).spot_orient);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    % plot spot position
    plot(Xc, Yc, 'xr', 'MarkerSize',10, 'LineWidth', 1.5)
    hold on
    % plot ellipses fitted to spots
    plot(x,y,'magenta','Linewidth',1.5)
end


% --- Executes during object creation, after setting all properties.
function frameslider_CreateFcn(hObject, ~, ~)
% hObject    handle to frameslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% stop waiting when user clicks on the button
if isequal(get(handles.figure1, 'waitstatus'),'waiting')
    handles.output = get(handles.frameslider, 'value');
    uiresume(handles.figure1);
else
    delete(handles.figure1);
end

guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% close the figure
uiresume(hObject);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% stop waiting when user clicks on the button
if isequal(get(handles.figure1, 'waitstatus'),'waiting')
    handles.output = 0;
    uiresume(handles.figure1);
else
    delete(handles.figure1); 
end

guidata(hObject, handles);
