function varargout = gui_shift(varargin)
% GUI_SHIFT MATLAB code for gui_shift.fig
%      GUI_SHIFT, by itself, creates a new GUI_SHIFT or raises the existing
%      singleton*.
%
%      H = GUI_SHIFT returns the handle to a new GUI_SHIFT or the handle to
%      the existing singleton*.
%
%      GUI_SHIFT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SHIFT.M with the given input arguments.
%
%      GUI_SHIFT('Property','Value',...) creates a new GUI_SHIFT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_shift_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_shift_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_shift

% Last Modified by GUIDE v2.5 11-Jul-2019 18:20:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_shift_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_shift_OutputFcn, ...
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


% --- Executes just before gui_shift is made visible.
function gui_shift_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_shift (see VARARGIN)

% Input values
filename = varargin{1}{1}; % path to fluorescence image
handles.imagedata = imread(filename);
handles.meshes = varargin{1}{2};

% Choose default command line output for gui_shift
handles.output_x = str2double(get(handles.input_xshift, 'String'));
handles.output_y = str2double(get(handles.input_yshift, 'String'));

% zero initial shifts in x and y direction (x: distance to left edge of image, y:
% distance to top of image)
xshift = 0;
yshift = 0;

% plot the shifted data
set(gcf,'toolbar','figure');
imshow(imadjust(handles.imagedata), [],'InitialMagnification','fit');
hold on
for ncell = 1:length(handles.meshes)
    mesh = handles.meshes{ncell}.mesh;
    plot(mesh(:,1)+xshift,mesh(:,2)+yshift,'Color','cyan','MarkerSize',10);
    hold on
    plot(mesh(:,3)+xshift,mesh(:,4)+yshift,'Color','cyan','MarkerSize',10);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_shift wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_shift_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output_x;
varargout{2} = handles.output_y;

% delete the figure
delete(handles.figure1);



function input_xshift_Callback(hObject, ~, handles)
% hObject    handle to input_xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_xshift as text
%        str2double(get(hObject,'String')) returns contents of input_xshift as a double

xshift = str2double(get(hObject,'String'));
yshift = str2double(get(handles.input_yshift,'String'));

% Choose default command line output for gui_shift
handles.output_x = xshift;
handles.output_y = yshift;

% plot the shifted data
imshow(imadjust(handles.imagedata), [],'InitialMagnification','fit');
hold on
for ncell = 1:length(handles.meshes)
    mesh = handles.meshes{ncell}.mesh;
    plot(mesh(:,1)+xshift,mesh(:,2)+yshift,'Color','cyan','MarkerSize',10);
    hold on
    plot(mesh(:,3)+xshift,mesh(:,4)+yshift,'Color','cyan','MarkerSize',10);
end


% --- Executes during object creation, after setting all properties.
function input_xshift_CreateFcn(hObject, ~, ~)
% hObject    handle to input_xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_yshift_Callback(hObject, ~, handles)
% hObject    handle to input_yshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_yshift as text
%        str2double(get(hObject,'String')) returns contents of input_yshift as a double

xshift = str2double(get(handles.input_xshift,'String'));
yshift = str2double(get(hObject,'String'));

% plot the shifted data
imshow(imadjust(handles.imagedata), [],'InitialMagnification','fit');
hold on
for ncell = 1:length(handles.meshes)
    mesh = handles.meshes{ncell}.mesh;
    plot(mesh(:,1)+xshift,mesh(:,2)+yshift,'Color','cyan','MarkerSize',10);
    hold on
    plot(mesh(:,3)+xshift,mesh(:,4)+yshift,'Color','cyan','MarkerSize',10);
end


% --- Executes during object creation, after setting all properties.
function input_yshift_CreateFcn(hObject, ~, ~)
% hObject    handle to input_yshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% stop waiting when user clicks on the button
if isequal(get(handles.figure1, 'waitstatus'),'waiting')
    handles.output_x = str2double(get(handles.input_xshift,'String'));
    handles.output_y = str2double(get(handles.input_yshift,'String'));
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

% Hint: delete(hObject) closes the figure
uiresume(hObject);
