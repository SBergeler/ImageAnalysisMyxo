function varargout = gui_cell_tracking(varargin)
% GUI_CELL_TRACKING MATLAB code for gui_cell_tracking.fig
%      GUI_CELL_TRACKING, by itself, creates a new GUI_CELL_TRACKING or raises the existing
%      singleton*.
%
%      H = GUI_CELL_TRACKING returns the handle to a new GUI_CELL_TRACKING or the handle to
%      the existing singleton*.
%
%      GUI_CELL_TRACKING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CELL_TRACKING.M with the given input arguments.
%
%      GUI_CELL_TRACKING('Property','Value',...) creates a new GUI_CELL_TRACKING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_cell_tracking_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_cell_tracking_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_cell_tracking

% Last Modified by GUIDE v2.5 11-Jul-2019 10:16:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_cell_tracking_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_cell_tracking_OutputFcn, ...
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


% --- Executes just before gui_cell_tracking is made visible.
function gui_cell_tracking_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_cell_tracking (see VARARGIN)

filenames = varargin{1}{1}; % paths to phase contrast images
handles.imagedata = cellfun(@imread, filenames, 'uniformoutput', false);
handles.meshes = varargin{1}{2};
handles.cell_ind = varargin{1}{3};
handles.boxes = varargin{1}{4};

% Display the first one 
image_PH = handles.imagedata{1};

this_area = handles.boxes(handles.cell_ind,:); 
min_x = this_area(1);
max_x = this_area(2);
min_y = this_area(3);
max_y = this_area(4);
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

imshow(imadjust(zoomed_fluor));
xlabel(['frame ' num2str(1)]);
if(handles.cell_ind <= length(handles.meshes{1}) && ~isempty(handles.meshes{1}{handles.cell_ind})) 
    hold on
    % plot cell contour
    mesh = handles.meshes{1}{handles.cell_ind}.mesh;
    plot(mesh(:,1)-min_x+1,mesh(:,2)-min_y+1,'Color','cyan','MarkerSize',10);
    hold on
    plot(mesh(:,3)-min_x+1,mesh(:,4)-min_y+1,'Color','cyan','MarkerSize',10);
end

% Update the slider to accomodate all of the images
set(handles.frameslider, 'Min', 1, 'Max', numel(filenames), ...
    'SliderStep', [1 1]/(numel(filenames) - 1), 'Value', 1)

% Choose default command line output for gui_cell_tracking
handles.output = get(handles.frameslider, 'value');
handles.output_delete_daughters = false;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_cell_tracking wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_cell_tracking_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% chose current frame value as command line output
varargout{1} = handles.output; %get(handles.frameslider, 'value');
varargout{2} = handles.output_delete_daughters; 

% delete the figure
delete(handles.figure1);


% --- Executes on slider movement.
function frameslider_Callback(hObject, ~, handles)
% hObject    handle to frameslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% get the value of the current frame
slideval = get(hObject, 'Value');
frame = round(slideval);

image_PH = handles.imagedata{frame};

this_area = handles.boxes(handles.cell_ind,:);
min_x = this_area(1);
max_x = this_area(2);
min_y = this_area(3);
max_y = this_area(4);
zoomed_PH = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

imshow(imadjust(zoomed_PH))
xlabel(['frame ' num2str(frame)]);
if(frame <= length(handles.meshes) && handles.cell_ind <= length(handles.meshes{frame}))
    hold on
    % plot cell contour
    if(~isempty(handles.meshes{frame}{handles.cell_ind}))
        mesh = handles.meshes{frame}{handles.cell_ind}.mesh;
        plot(mesh(:,1)-min_x+1,mesh(:,2)-min_y+1,'Color','cyan','MarkerSize',10);
        hold on
        plot(mesh(:,3)-min_x+1,mesh(:,4)-min_y+1,'Color','cyan','MarkerSize',10);
    end
    % plot cell contour of daughter cells (if existent)
    indd = [];
    for i = 1:length(handles.meshes{frame})
        if ~isempty(handles.meshes{frame}{i}) && ~isempty(handles.meshes{frame}{i}.ancestors)
            if handles.meshes{frame}{i}.ancestors == handles.cell_ind
                indd = [indd, i];
            end
        end
    end
    if ~isempty(indd)
        for i = 1:length(indd)
            mesh = handles.meshes{frame}{indd(i)}.mesh;
            hold on 
            plot(mesh(:,1)-min_x+1,mesh(:,2)-min_y+1,'Color','red','MarkerSize',10);
            hold on
            plot(mesh(:,3)-min_x+1,mesh(:,4)-min_y+1,'Color','red','MarkerSize',10);
        end
    end
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
function figure1_CloseRequestFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% close the figure
% uiresume(hObject);

% stop waiting when user clicks on the button
if isequal(get(handles.figure1, 'waitstatus'),'waiting')
    handles.output = length(handles.meshes);
    uiresume(handles.figure1);
else
    delete(handles.figure1); 
end

guidata(hObject, handles);



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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, ~, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output_delete_daughters = true;

guidata(hObject, handles);
