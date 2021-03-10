function varargout = gui_spot_detection_timelapse(varargin)
% GUI_SPOT_DETECTION_TIMELAPSE MATLAB code for gui_spot_detection_timelapse.fig
%      GUI_SPOT_DETECTION_TIMELAPSE, by itself, creates a new GUI_SPOT_DETECTION_TIMELAPSE or raises the existing
%      singleton*.
%
%      H = GUI_SPOT_DETECTION_TIMELAPSE returns the handle to a new GUI_SPOT_DETECTION_TIMELAPSE or the handle to
%      the existing singleton*.
%
%      GUI_SPOT_DETECTION_TIMELAPSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SPOT_DETECTION_TIMELAPSE.M with the given input arguments.
%
%      GUI_SPOT_DETECTION_TIMELAPSE('Property','Value',...) creates a new GUI_SPOT_DETECTION_TIMELAPSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_spot_detection_timelapse_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_spot_detection_timelapse_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_spot_detection_timelapse

% Last Modified by GUIDE v2.5 03-Dec-2019 16:15:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_spot_detection_timelapse_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_spot_detection_timelapse_OutputFcn, ...
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


% --- Executes just before gui_spot_detection_timelapse is made visible.
function gui_spot_detection_timelapse_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_spot_detection_timelapse (see VARARGIN)

filenames = varargin{1}{1}; % paths to fluorescence images
handles.imagedata = cellfun(@imread, filenames, 'uniformoutput', false);
handles.matpath = varargin{1}{2};
handles.xshift = varargin{1}{3};
handles.yshift = varargin{1}{4};
handles.cell_id = varargin{1}{5};
handles.box = varargin{1}{6};

handles.masks = load(handles.matpath);

handles.frame = 1;

set(0, 'CurrentFigure', handles.figure1)
set(gcf,'toolbar','figure'); % add toolbar to figure such that one can zoom in and out

image_PH = handles.imagedata{1};
 
min_x = handles.box(1);
max_x = handles.box(2);
min_y = handles.box(3);
max_y = handles.box(4);
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

imshow(imadjust(zoomed_fluor));
xlabel(['cellid ' num2str(handles.cell_id) ', frame ' num2str(1)]);
hold on

% write the lower and upper intensity threshold used by imadjust
% (approximate values) as static texts
low_int_adj = double(min(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
high_int_adj = double(max(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
set(handles.high_int_threshold, 'String', num2str(high_int_adj));
set(handles.low_int_threshold, 'String', num2str(low_int_adj));

ind = (handles.masks.cellList.cellId{1} == handles.cell_id);
set(handles.output_text,'String','Orientation of spot(s) relative to cell:')
if any(ind)
    this_cell = handles.masks.cellList.meshData{1}{ind};
    mesh = this_cell.mesh;
    % change the meshes according to the shifts
    mesh(:,[1 3]) = mesh(:,[1 3]) + handles.xshift;
    mesh(:,[2 4]) = mesh(:,[2 4]) + handles.yshift;
    centerline = [mean([mesh(:,1) mesh(:,3)],2) mean([mesh(:,2) mesh(:,4)],2)];
    % plot the cell surrounding
    plot(mesh(:,1)-min_x+1,mesh(:,2)-min_y+1,'Color','cyan','LineWidth',1.5)
    hold on
    plot(mesh(:,3)-min_x+1,mesh(:,4)-min_y+1,'Color','cyan','LineWidth',1.5)
    hold on
    if isfield(this_cell,'spots_matlab')
        this_cell_spots = this_cell.spots_matlab;
        % plot ellipses fitted to spots
        t = linspace(0,2*pi,50);
        rel_orient = [];
        if ~isempty(this_cell_spots(1).l)
            for n = 1:length(this_cell_spots) % for all spots
                a = this_cell_spots(n).spot_major_axis/2;
                b = this_cell_spots(n).spot_minor_axis/2;
                Xc = this_cell_spots(n).x-min_x+1;
                Yc = this_cell_spots(n).y-min_y+1;
                phi = deg2rad(-this_cell_spots(n).spot_orient);
                x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
                y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
                % plot spot position
                plot(Xc, Yc, 'xr', 'MarkerSize',10, 'LineWidth', 1.5)
                hold on
                % plot spot outline
                %cellfun(@(x) plot(x(:,2),x(:,1),'Color','green','LineWidth',1.5),boundaries)
                %hold on
                % plot ellipses fitted to spots
                plot(x,y,'magenta','Linewidth',1.5)
                hold on
                % plot major axis of ellipse
                % the extremal points of the ellipse are the points at t =
                % 0 and t = pi
                x1 = Xc + a*cos(phi);
                y1 = Yc + a*sin(phi);
                x2 = Xc - a*cos(phi);
                y2 = Yc - a*sin(phi);
                line([x1,x2],[y1,y2],'Color','magenta','LineWidth',1.5);
                hold on
                spot_position = this_cell_spots(n).positions;
                centerline_selected = centerline(max((spot_position-2),1):min((spot_position+2),length(centerline)),:);
                % plot centerline to show the orientation of the cell compared
                % to the spot
                plot(centerline_selected(:,1),centerline_selected(:,2),'xw')
                rel_orient = [rel_orient this_cell_spots(n).rel_orient];
            end
        end
    end
    set(handles.output_text,'String',['Orientation of spot(s) relative to cell:', newline, ...
        num2str(rel_orient)])
end

% Update the slider to accomodate all of the images
ntframes = length(handles.masks.cellList.meshData);
set(handles.frameslider, 'Min', 1, 'Max', ntframes, ...
    'SliderStep', [1 1]/(ntframes - 1), 'Value', 1)

% Choose default command line output for gui_spot_detection_timelapse
handles.output = false; %get(handles.frameslider, 'value');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_spot_detection_timelapse wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_spot_detection_timelapse_OutputFcn(~, ~, handles) 
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

% get the value of the current frame
frame = round(get(hObject, 'Value'));

% add frame to handle to be able to use the value also in the
% other functions
handles.frame = frame;

image_PH = handles.imagedata{frame};
 
min_x = handles.box(1);
max_x = handles.box(2);
min_y = handles.box(3);
max_y = handles.box(4);
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

cla(handles.figure1) % clear the previous graph
imshow(imadjust(zoomed_fluor));
xlabel(['cellid ' num2str(handles.cell_id) ', frame ' num2str(frame)]);
hold on

% write the lower and upper intensity threshold used by imadjust
% (approximate values) as static texts
low_int_adj = double(min(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
high_int_adj = double(max(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
set(handles.high_int_threshold, 'String', num2str(high_int_adj));
set(handles.low_int_threshold, 'String', num2str(low_int_adj));

set(handles.output_text,'String','Orientation of spot(s) relative to cell:')

% index of the cell with cell id 'cell_id' in the current frame
ind = (handles.masks.cellList.cellId{frame} == handles.cell_id);
if any(ind) % when the cell is detected in the current frame
    this_cell = handles.masks.cellList.meshData{frame}{ind};
    mesh = this_cell.mesh;
    % change the meshes according to the shifts
    mesh(:,[1 3]) = mesh(:,[1 3]) + handles.xshift;
    mesh(:,[2 4]) = mesh(:,[2 4]) + handles.yshift;
    centerline = [mean([mesh(:,1) mesh(:,3)],2) mean([mesh(:,2) mesh(:,4)],2)];
    % plot the cell surrounding
    plot(mesh(:,1)-min_x+1,mesh(:,2)-min_y+1,'Color','cyan','LineWidth',1.5)
    hold on
    plot(mesh(:,3)-min_x+1,mesh(:,4)-min_y+1,'Color','cyan','LineWidth',1.5)
    hold on
    if isfield(this_cell,'spots_matlab')
        this_cell_spots = this_cell.spots_matlab;
        % plot ellipses fitted to spots
        t = linspace(0,2*pi,50);
        rel_orient = [];
        if ~isempty(this_cell_spots(1).l)
            for n = 1:length(this_cell_spots) % for all spots
                a = this_cell_spots(n).spot_major_axis/2;
                b = this_cell_spots(n).spot_minor_axis/2;
                Xc = this_cell_spots(n).x-min_x+1;
                Yc = this_cell_spots(n).y-min_y+1;
                phi = deg2rad(-this_cell_spots(n).spot_orient);
                x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
                y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
                % plot spot position
                plot(Xc, Yc, 'xr', 'MarkerSize',10, 'LineWidth', 1.5)
                hold on
                % plot spot outline
                %cellfun(@(x) plot(x(:,2),x(:,1),'Color','green','LineWidth',1.5),boundaries)
                %hold on
                % plot ellipses fitted to spots
                plot(x,y,'magenta','Linewidth',1.5)
                hold on
                % plot major axis of ellipse
                % the extremal points of the ellipse are the points at t =
                % 0 and t = pi
                x1 = Xc + a*cos(phi);
                y1 = Yc + a*sin(phi);
                x2 = Xc - a*cos(phi);
                y2 = Yc - a*sin(phi);
                line([x1,x2],[y1,y2],'Color','magenta','LineWidth',1.5);
                hold on
                spot_position = this_cell_spots(n).positions;
                centerline_selected = centerline(max((spot_position-2),1):min((spot_position+2),length(centerline)),:);
                % plot centerline to show the orientation of the cell compared
                % to the spot
                plot(centerline_selected(:,1),centerline_selected(:,2),'xw')
                rel_orient = [rel_orient this_cell_spots(n).rel_orient];
            end
        end
    end
    set(handles.output_text,'String',['Orientation of spot(s) relative to cell:', newline, ...
        num2str(rel_orient)])
end

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function frameslider_CreateFcn(hObject, ~, ~)
% hObject    handle to frameslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% close the figure
uiresume(hObject);



function low_int_threshold_Callback(hObject, ~, handles)
% hObject    handle to low_int_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of low_int_threshold as text
%        str2double(get(hObject,'String')) returns contents of low_int_threshold as a double
low_int_threshold = str2double(get(hObject,'String'))/100;
high_int_threshold = str2double(get(handles.high_int_threshold,'String'))/100;

% plot the data with the updated intensity thresholds
image_PH = handles.imagedata{handles.frame};
 
min_x = handles.box(1);
max_x = handles.box(2);
min_y = handles.box(3);
max_y = handles.box(4);
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

% find the previous plot of the image and replace it with the new one
im_old = findobj(handles.axes1.Children,'type','image');
delete(im_old);
imshow(imadjust(zoomed_fluor,[low_int_threshold high_int_threshold],[]));
h = get(gca,'Children');
set(gca,'Children',[h(2:end); h(1)])

% --- Executes during object creation, after setting all properties.
function low_int_threshold_CreateFcn(hObject, ~, ~)
% hObject    handle to low_int_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function high_int_threshold_Callback(hObject, ~, handles)
% hObject    handle to high_int_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of high_int_threshold as text
%        str2double(get(hObject,'String')) returns contents of high_int_threshold as a double
high_int_threshold = str2double(get(hObject,'String'))/100;
low_int_threshold = str2double(get(handles.low_int_threshold,'String'))/100;

% plot the data with the updated intensity thresholds
image_PH = handles.imagedata{handles.frame};
 
min_x = handles.box(1);
max_x = handles.box(2);
min_y = handles.box(3);
max_y = handles.box(4);
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

% find the previous plot of the image and replace it with the new one
im_old = findobj(handles.axes1.Children,'type','image');
delete(im_old);
imshow(imadjust(zoomed_fluor,[low_int_threshold high_int_threshold],[]));
h = get(gca,'Children');
set(gca,'Children',[h(2:end); h(1)])


% --- Executes during object creation, after setting all properties.
function high_int_threshold_CreateFcn(hObject, ~, ~)
% hObject    handle to high_int_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reset_button.
function reset_button_Callback(~, ~, handles)
% hObject    handle to reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% plot the data with the original intensity thresholds (from imadjust)
image_PH = handles.imagedata{handles.frame};
 
min_x = handles.box(1);
max_x = handles.box(2);
min_y = handles.box(3);
max_y = handles.box(4);
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

% find the previous plot of the image and replace it with the new one
im_old = findobj(handles.axes1.Children,'type','image');
delete(im_old);
imshow(imadjust(zoomed_fluor));
h = get(gca,'Children');
set(gca,'Children',[h(2:end); h(1)])

low_int_adj = double(min(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
high_int_adj = double(max(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
set(handles.high_int_threshold, 'String', num2str(high_int_adj));
set(handles.low_int_threshold, 'String', num2str(low_int_adj));


% --- Executes on button press in pushbuttonstop.
function pushbuttonstop_Callback(hObject, ~, handles)
% hObject    handle to pushbuttonstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = true;
guidata(hObject, handles);
% stop waiting when user clicks on the button
if isequal(get(handles.figure1, 'waitstatus'),'waiting')
    uiresume(handles.figure1);
else
    delete(handles.figure1); 
end


% --- Executes on button press in pushbuttonstopthiscell.
function pushbuttonstopthiscell_Callback(~, ~, handles)
% hObject    handle to pushbuttonstopthiscell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% stop waiting when user clicks on the button
if isequal(get(handles.figure1, 'waitstatus'),'waiting')
    uiresume(handles.figure1);
else
    delete(handles.figure1); 
end
