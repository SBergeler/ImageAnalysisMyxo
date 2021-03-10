function varargout = gui_spot_detection_indep_frames(varargin)
% GUI_SPOT_DETECTION_INDEP_FRAMES MATLAB code for gui_spot_detection_indep_frames.fig
%      GUI_SPOT_DETECTION_INDEP_FRAMES, by itself, creates a new GUI_SPOT_DETECTION_INDEP_FRAMES or raises the existing
%      singleton*.
%
%      H = GUI_SPOT_DETECTION_INDEP_FRAMES returns the handle to a new GUI_SPOT_DETECTION_INDEP_FRAMES or the handle to
%      the existing singleton*.
%
%      GUI_SPOT_DETECTION_INDEP_FRAMES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SPOT_DETECTION_INDEP_FRAMES.M with the given input arguments.
%
%      GUI_SPOT_DETECTION_INDEP_FRAMES('Property','Value',...) creates a new GUI_SPOT_DETECTION_INDEP_FRAMES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_spot_detection_indep_frames_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_spot_detection_indep_frames_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_spot_detection_indep_frames

% Last Modified by GUIDE v2.5 06-Sep-2019 20:37:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_spot_detection_indep_frames_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_spot_detection_indep_frames_OutputFcn, ...
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


% --- Executes just before gui_spot_detection_indep_frames is made visible.
function gui_spot_detection_indep_frames_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_spot_detection_indep_frames (see VARARGIN)

filenames = varargin{1}{1}; % paths to fluorescence images
handles.imagedata = cellfun(@imread, filenames, 'uniformoutput', false);
handles.matpath = varargin{1}{2};
handles.xshift = varargin{1}{3};
handles.yshift = varargin{1}{4};

handles.masks = load(handles.matpath);

handles.frame = 1;
handles.cell_no = 1;

set(0, 'CurrentFigure', handles.figure1)
set(gcf,'toolbar','figure'); % add toolbar to figure such that one can zoom in and out

this_cell = handles.masks.cellList.meshData{1}{1};
this_area = this_cell.box; % [x, y, width, height]

image_PH = handles.imagedata{1};
 
min_x = this_area(1) + handles.xshift;
max_x = this_area(1) + this_area(3) + handles.xshift;
min_y = this_area(2) + handles.yshift;
max_y = this_area(2) + this_area(4) + handles.yshift;
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

imshow(imadjust(zoomed_fluor));
xlabel(['frame ' num2str(1) ', cellid ' num2str(handles.masks.cellList.cellId{1}(1))]);
hold on

% write the lower and upper intensity threshold used by imadjust
% (approximate values) as static texts
low_int_adj = double(min(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
high_int_adj = double(max(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
set(handles.high_int_threshold, 'String', num2str(high_int_adj));
set(handles.low_int_threshold, 'String', num2str(low_int_adj));

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

% Update the slider to accomodate all of the images
ntcells = sum(cellfun(@(x) length(x), handles.masks.cellList.meshData));
handles.cells_per_frame = cumsum(cellfun(@(x) length(x), handles.masks.cellList.meshData));
set(handles.frameslider, 'Min', 1, 'Max', ntcells, ...
    'SliderStep', [1 1]/(ntcells - 1), 'Value', 1)

% Choose default command line output for gui_spot_detection_indep_frames
handles.output = get(handles.frameslider, 'value');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_spot_detection_indep_frames wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_spot_detection_indep_frames_OutputFcn(~, ~, handles) 
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

% get the value of the current frame and cell
slideval = round(get(hObject, 'Value'));
frame = 1;
cell_no = slideval;
for i = 1:(length(handles.cells_per_frame)-1)
    if slideval > handles.cells_per_frame(i)
        frame = i+1;
        cell_no = slideval - handles.cells_per_frame(i);
    end
end

% add frame and cell_no to handle to be able to use those values also in the
% other functions
handles.frame = frame;
handles.cell_no = cell_no;

this_cell = handles.masks.cellList.meshData{frame}{cell_no};
this_area = this_cell.box; % [x, y, width, height]

% update the plot limits to avoid cropping the plot
handles.axes1.XLim = [0.5, this_area(3)+1.5];
handles.axes1.YLim = [0.5, this_area(4)+1.5];

image_PH = handles.imagedata{frame};
 
min_x = this_area(1) + handles.xshift;
max_x = this_area(1) + this_area(3) + handles.xshift;
min_y = this_area(2) + handles.yshift;
max_y = this_area(2) + this_area(4) + handles.yshift;
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

cla(handles.figure1) % clear the previous graph
imshow(imadjust(zoomed_fluor));
xlabel(['frame ' num2str(frame) ', cellid ' num2str(handles.masks.cellList.cellId{frame}(cell_no))]);
hold on

% write the lower and upper intensity threshold used by imadjust
% (approximate values) as static texts
low_int_adj = double(min(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
high_int_adj = double(max(zoomed_fluor(:)))/double(max(imadjust(zoomed_fluor(:))))*100;
set(handles.high_int_threshold, 'String', num2str(high_int_adj));
set(handles.low_int_threshold, 'String', num2str(low_int_adj));

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


% --- Executes on button press in delete_spot.
function delete_spot_Callback(hObject, eventdata, handles)
% hObject    handle to delete_spot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in next_spot.
function next_spot_Callback(hObject, eventdata, handles)
% hObject    handle to next_spot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function low_int_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to low_int_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of low_int_threshold as text
%        str2double(get(hObject,'String')) returns contents of low_int_threshold as a double
low_int_threshold = str2double(get(hObject,'String'))/100;
high_int_threshold = str2double(get(handles.high_int_threshold,'String'))/100;

% plot the data with the updated intensity thresholds
image_PH = handles.imagedata{handles.frame};

this_cell = handles.masks.cellList.meshData{handles.frame}{handles.cell_no};
this_area = this_cell.box;
 
min_x = this_area(1) + handles.xshift;
max_x = this_area(1) + this_area(3) + handles.xshift;
min_y = this_area(2) + handles.yshift;
max_y = this_area(2) + this_area(4) + handles.yshift;
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

% find the previous plot of the image and replace it with the new one
im_old = findobj(handles.axes1.Children,'type','image');
delete(im_old);
imshow(imadjust(zoomed_fluor,[low_int_threshold high_int_threshold],[]));
h = get(gca,'Children');
set(gca,'Children',[h(2:end); h(1)])

% --- Executes during object creation, after setting all properties.
function low_int_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to low_int_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function high_int_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to high_int_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of high_int_threshold as text
%        str2double(get(hObject,'String')) returns contents of high_int_threshold as a double
high_int_threshold = str2double(get(hObject,'String'))/100;
low_int_threshold = str2double(get(handles.low_int_threshold,'String'))/100;

% plot the data with the updated intensity thresholds
image_PH = handles.imagedata{handles.frame};

this_cell = handles.masks.cellList.meshData{handles.frame}{handles.cell_no};
this_area = this_cell.box;
 
min_x = this_area(1) + handles.xshift;
max_x = this_area(1) + this_area(3) + handles.xshift;
min_y = this_area(2) + handles.yshift;
max_y = this_area(2) + this_area(4) + handles.yshift;
zoomed_fluor = image_PH(max(1,min_y):min(size(image_PH,1),max_y),max(1,min_x):min(size(image_PH,2),max_x));

% find the previous plot of the image and replace it with the new one
im_old = findobj(handles.axes1.Children,'type','image');
delete(im_old);
imshow(imadjust(zoomed_fluor,[low_int_threshold high_int_threshold],[]));
h = get(gca,'Children');
set(gca,'Children',[h(2:end); h(1)])


% --- Executes during object creation, after setting all properties.
function high_int_threshold_CreateFcn(hObject, eventdata, handles)
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

this_cell = handles.masks.cellList.meshData{handles.frame}{handles.cell_no};
this_area = this_cell.box;
 
min_x = this_area(1) + handles.xshift;
max_x = this_area(1) + this_area(3) + handles.xshift;
min_y = this_area(2) + handles.yshift;
max_y = this_area(2) + this_area(4) + handles.yshift;
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
