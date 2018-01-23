function varargout = ViewCellTracks(varargin)
% VIEWCELLTRACKS MATLAB code for ViewCellTracks.fig
%      VIEWCELLTRACKS, by itself, creates a new VIEWCELLTRACKS or raises the existing
%      singleton*.
%
%      H = VIEWCELLTRACKS returns the handle to a new VIEWCELLTRACKS or the handle to
%      the existing singleton*.
%
%      VIEWCELLTRACKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWCELLTRACKS.M with the given input arguments.
%
%      VIEWCELLTRACKS('Property','Value',...) creates a new VIEWCELLTRACKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViewCellTracks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViewCellTracks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewCellTracks

% Last Modified by GUIDE v2.5 23-Jan-2018 13:51:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ViewCellTracks_OpeningFcn, ...
                   'gui_OutputFcn',  @ViewCellTracks_OutputFcn, ...
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


% --- Executes just before ViewCellTracks is made visible.
function ViewCellTracks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ViewCellTracks (see VARARGIN)

% Choose default command line output for ViewCellTracks
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ViewCellTracks wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ViewCellTracks_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


handles.currT = int32(get(handles.slider1, 'Value'));
% axes(handles.axes1);hold off
UpdateImage(handles,handles.currT);
xyt_curr = handles.fulltrack ;
colormap = handles.colormap;
handles.tracktoplot = str2double(get(handles.ReadTrackID, 'String'));
if  handles.currT<=size(xyt_curr,1)% not all cells are tracked all the way
hplot2 = plot(xyt_curr(handles.currT,1),xyt_curr(handles.currT,2)...
    ,'p','MarkerFaceColor',colormap(handles.trackcolor,:),'MarkerEdgeColor','k','MarkerSize',8,'LineWidth',1);
title(['Current time point  ' num2str(xyt_curr(handles.currT,3)*handles.delta_t/60)   'hrs since start; Frame(' num2str(xyt_curr(handles.currT,3)) ')']);
handles.residualdat = hplot2;
else % if the slider is pressed, and the track lost, only the iage will be shown
title(['This Track ID was lost. Current Frame(' num2str(handles.currT) ')']);
handles.residualdat = [];
end
if isfinite(handles.tracktoplot)
%get rid of the plot of the older track
delete(handles.oldlabel);
end
%handles.residualdat.Visible = 'off';
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function ReadTrackID_Callback(hObject, eventdata, handles)
% hObject    handle to ReadTrackID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ReadTrackID as text
%        str2double(get(hObject,'String')) returns contents of ReadTrackID as a double
handles.PlotTrackID.Visible = 'On';
handles.tracktoplot = str2double(get(handles.ReadTrackID, 'String'));
if handles.tracktoplot > handles.total_tracks
    disp(['There are ' num2str(handles.total_tracks) ' tracks in this image']);
    handles.PlotTrackID.Visible = 'Off';
    return
end
if handles.tracktoplot == 0
    disp('Track ID must be greater than 0');
    handles.PlotTrackID.Visible = 'Off';
    return
end
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function ReadTrackID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ReadTrackID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotTrackID.
function PlotTrackID_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTrackID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% TODO: here try to clear hobject or handles, since after couple tracks,
% gui gets very slow
colormap = prism;
randcolor = randi(size(colormap,1));% select the track label color once
xyt = handles.tracks;
handles.trackcolor = randcolor;
handles.colormap = prism;

if handles.tracktoplot > 0 
firstT = xyt(handles.tracktoplot).dat(1,3);
handles.tracktoplot = str2double(get(handles.ReadTrackID, 'String'));
handles.fulltrack = xyt(handles.tracktoplot).dat;
axes(handles.axes1); hold on
% anytime a new track id is chosen, set the slider to initial time point
set(handles.slider1, 'Value', 1); 
UpdateImage(handles,firstT);

if ~isempty(xyt(handles.tracktoplot).dat)
hdata  = plot(xyt(handles.tracktoplot).dat(firstT,1),xyt(handles.tracktoplot).dat(firstT,2)...
    ,'p','MarkerFaceColor',colormap(handles.trackcolor,:),'MarkerEdgeColor','k','MarkerSize',8,'LineWidth',1);hold on
title(['Current time point  ' num2str(xyt(handles.tracktoplot).dat(firstT,3)*handles.delta_t/60)   'hrs since start']);
handles.oldlabel = hdata;
else
    disp('Choose another trackID, this track is empty')
end
end
guidata(hObject, handles);

% --- Executes on button press in loadmaxprojections.
function loadmaxprojections_Callback(hObject, eventdata, handles)
% hObject    handle to loadmaxprojections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);hold off
cla reset;
handles.axes1.YTickLabel=[];
handles.axes1.XTickLabel=[];
box on
image.DirectoryProjections = { {'uigetdir(''.'')'} };
image.MicroscopePosition_Andor = 0;
image.channel = 0;
image.delta_t_in_minutes = 0;
%image.contrastLims=[0 0];
image = StructDlg(image);
handles.delta_t = image.delta_t_in_minutes;
handles.chan = image.channel;
handles.dir = image.DirectoryProjections;
handles.pos = image.MicroscopePosition_Andor;
%handles.contrastLims=image.contrastLims;% user will input initial contrast limits to show the image
handles.currT = 1;
%here load the correct channel max projections data
[multi_img,nt] = loadProjections(handles);
 handles.szt =nt;
 handles.allprojections =multi_img; 
% now display the images
imshow(multi_img(:,:,1),[0 3500],'Parent',handles.axes1);% todo: deal with contrast limits smarter
set(handles.slider1, 'Min', 1); 
set(handles.slider1, 'Max', nt);
set(handles.slider1, 'Value', 1); %init at timept 1
set(handles.slider1, 'SliderStep', [1 / (size(handles.allprojections,3)) ,...
    10 / (size(handles.allprojections,3)) ]);%
%here make the loadtracks button visible
handles.loadtrackingdata.Visible = 'On';
guidata(hObject, handles);

% --- Executes on button press in loadtrackingdata.
function loadtrackingdata_Callback(hObject, eventdata, handles)
% hObject    handle to loadtrackingdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tracks.Ilastikfile_TrackingData_h5 ={ {'uigetfile(''.'')'} };
%tracks.trackedcelltype =0;
tracks = StructDlg(tracks);
handles.ilastikfile = tracks.Ilastikfile_TrackingData_h5;
%handles.trackedcelltype= tracks.trackedcelltype ;
[coordintime] = extracktIlastikTracks(handles);
handles.tracks = coordintime;
handles.total_tracks = size(coordintime,2);
handles.DisplayTotalTracksField.Visible = 'On';
set(handles.DisplayTotalTracksField,'String',num2str(handles.total_tracks));
handles.TotalTracksStaticText.Visible = 'On';
handles.ReadTrackID.Visible = 'On';
handles.ClearButton.Visible = 'On';
guidata(hObject, handles);

function [multi_img,nt] = loadProjections(handles)
ff1 = readAndorDirectory(handles.dir); 
nucmoviefile = getAndorFileName(ff1,ff1.p(handles.pos),[],[],ff1.w(handles.chan));%getAndorFileName(files,pos,time,z,w)
nreader = bfGetReader(nucmoviefile);
nt = nreader.getSizeT;% get the number of frames
disp(['OPENED ' nucmoviefile ]);
proj = bfopen(nucmoviefile);
multi_img = zeros(size(proj{1}{1},1),size(proj{1}{1},2),nt);% initialize to image size and the second dim is the number of time points
for jj=1:nt
multi_img(:,:,jj) = (proj{1}{jj});
end

function UpdateImage(handles,time)
imshow(handles.allprojections(:,:,time),[0 3500],'Parent',handles.axes1);


function [coordintime] = extracktIlastikTracks(handles)
% here fully use Ilastik Automatic tracking 
clear coordintime
handles.pos = 1;%1
handles.chan = 1;% ff.w(chan): 0 - cfp cells; 1- nuc marker of other cell type
[coordintime] = process_ilastik_trackAN(handles.ilastikfile);

function DisplayTotalTracksField_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayTotalTracksField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DisplayTotalTracksField as text
%        str2double(get(hObject,'String')) returns contents of DisplayTotalTracksField as a double


% --- Executes during object creation, after setting all properties.
function DisplayTotalTracksField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplayTotalTracksField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in ClearButton.
function ClearButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);hold off
cla reset;
handles.axes1.YTickLabel=[];
handles.axes1.XTickLabel=[];
box on
set(handles.slider1, 'Value', 1); 
set(handles.ReadTrackID,'String','Track ID');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function TotalTracksStaticText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TotalTracksStaticText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% TODO: look into lagging after looking at couple tracks...seems that a
% little slower response, need to clear handles?
% TODO: maybe add the feature where you can input the good trackID and save
% it back for analysis
% TODO: other channel input (image) or other segmented cells input(plot
% other cell type), although this will only be relevan to sorting