function varargout = ShowStack(varargin)
% Function to visualize a 3D data stack. 
% Usage: ShowStack(matrix);
%        Input: A 3-D matrix of arrangement X,Y,Frames to scroll through.  
%        If input is a 4-D matrix,'Showstack' will average over the last dimension.
%
%       GUI inputs:
%       Colorrange is given by min and max color. 
%       The slider at the figure bottom is used to scroll through frames in the stack. 
%       Filter length determines the size of a filter that can be applied
%       to the image. The filter type is either a gaussian or a box filter
%       and can be chosen with the drop-out menu at the bottom left. 
%       If chosen filter is gaussian, the sigma of the kernel is determined
%       by 'Filter sigma'.
%       Smoothed frames can be saved to base workspace using the 'Save to
%       workspace' button in the top left button. Depending on filter
%       length, the smoothed frame will be slightly smaller as the original
%       frame because the code only uses the 'valid' part of the frame
%       without zero-padded edges.
%
% If you experience unexpected behavior of the code or have suggestions, 
% leave a comment and I will try to improve it in the future.
%
% Showstack was created, using GUIDE in Matlab2014b. 
% sm 5/8/2016

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ShowStack_OpeningFcn, ...
                   'gui_OutputFcn',  @ShowStack_OutputFcn, ...
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


% --- Executes just before ShowStack is made visible.
function ShowStack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ShowStack (see VARARGIN)

% Choose default command line output for ShowStack
handles.output = hObject;

if isempty(varargin)
    error('No input data found. Data should be a 3- or 4-D frame stack')
end

if length(size(varargin{1})) == 4
    handles.UserData.Frames = squeeze(nanmean(varargin{1},4)); %4-D array. Should be X,Y,Frames,Trials and will be averaged over last dimension.
elseif length(size(varargin{1})) == 3
    handles.UserData.Frames = varargin{1}; %3-D array. Should be X,Y,Frames/Trials
else
    error('Wrong input. Data should be a 3- or 4-D frame stack')
end

fLength = round(str2double(get(handles.FrameSmth,'string'))); %determines degree of smoothing.
data = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,1)),fLength,str2double(handles.FrameSigma.String),get(handles.FilterType,'value')); %get first frame to display
imshow(data,[str2double(get(handles.MinColor,'String')) str2double(get(handles.MaxColor,'String'))],'parent',handles.FrameImage); %create image object for preview
axis(handles.FrameImage,'square');colormap(handles.FrameImage,'jet');colorbar;
set(handles.CurrentFrame,'string',['Showing Frame 1 / ' num2str(size(handles.UserData.Frames,3))])

set(handles.FrameSlider,'Min',1);
set(handles.FrameSlider,'Max',size(handles.UserData.Frames,3));
set(handles.FrameSlider,'Value',1);
set(handles.FrameSlider,'SliderStep',[1/(size(handles.UserData.Frames,3)-1) 1/(size(handles.UserData.Frames,3)-1)]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ShowStack wait for user response (see UIRESUME)
% uiwait(handles.ShowStack);


% --- Outputs from this function are returned to the command line.
function varargout = ShowStack_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function FrameSlider_Callback(hObject, eventdata, handles)
% hObject    handle to FrameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(hObject,'Value',round(get(hObject,'Value')));
fLength = round(str2double(get(handles.FrameSmth,'string'))); %determines degree of smoothing.
if fLength == 1 %no smooth
    data = squeeze(handles.UserData.Frames(:,:,get(hObject,'Value'))); %get current frame to display
else
    data = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,get(hObject,'Value'))),fLength,str2double(handles.FrameSigma.String),get(handles.FilterType,'value')); %smooth current frame
end

imshow(data,[str2double(get(handles.MinColor,'String')) str2double(get(handles.MaxColor,'String'))],'parent',handles.FrameImage); %create image object for preview
axis(handles.FrameImage,'square');colormap(handles.FrameImage,'jet');colorbar;
set(handles.CurrentFrame,'string',['Showing Frame ' num2str(get(hObject,'Value')) ' / ' num2str(size(handles.UserData.Frames,3))])

% --- Executes during object creation, after setting all properties.
function FrameSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function MaxColor_Callback(hObject, eventdata, handles)
% hObject    handle to MaxColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxColor as text
%        str2double(get(hObject,'String')) returns contents of MaxColor as a double

if isnan(str2double(get(hObject,'string')))
   set(hObject,'string','0.01')
end
caxis(handles.FrameImage,[str2double(get(handles.MinColor,'String')) str2double(get(handles.MaxColor,'String'))]);

% --- Executes during object creation, after setting all properties.
function MaxColor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MinColor_Callback(hObject, eventdata, handles)
% hObject    handle to MinColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinColor as text
%        str2double(get(hObject,'String')) returns contents of MinColor as a double

if isnan(str2double(get(hObject,'string')))
   set(hObject,'string','0')
end
caxis(handles.FrameImage,[str2double(get(handles.MinColor,'String')) str2double(get(handles.MaxColor,'String'))]);


% --- Executes during object creation, after setting all properties.
function MinColor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FrameSmth_Callback(hObject, eventdata, handles)
% hObject    handle to FrameSmth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameSmth as text
%        str2double(get(hObject,'String')) returns contents of FrameSmth as a double
if isnan(str2double(get(hObject,'string')))
    set(hObject,'Value',round(get(hObject,'Value'))); %make sure value is an integer
    bin = round(str2double(get(handles.FrameSmth,'string'))); %determines degree of smoothing.
    data = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,get(hObject,'Value'))),fLength,str2double(handles.FrameSigma.String),get(handles.FilterType,'value')); %smooth current frame
    imshow(data,[str2double(get(handles.MinColor,'String')) str2double(get(handles.MaxColor,'String'))],'parent',handles.FrameImage); %create image object for preview
    axis(handles.FrameImage,'square');colormap(handles.FrameImage,'jet');colorbar;
    set(handles.CurrentFrame,'string',['Showing Frame ' num2str(get(hObject,'Value')) ' / ' num2str(size(handles.UserData.Frames,3))])

   set(hObject,'string','0')
end

% --- Executes during object creation, after setting all properties.
function FrameSmth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameSmth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in FilterType.
function FilterType_Callback(hObject, eventdata, handles)
% hObject    handle to FilterType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FilterType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FilterType


% --- Executes during object creation, after setting all properties.
function FilterType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilterType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FrameSigma_Callback(hObject, eventdata, handles)
% hObject    handle to FrameSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameSigma as text
%        str2double(get(hObject,'String')) returns contents of FrameSigma as a double

if isnan(str2double(get(hObject,'string')))
    set(hObject,'string','1.76');
    disp([get(hObject,'string') ' is not a valid input for sigma. Set back to default (1.76)'])
end

% --- Executes during object creation, after setting all properties.
function FrameSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SaveFrame.
function SaveFrame_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cFrame = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,get(hObject,'Value'))),round(str2double(get(handles.FrameSmth,'string'))),str2double(handles.FrameSigma.String),get(handles.FilterType,'value')); %smooth current frame
assignin('base','SavedFrame',cFrame);


function DataOut = ApplyFilter2(DataIn,fLength,sigma,type)
% short code to apply either gaussian or box filter. 
% Usage: DataOut = ApplyGauss(DataIn,fLength,sigma,type)
% Input:    DataIn = 2D matrix on which the filter should be applied
%           fLength = filter length.
%           sigma = sigma of the gaussian kernel. Standard is 1.76 for a single standard deviation.
%           type: type of filter. 1 for gaussian, 2 for box filter
% Output:   DataOut = filtered version of DataIn.

if fLength > 1
    if type == 1
        [x,y]=meshgrid(-fLength:fLength,-fLength:fLength); % create gaussian filter based on flength and sigma
        Kernel= exp(-(x.^2+y.^2)/(2*sigma*sigma))/(2*pi*sigma*sigma);
    elseif type == 2
        Kernel = ones(fLength,fLength); % Create box filter based on flength
        Kernel = Kernel ./ fLength^2;
    end
    DataOut = conv2(double(DataIn), Kernel,'valid'); %convolve original matrix with filter
else
    DataOut = DataIn; %don't apply filter if fLength is equal or lower then 1
end
