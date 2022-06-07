function varargout = ortho_export(varargin)
% ORTHO_EXPORT M-file for ortho_export.fig
%      ORTHO_EXPORT, by itself, creates a new ORTHO_EXPORT or raises the existing
%      singleton*.
%
%      H = ORTHO_EXPORT returns the handle to a new ORTHO_EXPORT or the handle to
%      the existing singleton*.
%
%      ORTHO_EXPORT('Property','Value',...) creates a new ORTHO_EXPORT using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ortho_export_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ORTHO_EXPORT('CALLBACK') and ORTHO_EXPORT('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ORTHO_EXPORT.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help ortho_export

% Last Modified by GUIDE v2.5 26-Jul-2005 22:12:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ortho_export_OpeningFcn, ...
                   'gui_OutputFcn',  @ortho_export_OutputFcn, ...
                   'gui_LayoutFcn',  @ortho_export_LayoutFcn, ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ortho_export is made visible.
function ortho_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ortho_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% if exist('~/matlab/img/LuisOrtho.wav')
%     snd = wavread('~/matlab/img/LuisOrtho.wav');
%     sound(snd,44100);
% end
% UIWAIT makes ortho_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% these are the boolean variables to determine if a field has been filled
% in or an option chosen ...

global args

    args.ROIsize = 0;
    args.ROItype = 'cube';
    args.threshold = 0;
    args.onsets = [];
    args.window = 20;
    args.spm_file = [];
    args.spm_file2 = [];
    args.anat_file = [];
    args.tseries_path = [];
    args.tseries_file = [];
    args.tseries_file2 = [];
    args.doDetrend = 0;
    args.doGfilter = 0;
    args.doFFT = 0;
    args.ignore_origin = 0;
    args.wscale = [];
    args.interact = 1;
    args.xyz=[];
    args.mask_file = [];
    args.output_name = 'Ortho';
    args.voxFile = [];


% --- Outputs from this function are returned to the command line.
function varargout = ortho_export_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in GO_button.
function GO_button_Callback(hObject, eventdata, handles)
% hObject    handle to GO_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global args
args
ortho2005(args)

return

% --- Executes on button press in STOP_button.
function STOP_button_Callback(hObject, eventdata, handles)
% hObject    handle to STOP_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global myfig Evfig histFig FFTfig ACTIVE

ACTIVE=0;

if exist('myfig')
    close(myfig)
end
if exist('EvFig')
    close(Evfig)
end
if exist('histFig')
    close(histFig)
end
if exist('FFTfig')
    close(FFTfig)
end

return

% --- Executes during object creation, after setting all properties.
function ROItype_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROItype_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in ROItype_list.
function ROItype_list_Callback(hObject, eventdata, handles)
% hObject    handle to ROItype_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ROItype_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROItype_list
    global args
    ROItype = get(findobj('Tag','ROItype_list'), 'Value')
    switch ROItype
        case 1
            args.ROItype = 'cube'
        case 2
            args.ROItype = 'sphere'
            
        case 3
            args.ROItype = 'voxFile'
            
        case 4
            args.ROItype = 'maskFile'
    end
        

% --- Executes during object creation, after setting all properties.
function NN_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NN_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function NN_field_Callback(hObject, eventdata, handles)
% hObject    handle to NN_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NN_field as text
%        str2double(get(hObject,'String')) returns contents of NN_field as a double
global args
args.ROIsize = str2num(get(gco, 'String'))

% --- Executes on button press in AnatFile_button.
function AnatFile_button_Callback(hObject, eventdata, handles)
% hObject    handle to AnatFile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img','Select  *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','AnatFile_field');
    set(handle, 'String',name);
    args.anat_file = name;
return
    
    % --- Executes on button press in SM1_button.
function SM1_button_Callback(hObject, eventdata, handles)
% hObject    handle to SM1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','SM1_field');
    set(handle, 'String',name);
    args.spm_file = name;
return
 

% --- Executes on button press in SM2_button.
function SM2_button_Callback(hObject, eventdata, handles)
% hObject    handle to SM2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','SM2_field');
    set(handle, 'String',name);
    args.spm_file2 = name;
return

% --- Executes on button press in TS1_button.
function TS1_button_Callback(hObject, eventdata, handles)
% hObject    handle to TS1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','TS1_field');
    set(handle, 'String',name);
    args.tseries_file = name;
return

% --- Executes on button press in TS2_button.
function TS2_button_Callback(hObject, eventdata, handles)
% hObject    handle to TS2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','TS2_field');
    set(handle, 'String',name);
    args.tseries2_file = name;
return


% --- Executes on button press in TS3_button.
function TS3_button_Callback(hObject, eventdata, handles)
% hObject    handle to TS3_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','TS3_field');
    set(handle, 'String',name);

return


% --- Executes during object creation, after setting all properties.
function Th1_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Th1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Th1_field_Callback(hObject, eventdata, handles)
% hObject    handle to Th1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Th1_field as text
%        str2double(get(hObject,'String')) returns contents of Th1_field as a double
    global args
    str = get(findobj('Tag','Th1_field'), 'String');
    args.threshold = str2num(str);
return
    
% --- Executes during object creation, after setting all properties.
function Th2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Th2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Th2_field_Callback(hObject, eventdata, handles)
% hObject    handle to Th2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Th2_field as text
%        str2double(get(hObject,'String')) returns contents of Th2_field as a double


% --- Executes during object creation, after setting all properties.
function ons1_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ons1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ons1_field_Callback(hObject, eventdata, handles)
% hObject    handle to ons1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ons1_field as text
%        str2double(get(hObject,'String')) returns contents of ons1_field as a double
global args

str = get(hObject,'String')
if length(str)>0 
    args.onsets = str;
end
return

% --- Executes during object creation, after setting all properties.
function on2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to on2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --- Executes during object creation, after setting all properties.
function ons3_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ons3_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ons3_field_Callback(hObject, eventdata, handles)
% hObject    handle to ons3_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ons3_field as text
%        str2double(get(hObject,'String')) returns contents of ons3_field as a double
    

% --- Executes during object creation, after setting all properties.
function TW_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TW_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function TW_field_Callback(hObject, eventdata, handles)
% hObject    handle to TW_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TW_field as text
%        str2double(get(hObject,'String')) returns contents of TW_field as a double
   global args
    str = get(hObject,'String')
    args.window = str2num(str);
    return

% --- Executes during object creation, after setting all properties.
function anat_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anat_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function SM2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SM2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function SM2_field_Callback(hObject, eventdata, handles)
% hObject    handle to SM2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SM2_field as text
%        str2double(get(hObject,'String')) returns contents of SM2_field as a double
global args
args.spm_file2 = get(hObject,'String');

return

% --- Executes during object creation, after setting all properties.
function SM1_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SM1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function SM1_field_Callback(hObject, eventdata, handles)
% hObject    handle to SM1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SM1_field as text
%        str2double(get(hObject,'String')) returns contents of SM1_field as a double
global args
args.spm_file = get(hObject,'String');
return

% --- Executes during object creation, after setting all properties.
function TS3_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TS3_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function TS3_field_Callback(hObject, eventdata, handles)
% hObject    handle to TS3_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TS3_field as text
%        str2double(get(hObject,'String')) returns contents of TS3_field as a double



% --- Executes during object creation, after setting all properties.
function TS1_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TS1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function TS1_field_Callback(hObject, eventdata, handles)
% hObject    handle to TS1_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TS1_field as text
%        str2double(get(hObject,'String')) returns contents of TS1_field as a double
   global args

    args.tseries_file = get(hObject,'String');

return

% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes on button press in VF_button.
function VF_button_Callback(hObject, eventdata, handles)
% hObject    handle to VF_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global args
    [name path] = uigetfile('*','Select  file with ROI');
	name = strcat(path,name)
    handle = findobj('Tag','VF_field');
    set(handle, 'String',name);
    if args.ROItype=='mask_file'
        args.mask_file = name;
    else
        args.voxFile = name;
    end
    
return
    
% --- Executes during object creation, after setting all properties.
function VF_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VF_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function VF_field_Callback(hObject, eventdata, handles)
% hObject    handle to VF_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VF_field as text
%        str2double(get(hObject,'String')) returns contents of VF_field as a double
   global args
    str = get(hObject,'String')
    args.anat_file = str;

% --- Executes during object creation, after setting all properties.
function ons2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ons2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ons2_field_Callback(hObject, eventdata, handles)
% hObject    handle to ons2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ons2_field as text
%        str2double(get(hObject,'String')) returns contents of ons2_field as a double

% --- Executes during object creation, after setting all properties.
function AnatFile_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnatFile_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function AnatFile_field_Callback(hObject, eventdata, handles)
% hObject    handle to AnatFile_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AnatFile_field as text
%        str2double(get(hObject,'String')) returns contents of AnatFile_field as a double
    global args
    str = get(hObject,'String')
    args.anat_file = str;

return


% --- Executes during object creation, after setting all properties.
function TS2_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TS2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function TS2_field_Callback(hObject, eventdata, handles)
% hObject    handle to TS2_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TS2_field as text
%        str2double(get(hObject,'String')) returns contents of TS2_field as a double

return


% --- Executes on button press in detrend_cbx.
function detrend_cbx_Callback(hObject, eventdata, handles)
% hObject    handle to detrend_cbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detrend_cbx
global args
args.doDetrend = get(hObject,'Value');
return

% --- Executes on button press in GFilt_cbx.
function GFilt_cbx_Callback(hObject, eventdata, handles)
% hObject    handle to GFilt_cbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GFilt_cbx
global args
args.doGfilter = get(hObject,'Value');
return


% --- Executes on button press in FTbox.
function FTbox_Callback(hObject, eventdata, handles)
% hObject    handle to FTbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FTbox
global args
args.doFFT = get(hObject,'Value')
return



function output_name_Callback(hObject, eventdata, handles)
% hObject    handle to output_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_name as text
%        str2double(get(hObject,'String')) returns contents of output_name as a double
global args
args.output_name = get(hObject,'String');
return

% --- Executes during object creation, after setting all properties.
function output_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Creates and returns a handle to the GUI figure. 
function h1 = ortho_export_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end
load ortho_export.mat


appdata = [];
appdata.GUIDEOptions = mat{1};
appdata.UsedByGUIData_m = struct(...
    'figure1', [], ...
    'text1', [], ...
    'pushbutton1', 102.001098632812, ...
    'output', []);
appdata.lastValidTag = 'figure1';
appdata.GUIDELayoutEditor = [];

h1 = figure(...
'Units','characters',...
'Color',[0.702000021934509 0.702000021934509 0.702000021934509],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'DockControls','off',...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','ortho',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[108.6 20.8205128205128 117.833333333333 32.3333333333333],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[],...
'Behavior',get(0,'defaultfigureBehavior'),...
'Visible','on',...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'GO_button';

h2 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
'Callback','ortho_export(''GO_button_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'FontSize',10,...
'ListboxTop',0,...
'Position',[0.512022630834512 0.0283505154639175 0.199434229137199 0.134020618556701],...
'String','GO!',...
'Tag','GO_button',...
'UserData',[],...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'STOP_button';

h3 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
'Callback','ortho_export(''STOP_button_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'FontSize',10,...
'ListboxTop',0,...
'Position',[0.738330975954738 0.0309278350515464 0.185289957567185 0.134020618556701],...
'String','CLEAR DISPLAY',...
'Tag','STOP_button',...
'UserData',[],...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'text5';

h4 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'FontSize',10,...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[0.401697312588402 0.762886597938144 0.0848656294200849 0.0592783505154639],...
'String','Thresholds',...
'Style','text',...
'Tag','text5',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'NN_field';

h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''NN_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[107.333333333333 23.5833333333333 8.16666666666667 1.58333333333333],...
'String','0',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''NN_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','NN_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'text16';

h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'Position',[88.1666666666667 23.75 20.1666666666667 1.16666666666667],...
'String','ROI size (vox or mm)',...
'Style','text',...
'Tag','text16',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'text20';

h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'Position',[61.3333333333333 6.5 18.5 1.08333333333333],...
'String','Time Window',...
'Style','text',...
'Tag','text20',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'AnatFile_button';

h8 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.899999976158142 0.899999976158142 0.899999976158142],...
'Callback','ortho_export(''AnatFile_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[3.66666666666667 27.1666666666667 14.6666666666667 2],...
'String','Anatomical',...
'Tag','AnatFile_button',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'SM1_button';

h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 0.699999988079071 0.5],...
'Callback','ortho_export(''SM1_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[3.66666666666667 23.0833333333333 14.3333333333333 1.83333333333333],...
'String','Stats map 1',...
'Tag','SM1_button',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'SM2_button';

h10 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.67843137254902 0.92156862745098 1],...
'Callback','ortho_export(''SM2_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[3.66666666666667 20.75 14.1666666666667 2],...
'String','Stats map 2',...
'Tag','SM2_button',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'TS1_button';

h11 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.899999976158142 0.899999976158142 0.800000011920929],...
'Callback','ortho_export(''TS1_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[3.66666666666667 13.4166666666667 15.8333333333333 1.91666666666667],...
'String','Time Series 1',...
'Tag','TS1_button',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'TS2_button';

h12 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''TS2_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ForegroundColor',[0.400000005960464 0.400000005960464 0.400000005960464],...
'Position',[3.66666666666667 11.1666666666667 15.5 1.91666666666667],...
'String','Time Series 2',...
'Tag','TS2_button',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'TS3_button';

h13 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''TS3_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ForegroundColor',[0.400000005960464 0.400000005960464 0.400000005960464],...
'Position',[3.66666666666667 8.83333333333333 15.5 1.83333333333333],...
'String','Time Series 3',...
'Tag','TS3_button',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'Th1_field';

h14 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''Th1_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[48.1666666666667 23.5 7.5 1.41666666666667],...
'String','3.0',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''Th1_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','Th1_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'Th2_field';

h15 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''Th2_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[48 21.1666666666667 7.5 1.41666666666667],...
'String','3.0',...
'Style','edit',...
'Value',3,...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''Th2_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','Th2_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'ons1_field';

h16 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''ons1_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[79.5 13 31 1.91666666666667],...
'String','0',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''ons1_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','ons1_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'ons2_field';

h17 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''ons2_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ForegroundColor',[1 1 1],...
'Position',[79.5 10.8333333333333 31 1.91666666666667],...
'String','not yet',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''ons2_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','ons2_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'ons3_field';

h18 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''ons3_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ForegroundColor',[1 1 1],...
'Position',[79.5 9.08333333333333 31 1.91666666666667],...
'String','not yet',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''ons3_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','ons3_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'TW_field';

h19 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''TW_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[77.1666666666667 6.41666666666667 8.5 1.83333333333333],...
'String','20',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''TW_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','TW_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'AnatFile_field';

h20 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''AnatFile_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[22.1666666666667 27.1666666666667 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''AnatFile_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','AnatFile_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'SM2_field';

h21 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''SM2_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[22.1666666666667 20.75 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''SM2_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','SM2_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'SM1_field';

h22 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''SM1_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[22.1666666666667 23 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''SM1_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','SM1_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'TS3_field';

h23 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''TS3_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ForegroundColor',[1 1 1],...
'Position',[22.1666666666667 8.91666666666667 22.6666666666667 2],...
'String','not yet',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''TS3_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','TS3_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'TS1_field';

h24 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''TS1_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[22.1666666666667 13.5 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''TS1_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','TS1_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'TS2_field';

h25 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''TS2_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ForegroundColor',[1 1 1],...
'Position',[22.1666666666667 11.25 22.6666666666667 2],...
'String','not yet',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''TS2_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','TS2_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'VF_button';

h26 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''VF_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[88.1666666666667 27.6666666666667 26.8333333333333 1.75],...
'String','Mask Image / Vox. File ',...
'Tag','VF_button',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'VF_field';

h27 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''VF_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[88.1666666666667 25.5 27.3333333333333 1.75],...
'String','',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''VF_field_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','VF_field',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'detrend_cbx';

h28 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''detrend_cbx_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[39.5 6.33333333333333 10.5 1.91666666666667],...
'String','Detrend',...
'Style','checkbox',...
'Tag','detrend_cbx',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'GFilt_cbx';

h29 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''GFilt_cbx_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[22 6.33333333333333 13.6666666666667 2],...
'String','Lo-Pass Filt.',...
'Style','checkbox',...
'Tag','GFilt_cbx',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'text27';

h30 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'Position',[63.8333333333333 12.8333333333333 12.6666666666667 1.33333333333333],...
'String','Onsets 1',...
'Style','text',...
'Tag','text27',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'text28';

h31 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ForegroundColor',[0.400000005960464 0.400000005960464 0.400000005960464],...
'Position',[63.5 11 12.6666666666667 1.33333333333333],...
'String','Onsets 2',...
'Style','text',...
'Tag','text28',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'text29';

h32 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ForegroundColor',[0.400000005960464 0.400000005960464 0.400000005960464],...
'Position',[63.5 9.16666666666667 12.6666666666667 1.33333333333333],...
'String','Onsets 3',...
'Style','text',...
'Tag','text29',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel1';

h33 = uibuttongroup(...
'Parent',h1,...
'Units','characters',...
'FontSize',15,...
'Title','Time Series',...
'Tag','uipanel1',...
'Behavior',struct(),...
'Clipping','on',...
'Position',[1.5 5.91666666666667 56.8333333333333 11.5],...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'FTbox';

h34 = uicontrol(...
'Parent',h33,...
'Units','characters',...
'Callback','ortho_export(''FTbox_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[4.5 0.333333333333333 11 2.16666666666667],...
'String','Plot FFT',...
'Style','checkbox',...
'Tag','FTbox',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel2';

h35 = uibuttongroup(...
'Parent',h1,...
'Units','characters',...
'FontSize',15,...
'Title','Display',...
'Tag','uipanel2',...
'Behavior',struct(),...
'Clipping','on',...
'Position',[1.5 19.8333333333333 56.8333333333333 11.5],...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel3';

h36 = uibuttongroup(...
'Parent',h1,...
'Units','characters',...
'FontSize',15,...
'Title','ROI',...
'Tag','uipanel3',...
'Behavior',struct(),...
'Clipping','on',...
'Position',[59.5 22.5 56.8333333333333 8.83333333333333],...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel4';

h37 = uibuttongroup(...
'Parent',h1,...
'Units','characters',...
'FontSize',15,...
'Title','Timing (scan units)',...
'Tag','uipanel4',...
'Behavior',struct(),...
'Clipping','on',...
'Position',[59.3333333333333 5.91666666666667 57.3333333333333 11.5],...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'ROItype_list';

h38 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_export(''ROItype_list_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[61.5 23.8333333333333 25.6666666666667 5.66666666666667],...
'String',{  'Cube'; 'Thresholded sphere'; 'Voxel File'; 'Mask file' },...
'Style','listbox',...
'Value',1,...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''ROItype_list_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','ROItype_list',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'text30';

h39 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'Position',[60 19.6666666666667 17.5 1.41666666666667],...
'String','Output Name Prefix',...
'Style','text',...
'Tag','text30',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'output_name';

h40 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.699999988079071 0.699999988079071 0.699999988079071],...
'Callback','ortho_export(''output_name_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[79.3333333333333 19.6666666666667 33.1666666666667 1.41666666666667],...
'String','',...
'Style','edit',...
'CreateFcn', {@local_CreateFcn, 'ortho_export(''output_name_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','output_name',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'text31';

h41 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontAngle','italic',...
'FontSize',17,...
'FontWeight','bold',...
'ForegroundColor',[0.313725490196078 0.313725490196078 0.313725490196078],...
'Position',[0.833333333333333 0.75 44.3333333333333 2.25],...
'String','Ortho 2005 (july version)',...
'Style','text',...
'Tag','text31',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );


hsingleton = h1;


% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   eval(createfcn);
end


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)

gui_StateFields =  {'gui_Name'
                    'gui_Singleton'
                    'gui_OpeningFcn'
                    'gui_OutputFcn'
                    'gui_LayoutFcn'
                    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);        
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % ORTHO_EXPORT
    % create the GUI
    gui_Create = 1;
elseif isequal(ishandle(varargin{1}), 1) && ispc && iscom(varargin{1}) && isequal(varargin{1},gcbo)
    % ORTHO_EXPORT(ACTIVEX,...)    
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif ischar(varargin{1}) && numargin>1 && isequal(ishandle(varargin{2}), 1)
    % ORTHO_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % ORTHO_EXPORT(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = 1;
end

if gui_Create == 0
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else
        feval(varargin{:});
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.
    
    % Do feval on layout code in m-file if it exists
    if ~isempty(gui_State.gui_LayoutFcn)
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);
        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen')
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        end
    end
    
    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    
    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig 
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA. If there is
        % user set GUI data already, keep that also.
        data = guidata(gui_hFigure);
        handles = guihandles(gui_hFigure);
        if ~isempty(handles)
            if isempty(data)
                data = handles;
            else
                names = fieldnames(handles);
                for k=1:length(names)
                    data.(char(names(k)))=handles.(char(names(k)));
                end
            end
        end
        guidata(gui_hFigure, data);
    end
    
    % If user specified 'Visible','off' in p/v pairs, don't make the figure
    % visible.
    gui_MakeVisible = 1;
    for ind=1:2:length(varargin)
        if length(varargin) == ind
            break;
        end
        len1 = min(length('visible'),length(varargin{ind}));
        len2 = min(length('off'),length(varargin{ind+1}));
        if ischar(varargin{ind}) && ischar(varargin{ind+1}) && ...
                strncmpi(varargin{ind},'visible',len1) && len2 > 1
            if strncmpi(varargin{ind+1},'off',len2)
                gui_MakeVisible = 0;
            elseif strncmpi(varargin{ind+1},'on',len2)
                gui_MakeVisible = 1;
            end
        end
    end
    
    % Check for figure param value pairs
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end
        try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end
    
    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});
    
    if ishandle(gui_hFigure)
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
        
        % Make figure visible
        if gui_MakeVisible
            set(gui_hFigure, 'Visible', 'on')
            if gui_Options.singleton 
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        rmappdata(gui_hFigure,'InGUIInitialization');
    end
    
    % If handle visibility is set to 'callback', turn it on until finished with
    % OutputFcn
    if ishandle(gui_hFigure)
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end
    
    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end
    
    if ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end    

function gui_hFigure = local_openfig(name, singleton)

% this application data is used to indicate the running mode of a GUIDE
% GUI to distinguish it from the design mode of the GUI in GUIDE.
setappdata(0,'OpenGuiWhenRunning',1);

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
try 
    gui_hFigure = openfig(name, singleton, 'auto');
catch
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
end
rmappdata(0,'OpenGuiWhenRunning');

