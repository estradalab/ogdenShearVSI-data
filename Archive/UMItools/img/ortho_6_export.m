function varargout = ortho_6_export(varargin)
% ORTHO_6_EXPORT M-file for ortho_6_export.fig
%      ORTHO_6_EXPORT, by itself, creates a new ORTHO_6_EXPORT or raises the existing
%      singleton*.
%
%      H = ORTHO_6_EXPORT returns the handle to a new ORTHO_6_EXPORT or the handle to
%      the existing singleton*.
%
%      ORTHO_6_EXPORT('Property','Value',...) creates a new ORTHO_6_EXPORT using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ortho_6_export_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ORTHO_6_EXPORT('CALLBACK') and ORTHO_6_EXPORT('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ORTHO_6_EXPORT.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help ortho_6_export

% Last Modified by GUIDE v2.5 02-Feb-2005 10:16:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ortho_6_export_OpeningFcn, ...
                   'gui_OutputFcn',  @ortho_6_export_OutputFcn, ...
                   'gui_LayoutFcn',  @ortho_6_export_LayoutFcn, ...
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


% --- Executes just before ortho_6_export is made visible.
function ortho_6_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ortho_6_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ortho_6_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% these are the boolean variables to determine if a field has been filled
% in or an option chosen ...
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
global doGfilter doDetrend doFFT
bAnat=0; bSM1=0; bSM2=0; 
bTS1=0; bTS2=0; bTS3=0; bVF=0; bOS1=0; bOS2=0; bOS3=0; bROI=0;
doGfilter = 0;
doDetrend = 0;
doFFT = 0;
warning off


% --- Outputs from this function are returned to the command line.
function varargout = ortho_6_export_OutputFcn(hObject, eventdata, handles)
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

global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI doFFT
disp('*Flags = bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 doFFT')

Flags =       [bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 doFFT]

fprintf('\n***** ******* *******')
if     Flags ==   [1    0    0    0    0    0    0   0    0   0  0]
    str=get(findobj('Tag','NN_field'), 'String')
    n=str2num(str);
    str=get(findobj('Tag','AnatFile_field'), 'String')
    
    fprintf('\n\n Calling orthov(%d,%s)...\n',n,str);
    fprintf('\n\n Exit by clicking outside of the plots');
    fprintf('\n');
    orthov(n,str);
    
elseif Flags ==   [1    1    0    0    0    0    0   0    0   0  0]
    str=get(findobj('Tag','Th1_field'), 'String')
    th=str2num(str);
    anat_file=get(findobj('Tag','AnatFile_field'), 'String')
    spm_file=get(findobj('Tag','SM1_field'), 'String')
    
    fprintf('\n\n Calling orthospm(%f,%s , %s)...',th,spm_file, anat_file);
    fprintf('\n Exit by clicking outside of the plots');
    fprintf('\n');
    orthospm(th,spm_file, anat_file);
    
elseif Flags ==   [0    0    0    1    0    0    0   0    0   0  0]
    str=get(findobj('Tag','NN_field'), 'String')
    n=str2num(str);
    
    str=get(findobj('Tag','TS1_field'), 'String')
    pos = find(str=='/');
    if isempty(pos)
        TSpath=pwd
        TSname=str
    else
        pos=pos(end);
        TSpath = str(1:pos-1)
        TSname = str(pos+1:end)
    end
    
    fprintf('\nCalling orthov2(%d,%s)...\n',n,str);
    fprintf('\n\n Exit by clicking outside of the plots');
    fprintf('\n RIGHT mouse click saves the time series');
    fprintf('\n');
    
    curDir=pwd;
    cd(TSpath);
    orthov2(n,str);
    cd(curDir);
    
elseif Flags ==   [1    0    0    1    0    0    0   0    0   0   0]
    str=get(findobj('Tag','NN_field'), 'String')
    n=str2num(str);
    str=get(findobj('Tag','TS1_field'), 'String')
    pos = find(str=='/');
    if isempty(pos)
        TSpath=pwd
        TSname=str
    else
        pos=pos(end);
        TSpath = str(1:pos-1)
        TSname = str(pos+1:end)
    end
    str2=get(findobj('Tag','AnatFile_field'), 'String')
    
    fprintf('\n Calling orthov2(%d,%s, %s)...\n',n,str,str2);
    fprintf('\nn Exit by clicking outside of the plots');
    fprintf('\n RIGHT mouse click saves the time series');
    fprintf('\n');
    
    curDir=pwd;
    cd(TSpath);
    orthov2(n,str, str2);
    cd(curDir);
    
elseif Flags ==   [0    0    0    1    0    0    0   0    0   0   1]
    str=get(findobj('Tag','NN_field'), 'String')
    n=str2num(str);
    str=get(findobj('Tag','TS1_field'), 'String')
    str=str(1:end-4);
    
    fprintf('\n Calling orthoft(%d,%s)...\n',n,str);
    fprintf('\n \nExit by clicking outside of the plots');
    fprintf('\n RIGHT mouse click saves the time series');
    fprintf('\n');
    
    pos = find(str=='/');
    if isempty(pos)
        TSpath=pwd
        TSname=str
    else
        pos=pos(end);
        TSpath = str(1:pos-1)
        TSname = str(pos+1:end)
    end
    
    curDir=pwd;
    cd(TSpath);
    orthoft(n,TSname);   
    cd(curDir);

elseif Flags ==   [1    1    0    1    0    0    0   1    0   0   0]
    str =get(findobj('Tag','NN_field'), 'String');
    n = str2num(str)
    str = get(findobj('Tag','Th1_field'), 'String');
    th = str2num(str)
    str = get(findobj('Tag','TS1_field'), 'String');
    pos = find(str=='/');
    if isempty(pos)
        TSpath=pwd
        TSname=str
    else
        pos=pos(end);
        TSpath = str(1:pos-1)
        TSname = str(pos+1:end)
    end
    ANAname = get(findobj('Tag','AnatFile_field'), 'String')
    SMname = get(findobj('Tag','SM1_field'), 'String')
    str = get(findobj('Tag','ons1_field'), 'String');
    eval(sprintf('ons = [%s]', str));
    str = get(findobj('Tag','TW_field'), 'String');
    window = str2num(str)
    ROItype = get(findobj('Tag','ROItype_list'), 'Value')
    switch ROItype
        case 1
            fprintf('\n Calling orthospm3(%d, %f, ons, %d, %s, %s, %s, %s)', ...
                n,th,ons,window, SMname,  ANAname, TSpath, TSname(1:end-8));
            fprintf('\n\n Exit by clicking outside of the plots');
            fprintf('\n RIGHT mouse click saves the time series');
            fprintf('\n\n');
            
            orthospm3(n,th,ons,window, SMname,  ANAname, TSpath, TSname(1:end-8));
            
        case 2
            fprintf('\n  Calling orthospm4b(%d, %f, ons, %d, %s, %s, %s, %s)',...
                n,th,ons,window, SMname,  ANAname, TSpath, TSname(1:end-8));    fprintf('\nExit by clicking outside of the plots');
            fprintf('\n\n Exit by clicking outside of the plots');
            fprintf('\n Use the RIGHT mouse button to look at the time series.  \nLeft button to move though the anatomy')          
            fprintf('\n\n');
            
            orthospm4b(n,th,ons,window, SMname,  ANAname, TSpath, TSname);
            
        case 3
            errormesg('Not implemented yet!')
    end
        
elseif Flags ==   [1    1    0    1    0    0    0   0    0   0   0]
    str =get(findobj('Tag','NN_field'), 'String');
    n = str2num(str)
    str = get(findobj('Tag','Th1_field'), 'String');
    th = str2num(str)
    str = get(findobj('Tag','TS1_field'), 'String');
    pos = find(str=='/');
    if isempty(pos)
        TSpath=pwd
        TSname=str
    else
        pos=pos(end);
        TSpath = str(1:pos-1)
        TSname = str(pos+1:end)
    end
    ANAname = get(findobj('Tag','AnatFile_field'), 'String')
    SMname = get(findobj('Tag','SM1_field'), 'String')
    
    fprintf('\nCalling orthospm2(%d, %f,  %s, %s, %s, %s)',...
        n,th, SMname,  ANAname, TSpath, TSname(1:end-8));
    fprintf('\n\n Exit by clicking outside of the plots');
    fprintf('\n RIGHT mouse click saves the time series');
    fprintf('\n');
    orthospm2(n,th, SMname,  ANAname, TSpath, TSname);
    

elseif Flags ==   [1    1    1    0    0    0    0   0    0   0  0]
    str=get(findobj('Tag','Th1_field'), 'String')
    th=str2num(str);
    str=get(findobj('Tag','Th2_field'), 'String')
    th2=str2num(str);

    anat_file=get(findobj('Tag','AnatFile_field'), 'String')
   
    spm_file=get(findobj('Tag','SM1_field'), 'String')
    spm_file2=get(findobj('Tag','SM2_field'), 'String')
    
    fprintf('\nCalling ovm(%f, %f, %s, %s, %s)',...
        th,th2, spm_file,spm_file2, anat_file);
    fprintf('\n\n Exit by clicking outside of the plots');
    fprintf('\n RIGHT mouse click saves the time series');
    fprintf('\n');
    
    ovm(th,th2, spm_file,spm_file2, anat_file);
    
else
    msgbox('This combination of inputs is not implemented')
end

return

% --- Executes on button press in STOP_button.
function STOP_button_Callback(hObject, eventdata, handles)
% hObject    handle to STOP_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
global myfig mainfig ftfig ACTIVE
ACTIVE=0;
close(myfig)
close(mainfig)
close(ftfig)
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


% --- Executes on button press in AnatFile_button.
function AnatFile_button_Callback(hObject, eventdata, handles)
% hObject    handle to AnatFile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
    [name path] = uigetfile('*.img','Select  *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','AnatFile_field');
    set(handle, 'String',name);
    bAnat = 1;
return
    
    % --- Executes on button press in SM1_button.
function SM1_button_Callback(hObject, eventdata, handles)
% hObject    handle to SM1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','SM1_field');
    set(handle, 'String',name);
    bSM1 = 1;
return
 

% --- Executes on button press in SM2_button.
function SM2_button_Callback(hObject, eventdata, handles)
% hObject    handle to SM2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','SM2_field');
    set(handle, 'String',name);
    bSM2 = 1;
return

% --- Executes on button press in TS1_button.
function TS1_button_Callback(hObject, eventdata, handles)
% hObject    handle to TS1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','TS1_field');
    set(handle, 'String',name);
    bTS1 = 1;
return

% --- Executes on button press in TS2_button.
function TS2_button_Callback(hObject, eventdata, handles)
% hObject    handle to TS2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','TS2_field');
    set(handle, 'String',name);
    bTS2 = 1;
return


% --- Executes on button press in TS3_button.
function TS3_button_Callback(hObject, eventdata, handles)
% hObject    handle to TS3_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
    [name path] = uigetfile('*.img','Select statistical *.img file');
	name = strcat(path,name)
    handle = findobj('Tag','TS3_field');
    set(handle, 'String',name);
    bTS3 = 1;
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
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
str = get(hObject,'String')
if length(str)>0 
    bOS1= 1;
else 
    bOS1=0;
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
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
str = get(hObject,'String')
if length(str)==0
    bSM2=0;
else
    bSM2=1;
end
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
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
str = get(hObject,'String')
if length(str)==0
    bSM1=0;
else
    bSM1=1;
end
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
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
bTS3=1;
return


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
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
str = get(hObject,'String');
if length(str)==0
    bTS1=0;
else
    bTS1=1;
end
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
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
bVF=1;
return


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
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
str = get(hObject,'String');
if length(str)>0
    bOS2 = 1;
else
    bOS2 = 0;
end
return

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
global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI
str = get(hObject,'String')
if length(str)>0 
    bAnat= 1;
else 
    bAnat=0;
end

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

global bAnat bSM1 bSM2 bTS1 bTS2 bTS3 bVF bOS1 bOS2 bOS3 bROI

str = get(hObject,'String')
if length(str)>0 , 
    bTS2=1;
else 
    bTS2=0;
end
return


% --- Executes on button press in detrend_cbx.
function detrend_cbx_Callback(hObject, eventdata, handles)
% hObject    handle to detrend_cbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detrend_cbx
global doDetrend
doDetrend = get(hObject,'Value');
return

% --- Executes on button press in GFilt_cbx.
function GFilt_cbx_Callback(hObject, eventdata, handles)
% hObject    handle to GFilt_cbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GFilt_cbx
global doGfilter
doGfilter = get(hObject,'Value');
return


% --- Executes on button press in FTbox.
function FTbox_Callback(hObject, eventdata, handles)
% hObject    handle to FTbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FTbox
global doFFT
doFFT = get(hObject,'Value')
return


% --- Creates and returns a handle to the GUI figure. 
function h1 = ortho_6_export_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end

h1 = figure(...
'Units','characters',...
'Color',[0.925490196078431 0.913725490196078 0.847058823529412],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','ortho_6',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',[108.6 25.2371794871795 117.666666666667 27.9166666666667],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[]);

setappdata(h1, 'GUIDEOptions',struct(...
'active_h', [], ...
'taginfo', struct(...
'figure', 2, ...
'pushbutton', 9, ...
'text', 30, ...
'edit', 21, ...
'frame', 5, ...
'radiobutton', 3, ...
'listbox', 2, ...
'togglebutton', 3, ...
'checkbox', 10), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'blocking', 0, ...
'lastSavedFile', 'C:\Documents and Settings\fmrilab\Desktop\hernan\matlab\ortho_6.m'));

setappdata(h1, 'UsedByGUIData_m',struct(...
'figure1', 100.00048828125, ...
'text1', 4.0008544921875, ...
'pushbutton1', 102.001098632813, ...
'output', 100.00048828125));


h2 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[0.514164305949009 0.66865671641791 0.456090651558074 0.241791044776119],...
'String',{  '' },...
'Style','frame',...
'Tag','frame3');


h3 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[0.0169971671388102 0.0925373134328358 0.485835694050992 0.817910447761194],...
'String',{  '' },...
'Style','frame',...
'Tag','frame2');


h4 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
'Callback','ortho_6_export(''GO_button_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'FontSize',10,...
'ListboxTop',0,...
'Position',[0.529745042492918 0.519402985074627 0.198300283286119 0.134328358208955],...
'String','GO!',...
'Tag','GO_button',...
'UserData',[]);


h5 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'BackgroundColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
'Callback','ortho_6_export(''STOP_button_Callback'',gcbo,[],guidata(gcbo))',...
'CData',[],...
'FontSize',10,...
'ListboxTop',0,...
'Position',[0.777620396600567 0.519402985074627 0.18413597733711 0.134328358208955],...
'String','CLEAR DISPLAY',...
'Tag','STOP_button',...
'UserData',[]);


h6 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'FontSize',10,...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[0.41643059490085 0.68955223880597 0.084985835694051 0.0597014925373134],...
'String','Thresholds',...
'Style','text',...
'Tag','text5');


h7 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'FontSize',15,...
'ListboxTop',0,...
'Position',[0.0410764872521247 0.871641791044776 0.155807365439093 0.0716417910447761],...
'String','Image Files',...
'Style','text',...
'Tag','text13');


h8 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'FontSize',15,...
'ListboxTop',0,...
'Position',[0.569405099150142 0.868656716417911 0.0977337110481586 0.0597014925373134],...
'String','ROI',...
'Style','text',...
'Tag','text14');


h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''ROItype_list_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'Position',[61.5 19.4166666666667 25.5 3.5],...
'String',{  'Num Neighbors'; 'Num Neighbors over threshold'; 'Voxel File' },...
'Style','listbox',...
'Value',1,...
'CreateFcn','ortho_6_export(''ROItype_list_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','ROItype_list');


h10 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[64 22.8333333333333 14.8333333333333 1],...
'String','Type of ROI',...
'Style','text',...
'Tag','text15');


h11 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''NN_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[87.6666666666667 20.5833333333333 8.5 1.83333333333333],...
'String','0',...
'Style','edit',...
'CreateFcn','ortho_6_export(''NN_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','NN_field');


h12 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[83 22.8333333333333 14 1.08333333333333],...
'String','# Neighbors',...
'Style','text',...
'Tag','text16');


h13 = uicontrol(...
'Parent',h1,...
'Units','normalized',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[0.514164305949009 0.0955223880597015 0.453257790368272 0.385074626865672],...
'String',{  '' },...
'Style','frame',...
'Tag','frame4');


h14 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',15,...
'ListboxTop',0,...
'Position',[66.5 12.4166666666667 24.8333333333333 1.41666666666667],...
'String','Timing (scan units)',...
'Style','text',...
'Tag','text17');


h15 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[61.3333333333333 3.58333333333333 18.5 1.08333333333333],...
'String','Time Window',...
'Style','text',...
'Tag','text20');


h16 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.899999976158142 0.899999976158142 0.899999976158142],...
'Callback','ortho_6_export(''AnatFile_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[3.66666666666667 21.6666666666667 14.6666666666667 2],...
'String','Anatomical',...
'Tag','AnatFile_button');


h17 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 0.699999988079071 0.5],...
'Callback','ortho_6_export(''SM1_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[3.66666666666667 17.5833333333333 14.3333333333333 1.83333333333333],...
'String','Stats map 1',...
'Tag','SM1_button');


h18 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.5 0.5 0.701960784313725],...
'Callback','ortho_6_export(''SM2_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[3.66666666666667 15.25 14.1666666666667 2],...
'String','Stats map 2',...
'Tag','SM2_button');


h19 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.899999976158142 0.899999976158142 0.800000011920929],...
'Callback','ortho_6_export(''TS1_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[3.66666666666667 10.0833333333333 15.8333333333333 1.91666666666667],...
'String','Time Series 1',...
'Tag','TS1_button');


h20 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_6_export(''TS2_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ForegroundColor',[0.400000005960464 0.400000005960464 0.400000005960464],...
'ListboxTop',0,...
'Position',[3.66666666666667 7.91666666666667 15.5 1.91666666666667],...
'String','Time Series 2',...
'Tag','TS2_button');


h21 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_6_export(''TS3_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ForegroundColor',[0.400000005960464 0.400000005960464 0.400000005960464],...
'ListboxTop',0,...
'Position',[3.66666666666667 5.83333333333333 15.5 1.83333333333333],...
'String','Time Series 3',...
'Tag','TS3_button');


h22 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''Th1_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[49.8333333333333 18 7.5 1.41666666666667],...
'String','3.0',...
'Style','edit',...
'CreateFcn','ortho_6_export(''Th1_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','Th1_field');


h23 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''Th2_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[49.6666666666667 15.8333333333333 7.5 1.41666666666667],...
'String','3.0',...
'Style','edit',...
'Value',3,...
'CreateFcn','ortho_6_export(''Th2_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','Th2_field');


h24 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''ons1_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[79.5 10.0833333333333 31 1.91666666666667],...
'String','0',...
'Style','edit',...
'CreateFcn','ortho_6_export(''ons1_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','ons1_field');


h25 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''ons2_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[79.5 7.91666666666667 31 1.91666666666667],...
'String','0',...
'Style','edit',...
'CreateFcn','ortho_6_export(''ons2_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','ons2_field');


h26 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''ons3_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[79.5 6.08333333333333 31 1.91666666666667],...
'String','0',...
'Style','edit',...
'CreateFcn','ortho_6_export(''ons3_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','ons3_field');


h27 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''TW_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[77.1666666666667 3.41666666666667 8.5 1.83333333333333],...
'String','20',...
'Style','edit',...
'CreateFcn','ortho_6_export(''TW_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','TW_field');


h28 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''AnatFile_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[22.1666666666667 21.6666666666667 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn','ortho_6_export(''AnatFile_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','AnatFile_field');


h29 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''SM2_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[22.1666666666667 15.25 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn','ortho_6_export(''SM2_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','SM2_field');


h30 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''SM1_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[22.1666666666667 17.5 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn','ortho_6_export(''SM1_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','SM1_field');


h31 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''TS3_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[22.1666666666667 5.91666666666667 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn','ortho_6_export(''TS3_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','TS3_field');


h32 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''TS1_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[22.1666666666667 10.0833333333333 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn','ortho_6_export(''TS1_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','TS1_field');


h33 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''TS2_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[22.1666666666667 8 22.6666666666667 2],...
'String','',...
'Style','edit',...
'CreateFcn','ortho_6_export(''TS2_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','TS2_field');


h34 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_6_export(''VF_button_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[99.8333333333333 22.8333333333333 11 1.83333333333333],...
'String','Vox. File',...
'Tag','VF_button');


h35 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','ortho_6_export(''VF_field_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[96.8333333333333 20.5 16.8333333333333 2],...
'String','',...
'Style','edit',...
'CreateFcn','ortho_6_export(''VF_field_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','VF_field');


h36 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_6_export(''detrend_cbx_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[102.833333333333 3.41666666666667 10.5 1.91666666666667],...
'String','Detrend',...
'Style','checkbox',...
'Tag','detrend_cbx');


h37 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_6_export(''GFilt_cbx_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[88.1666666666667 3.25 13.6666666666667 2],...
'String','Lo-Pass Filt.',...
'Style','checkbox',...
'Tag','GFilt_cbx');


h38 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','ortho_6_export(''FTbox_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[46 10.0833333333333 11 2.16666666666667],...
'String','Plot FFT',...
'Style','checkbox',...
'Tag','FTbox');


h39 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[63.8333333333333 9.83333333333333 12.6666666666667 1.33333333333333],...
'String','Onsets 1',...
'Style','text',...
'Tag','text27');


h40 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ForegroundColor',[0.400000005960464 0.400000005960464 0.400000005960464],...
'ListboxTop',0,...
'Position',[63.5 7.91666666666667 12.6666666666667 1.33333333333333],...
'String','Onsets 2',...
'Style','text',...
'Tag','text28');


h41 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ForegroundColor',[0.400000005960464 0.400000005960464 0.400000005960464],...
'ListboxTop',0,...
'Position',[63.5 6.08333333333333 12.6666666666667 1.33333333333333],...
'String','Onsets 3',...
'Style','text',...
'Tag','text29');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      ORTHO_6_EXPORT, by itself, creates a new ORTHO_6_EXPORT or raises the existing
%      singleton*.
%
%      H = ORTHO_6_EXPORT returns the handle to a new ORTHO_6_EXPORT or the handle to
%      the existing singleton*.
%
%      ORTHO_6_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ORTHO_6_EXPORT.M with the given input arguments.
%
%      ORTHO_6_EXPORT('Property','Value',...) creates a new ORTHO_6_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.5 $ $Date: 2003/07/17 18:28:28 $

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
        gui_Mfile = [getfield(gui_State, gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % ORTHO_6_EXPORT
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % ORTHO_6_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % ORTHO_6_EXPORT(...)
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

        % Generate HANDLES structure and store with GUIDATA
        guidata(gui_hFigure, guihandles(gui_hFigure));
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
        if ischar(varargin{ind}) & ischar(varargin{ind+1}) & ...
                strncmpi(varargin{ind},'visible',len1) & len2 > 1
            if strncmpi(varargin{ind+1},'off',len2)
                gui_MakeVisible = 0;
            elseif strncmpi(varargin{ind+1},'on',len2)
                gui_MakeVisible = 1;
            end
        end
    end
    
    % Check for figure param value pairs
    for index=1:2:length(varargin)
        if length(varargin) == index
            break;
        end
        try, set(gui_hFigure, varargin{index}, varargin{index+1}), catch, break, end
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

