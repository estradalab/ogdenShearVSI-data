function varargout = ortho_6(varargin)
% ORTHO_6 M-file for ortho_6.fig
%      ORTHO_6, by itself, creates a new ORTHO_6 or raises the existing
%      singleton*.
%
%      H = ORTHO_6 returns the handle to a new ORTHO_6 or the handle to
%      the existing singleton*.
%
%      ORTHO_6('Property','Value',...) creates a new ORTHO_6 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ortho_6_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ORTHO_6('CALLBACK') and ORTHO_6('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ORTHO_6.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help ortho_6

% Last Modified by GUIDE v2.5 04-Nov-2004 12:45:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ortho_6_OpeningFcn, ...
                   'gui_OutputFcn',  @ortho_6_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before ortho_6 is made visible.
function ortho_6_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ortho_6
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ortho_6 wait for user response (see UIRESUME)
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
function varargout = ortho_6_OutputFcn(hObject, eventdata, handles)
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
