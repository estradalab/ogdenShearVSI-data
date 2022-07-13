classdef DataExplorer < handle
    
    % DataExplorer
    % Version 2.0.0
    %
    % DataExplorer is the class file for the Volumetric Data Explorer app.
    % This app provides an interactive environment to explore higher
    % dimensional data using some of MATLAB's abilities for volumetric
    % visualization and animation. It was designed for data that was
    % measured in a 3D grid of data points, for example temperature or wind
    % speed taken at each point in a 3D space. If this data is taken over
    % time, it can also be animated. Any data that fits the form v =
    % f(x,y,z) or v = f(x,y,z,t) can be used.
    %
    % DataExplorer can be called from the command line to bring up the app:
    %
    % >> app = DataExplorer
    %
    % It can also be used in a functional form for inputting the data to be
    % visualized, v, the dimensional values, x y and z, and the time
    % vector, t. It can be used with the following syntax:
    %
    % >> app = DataExplorer(v)            % Specify data only
    % >> app = DataExplorer(x,y,z,v)      % Specify dimensions and data
    % >> app = DataExplorer(x,y,z,v,t)    % Specify dimensions, data and time
    % >> app = DataExplorer([],[],[],v,t) % Specify data and time, but not
    % dimensions
    %
    % Authored by Adam Filion
    % Copyright 2013-2016, The MathWorks, Inc.
    
    properties
        Figure % handle to figure
        % menu properties
        Data % holds the data to be visualized
        ImportMenu % holds handles to import data menu
        OptionsMenu % holds handles to Options menu
        HelpMenu % holds handles to Help menu
        DimValues % 1x3 cell array, holds X, Y, Z data if used
        TimeValues % holds time values if used
        FontSize % size of edit box fonts
        % Top level panel properties
        SlicePanel % handle to panel containing slice visualization
        ISOPanel % handle to panel containing isosurface visualization
        ControlPanel % handle to panel containing controls
        VBox % vertial box, HBox on top, ControlPanel on bottom
        HBox % horizontal box, contains SlicePanel on left and ISOPanel on right
        % Control panel content properties
        HBoxCon % horizontal box inside ControlPanel, from left to right contains DelayPanel, PlayButton, DisplayPanel, and SliderPanel
        DelayPanel % panel containing VBoxDelay
        DisplayPanel % panel containing VBoxDisplay
        SliderPanel % panel containing HBoxSlider
        HBoxSlider % horizontal box containing PlaySlider
        VBoxDelay % veritcal box containing DelayPlus, DelayText and DelayMinus
        VBoxDisplay % vertical box containing HBoxStartAt
        PlayVBox % vertical box containing PlayButton, AnimateButton, RecordButton
        PlaySlider % handle to slider control
        DelayPlus % handle to delay plus button
        DelayText % handle to edit box showing delay amount
        Delay % holds the delay value
        DelayMinus % handle to delay minus button
        HBoxStartAt % array of handles containing SampleText, SampleNum, TimeText and TimeNum
        SampleText % handle to text box for samples
        SampleNum % handle to edit box for starting sample number
        CurrentSample % holds the current sample being displayed
        TimeText % handle to text box for time
        TimeNum % handle to edit box for starting time number
        CurrentTime % holds the current time being displayed
        PlayButton % handle to PLAY button
        AnimateButton % handle to ANIMATE button
        RecordButton % handle to RECORD button
        % Slice panel content properties
        SliceGrid % handle to Grid used in SlicePanel, contains GridComp
        GridComp % 3x3 cell array of handles to different components in SliceGrid
        SliceAxis % handle to axis used for slice
        SliceValuePanel % 1x3 array of handles to panels used in GridComp{3,2}
        SliceValueVBox % 1x3 array of handles to verital boxes used in SliceValuePanel
        SliceValue % 1x3 array of handles to edit boxes used in SliceValueVBox
        SliceColorbar % handle to colorbar, if created
        SliceHandle % 3x1 array of handles to slices
        SliceSliderListener % 3x1 array of handles to listeners for slice callbacks
        % ISO panel content properties
        HBoxISO % horizontal box containing ISOPanel components
        ISOAxis % handle to axis for isosurface plot
        VBoxISO % handle to vertical box containing ISOPanel controls
        ISOLevelPanel % 1x2 array of handles to panels for ISO Levels
        VBoxISOLevel % 1x2 array of handles to vertical boxes used in ISOLevelPanel
        ISOLevelText % 1x2 array of handles to edit boxes used in VBoxISOLevel
        ISOAlphaPanel % 1x2 array of handles to panels for ISO Alphas
        VBoxISOAlpha % 1x2 array of handles to vertical boxes used in ISOLevelAlpha
        ISOAlphaText % 1x2 array of handles to edit boxes used in VBoxISOAlpha
        ISOGridEmpty % 1x4 array of handles to empty holding spots in VBoxISO
        PatchHandle % 1x2 array of handles to patch plots
        HBoxISOSlider % Horizontal box containing sliders for ISO panel
        ISOSlider % 1x2 array of handles to sliders in ISO panel
        ISOContainer % handle to undocumented container for ISO plot and colorbar
        ISOColorbar % handle to colorbar for ISO plot
        % Recording properties
        RecordFileName % name of file to record to
        RecordFileFormat % name of file format for recording, .avi, .mp4 or .mj2
        RecordFrameRate % frame rate of recording
        RecordMode % name of recording mode, either 'Specified' or 'Manual'
        RecordBy % name of record by selection, either 'Samples' or 'Slice/ISO'
        RecordArea % area of app to record
        RepeatOn % is repeat on?
        RepeatTimes % number of times to repeat the recording
        RecordSliceValuesExpr % expression for slice values
        RecordSliceValues % slice values to record
        RecordSliceValuesExprTemp % temporary expression for slice values
        RecordSliceValuesTemp % temporary slice values to record
        RecordISOValues % iso values to record
        RecordISOValuesExpr % expression to evaluate to generate iso values
        RecordISOValuesTemp % temporary iso values to record
        RecordISOValuesExprTemp % temporary expression to evaluate to generate iso values
        RecordSampleValues % sample values to record
        RecordSampleValuesExpr % expression to evaluate to generate sample values
        RecordSampleValuesTemp % sample values to record
        RecordSampleValuesExprTemp % expression to evaluate to generate sample values
        WriteObject % handle to avi file
        % Recording Window components
        RecordOptionFigure % handle to figure used for setting recording options
        RecordGrid % handle to grid for arranging panels in recording window UI
        RecordHBox % handle to horizontal box in Animation & Recording Options UI
        RecordVBox % 1x2 array of handles to vertical boxes used in RecordHBox
        RecordPanel % 3x2 array of handles to panels used in RecordVBox
        RecordButtonGroup % 1x3 array of handles to radio button groups in 'Recording Mode', 'Area to Record' and 'Animate By', in that order
        RecordStartRadio % 1x2 array of toggle buttons in 'Recording Mode' section
        StartStopHBox
        StartStopVBox
        StartStopTextBox
        RecordAreaRadio
        AnimateVBox
        AnimateHBox
        AnimateByRadio
        AnimateText
        RepeatOptions
        SamplesText
        SamplesValues
        SamplesVBox
        SampleSliceGrid
        SliceCheckbox
        SliceCheckboxVal
        SliceEdit
        SampleISOGrid
        ISOCheckbox
        ISOCheckboxVal
        ISOEdit
        FileOptionsVBox
        FileOptionsHBox
        FileOptionsText
        FileOptionsEdit
        FileOptionsButtonGroup
        FileOptionsRadio
        CloseWindowVBox
        CloseWindowHBox
        CloseWindowButton
        RecordEmpty
    end
    
    
    methods
        % Constructor
        function app = DataExplorer(varargin)
            % DataExplorer() is the constructor the DataExplorer class. It
            % creates the main app window, and can also take inputs to
            % programmatically accept data to visualize. See the app help
            % documentation for more info on programatic interface.
            
            % Need to check if GUI Layout Toolbox is available as this is
            % used in constructing the UI. First try for newer versions
            % that use new toolbox packaging and installation.
            installInfo = ver;
            containsTbx = arrayfun(@(x) strcmp(x.Name,'GUI Layout Toolbox'),installInfo,'UniformOutput',false);
            containsTbx = any(cell2mat(containsTbx));
            if ~containsTbx && ~verLessThan('matlab','8.4')
                error(['This app requires that the GUI Layout Toolbox first be installed from the MATLAB File Exchange: '...
                    'http://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox'])
            end
            % For older versions where toolbox files are simply put
            % somewhere on the path, check if one of the toolbox functions
            % is visible and use this as a proxy for if the toolbox is
            % available (not foolproof!)
            fcnFound = which('uiextras.HBox','-all');
            if isempty(fcnFound) && verLessThan('matlab','8.4')
                error(['It appears you may not have the GUI Layout Toolbox available and on the path. '...
                    'You can get it from the MATLAB File Exchange: '...
                    'http://www.mathworks.com/matlabcentral/fileexchange/27758-gui-layout-toolbox'])
            end
            % If multiple copies of GUI Layout Toolbox function is found,
            % warn that this may cause problems.
            if length(fcnFound) > 1
                warning(['It appears you may have multiple versions of the GUI Layout Toolbox available, '...
                    'or have shadowed it with custom functions. This may cause strange or erroneous behavior. '...
                    'Ensure that the proper versions are on top of your MATLAB path.'])
            end
            
            % create figure
            screensize = get(0,'ScreenSize');
            app.Figure = figure('Position',[(screensize(3)-1100)/2 (screensize(4)-600)/2-20 ...
                1100 600],'NumberTitle','off','Name','Volumetric Data Explorer','ResizeFcn',@app.ResizeFcnCB,...
                'HandleVisibility','off','CloseRequestFcn',@app.CloseRequestFcnCB,'Renderer','OpenGL');
            
            % find toolbar toggle tools and delete unwanted ones
            delete(findall(app.Figure,'Tag','Exploration.Brushing'));
            tt = findall(app.Figure,'Type','uitoggletool');delete(tt([1:3 6 9]));
            % find and delete toolbar push tools
            pt = findall(app.Figure,'Type','uipushtool');delete(pt([1:3,5:end]));
            % find and delete menus
            delete(findall(app.Figure,'Type','uimenu'));
            % custom data tip display
            dcm_obj = datacursormode(app.Figure);
            set(dcm_obj,'UpdateFcn',@app.DataTipCB)
            % create custom menu
            app.FontSize = 18;
            app.ImportMenu(1) = uimenu(app.Figure,'Label','Import Data');
            app.ImportMenu(2) = uimenu(app.ImportMenu(1),'Label','Import from Workspace','Callback',@app.ImportWorkspaceCB);
            app.ImportMenu(3) = uimenu(app.ImportMenu(1),'Label','Example Data');
            app.ImportMenu(4) = uimenu(app.ImportMenu(3),'Label','Oscillating Ellipsoid','Callback',@app.LoadEllipsoidCB);
            app.ImportMenu(5) = uimenu(app.ImportMenu(3),'Label','Fluid Flow','Callback',@app.LoadFluidFlowCB);
            app.OptionsMenu(1) = uimenu(app.Figure,'Label','Options');
            app.OptionsMenu(2) = uimenu(app.OptionsMenu(1),'Label','Save/Load options');
            app.OptionsMenu(3) = uimenu(app.OptionsMenu(2),'Label','Save current configuration','Callback',@app.SaveConfCB);
            app.OptionsMenu(4) = uimenu(app.OptionsMenu(2),'Label','Load saved configuration','Callback',@app.LoadConfCB);
            app.OptionsMenu(5) = uimenu(app.OptionsMenu(2),'Label','Delete saved configuration','Callback',@app.DeleteConfCB);
            app.OptionsMenu(6) = uimenu(app.OptionsMenu(2),'Label','Set current as default','Callback',@app.DefaultConfCB);
            app.OptionsMenu(7) = uimenu(app.OptionsMenu(2),'Label','Restore factory defaults','Callback',@app.RestoreConfCB);
            app.OptionsMenu(8) = uimenu(app.OptionsMenu(1),'Label','Link Rotation','Checked','on','Callback',@app.LinkRotationCB);
            app.OptionsMenu(9) = uimenu(app.OptionsMenu(1),'Label','Equalize Axes','Callback',@app.EqualAxesCB);
            app.OptionsMenu(10) = uimenu(app.OptionsMenu(1),'Label','Loop PLAY Button','Checked','off','Callback',@app.LoopAnimationCB);
            app.OptionsMenu(11) = uimenu(app.OptionsMenu(1),'Label','Animation & Recording Options','Checked','off','Callback',@app.RecordAnimationCB);
            app.OptionsMenu(12) = uimenu(app.OptionsMenu(1),'Label','Add Colorbars','Checked','off','Callback',@app.ColorbarCB);
            app.OptionsMenu(13) = uimenu(app.OptionsMenu(1),'Label','Set Colormap','Checked','off','Callback',@app.SetColormapCB);
            app.OptionsMenu(14) = uimenu(app.OptionsMenu(1),'Label','Set Font Size','Checked','off','Callback',@app.SetFontSizeCB);
            app.OptionsMenu(15) = uimenu(app.OptionsMenu(1),'Label','Set Axes Names','Checked','off','Callback',@app.SetAxesNamesCB);
            app.OptionsMenu(16) = uimenu(app.OptionsMenu(1),'Label','Beep on Warning','Checked','on','Callback',@app.BeepCB);
            app.OptionsMenu(17) = uimenu(app.OptionsMenu(1),'Label','Dynamic ISO Color','Checked','off','Callback',@app.DynamicISOCB);
            app.HelpMenu(1) = uimenu(app.Figure,'Label','Help');
            app.HelpMenu(2) = uimenu(app.HelpMenu(1),'Label','Open Help Documentation','Checked','off','Callback',@app.HelpCB);
            
            % create top level panels
            app.VBox = uiextras.VBox('Parent',app.Figure,'Padding',2,'Spacing',2);
            app.HBox = uiextras.HBox('Parent',app.VBox,'Padding',2,'Spacing',2);
            app.SlicePanel = uiextras.BoxPanel('Parent',app.HBox,'Title','Slices','FontSize',12);
            app.ISOPanel = uiextras.BoxPanel('Parent',app.HBox,'Title','Isosurfaces','FontSize',12);
            app.ControlPanel = uiextras.BoxPanel('Parent',app.VBox,'Title','Controls','FontSize',12);
            set(app.VBox, 'Sizes', [-1 150]);
            
            % create control panel contents
            app.HBoxCon = uiextras.HBox('Parent',app.ControlPanel,'Spacing',5,'Padding',2);
            app.DelayPanel = uiextras.BoxPanel('Parent',app.HBoxCon,'Title','Delay Animation (sec)');
            app.PlayVBox = uiextras.VBox('Parent',app.HBoxCon,'Spacing',2,'Padding',2);
            app.PlayButton = uicontrol('Parent',app.PlayVBox,'Style','togglebutton','String','PLAY','Value',0,'FontSize',24,'Callback',@app.PlayButtonCB);
            app.AnimateButton = uicontrol('Parent',app.PlayVBox,'Style','togglebutton','String','ANIMATE','Value',0,'FontSize',24,'Callback',@app.AnimateButtonCB);
            app.RecordButton = uicontrol('Parent',app.PlayVBox,'Style','togglebutton','String','RECORD','Value',0,'FontSize',24,'Callback',@app.RecordButtonCB);
            set(app.PlayVBox,'Sizes',[-1 -1 -1]);
            app.DisplayPanel = uiextras.BoxPanel('Parent',app.HBoxCon,'Title','Displaying:');
            app.PlaySlider = uicontrol('Parent',app.HBoxCon,'Style', 'slider','Min',1,'Max',50,'Value',1,'Callback', @app.PlaySliderCB);
            addlistener(app.PlaySlider,'ContinuousValueChange',@app.PlaySliderCB);
            app.CurrentSample = get(app.PlaySlider,'Value');
            set(app.HBoxCon, 'Sizes', [120 150 200 -1]);
            app.VBoxDelay = uiextras.VBox('Parent',app.DelayPanel,'Spacing',2,'Padding',2);
            app.VBoxDisplay = uiextras.VBox('Parent',app.DisplayPanel,'Spacing',2,'Padding',2);
            app.DelayPlus = uicontrol('Parent',app.VBoxDelay,'String','+','FontSize',18,'Callback',@app.DelayPlusCB);
            app.DelayText = uicontrol('Parent',app.VBoxDelay,'Style','edit','String','0','FontSize',app.FontSize,'Callback',@app.DelayTextCB,'BackgroundColor',[1 1 1]);
            app.Delay = str2double(get(app.DelayText,'String'));
            app.DelayMinus = uicontrol('Parent',app.VBoxDelay,'String','-','FontSize',18,'Callback',@app.DelayMinusCB);
            app.HBoxStartAt(1) = uiextras.HBox('Parent',app.VBoxDisplay);
            app.SampleText = uicontrol('Parent',app.HBoxStartAt(1),'Style','text','String','Current sample:','FontSize',16);
            app.SampleNum = uicontrol('Parent',app.HBoxStartAt(1),'Style','edit','String','1','FontSize',app.FontSize,'Callback',@app.SampleNumCB,'BackgroundColor',[1 1 1]);
            app.HBoxStartAt(2) = uiextras.HBox('Parent',app.VBoxDisplay);
            app.TimeText = uicontrol('Parent',app.HBoxStartAt(2),'Style','text','String','Current time:','FontSize',16);
            app.TimeNum = uicontrol('Parent',app.HBoxStartAt(2),'Style','edit','String','N/A','FontSize',app.FontSize,'Callback',@app.TimeNumCB,'BackgroundColor',[1 1 1]);
            
            % create slice panel contents
            app.SliceGrid = uiextras.Grid('Parent',app.SlicePanel,'Spacing',10,'Padding',10);
            app.GridComp = cell(3,2);
            app.GridComp{1,1} = uiextras.Empty('Parent',app.SliceGrid);
            app.GridComp{2,1} = uicontrol('Parent',app.SliceGrid,'Style','slider');
            app.GridComp{3,1} = uiextras.Empty('Parent',app.SliceGrid);
            app.GridComp{1,2} = uicontrol('Parent',app.SliceGrid,'Style','slider');
            app.GridComp{2,2} = uicontainer('Parent',app.SliceGrid);
            app.SliceAxis = axes('Parent',app.GridComp{2,2});grid(app.SliceAxis,'on');
            xlabel(app.SliceAxis,'X');ylabel(app.SliceAxis,'Y');zlabel(app.SliceAxis,'Z');
            app.GridComp{3,2} = uicontrol('Parent',app.SliceGrid,'Style','slider');
            app.GridComp{1,3} = uiextras.Empty('Parent',app.SliceGrid);
            app.GridComp{2,3} = uiextras.VBox('Parent',app.SliceGrid,'Spacing',2,'Padding',2);
            app.GridComp{3,3} = uiextras.Empty('Parent',app.SliceGrid);
            addlistener(app.GridComp{2,1},'ContinuousValueChange',@app.RedrawSlice);
            addlistener(app.GridComp{1,2},'ContinuousValueChange',@app.RedrawSlice);
            addlistener(app.GridComp{3,2},'ContinuousValueChange',@app.RedrawSlice);
            app.SliceValuePanel(1) = uiextras.BoxPanel('Parent',app.GridComp{2,3},'Title','X Slice At:');
            app.SliceValueVBox(1) = uiextras.VBox('Parent',app.SliceValuePanel(1),'Spacing',2,'Padding',2);
            app.SliceValue(1) = uicontrol('Parent',app.SliceValueVBox(1),'Style','edit','String','0','FontSize',app.FontSize,'Callback',@app.SliceValueCB,'BackgroundColor',[1 1 1]);
            app.SliceValuePanel(2) = uiextras.BoxPanel('Parent',app.GridComp{2,3},'Title','Y Slice At:');
            app.SliceValueVBox(2) = uiextras.VBox('Parent',app.SliceValuePanel(2),'Spacing',2,'Padding',2);
            app.SliceValue(2) = uicontrol('Parent',app.SliceValueVBox(2),'Style','edit','String','0','FontSize',app.FontSize,'Callback',@app.SliceValueCB,'BackgroundColor',[1 1 1]);
            app.SliceValuePanel(3) = uiextras.BoxPanel('Parent',app.GridComp{2,3},'Title','Z Slice At:');
            app.SliceValueVBox(3) = uiextras.VBox('Parent',app.SliceValuePanel(3),'Spacing',2,'Padding',2);
            app.SliceValue(3) = uicontrol('Parent',app.SliceValueVBox(3),'Style','edit','String','0','FontSize',app.FontSize,'Callback',@app.SliceValueCB,'BackgroundColor',[1 1 1]);
            set(app.SliceGrid, 'ColumnSizes', [30 -1 80], 'RowSizes', [30 -1 30] );
            
            % create isosurface panel contents
            app.HBoxISO = uiextras.HBox('Parent',app.ISOPanel,'Spacing',5,'Padding',2);
            app.ISOContainer = uicontainer('Parent',app.HBoxISO);
            app.ISOAxis = axes('Parent',app.ISOContainer);grid(app.ISOAxis,'on');
            xlabel(app.ISOAxis,'X');ylabel(app.ISOAxis,'Y');zlabel(app.ISOAxis,'Z');
            app.HBoxISOSlider = uiextras.HBox('Parent',app.HBoxISO,'Spacing',5,'Padding',2);
            app.VBoxISO = uiextras.VBox('Parent',app.HBoxISO,'Spacing',5,'Padding',2);
            app.ISOGridEmpty(1) = uiextras.Empty('Parent',app.VBoxISO);
            app.ISOLevelPanel(1) = uiextras.BoxPanel('Parent',app.VBoxISO,'Title','ISO Level 1');
            app.VBoxISOLevel(1) = uiextras.VBox('Parent',app.ISOLevelPanel(1),'Spacing',2,'Padding',2);
            app.ISOLevelText(1) = uicontrol('Parent',app.VBoxISOLevel(1),'Style','edit','String','2','FontSize',app.FontSize,'Callback',@app.ISOTextCB,'BackgroundColor',[1 1 1]);
            app.ISOAlphaPanel(1) = uiextras.BoxPanel('Parent',app.VBoxISO,'Title','ISO Alpha 1');
            app.VBoxISOAlpha(1) = uiextras.VBox('Parent',app.ISOAlphaPanel(1),'Spacing',2,'Padding',2);
            app.ISOAlphaText(1) = uicontrol('Parent',app.VBoxISOAlpha(1),'Style','edit','String','0.5','FontSize',app.FontSize,'Callback',@app.ISOTextCB,'BackgroundColor',[1 1 1]);
            app.ISOGridEmpty(2) = uiextras.Empty('Parent',app.VBoxISO);
            app.ISOGridEmpty(3) = uiextras.Empty('Parent',app.VBoxISO);
            app.ISOLevelPanel(2) = uiextras.BoxPanel('Parent',app.VBoxISO,'Title','ISO Level 2');
            app.VBoxISOLevel(2) = uiextras.VBox('Parent',app.ISOLevelPanel(2),'Spacing',2,'Padding',2);
            app.ISOLevelText(2) = uicontrol('Parent',app.VBoxISOLevel(2),'Style','edit','String','1','FontSize',app.FontSize,'Callback',@app.ISOTextCB,'BackgroundColor',[1 1 1]);
            app.ISOAlphaPanel(2) = uiextras.BoxPanel('Parent',app.VBoxISO,'Title','ISO Alpha 2');
            app.VBoxISOAlpha(2) = uiextras.VBox('Parent',app.ISOAlphaPanel(2),'Spacing',2,'Padding',2);
            app.ISOAlphaText(2) = uicontrol('Parent',app.VBoxISOAlpha(2),'Style','edit','String','0.7','FontSize',app.FontSize,'Callback',@app.ISOTextCB,'BackgroundColor',[1 1 1]);
            app.ISOGridEmpty(4) = uiextras.Empty('Parent',app.VBoxISO);
            app.ISOSlider(1) = uicontrol('Parent',app.HBoxISOSlider,'Style','slider');
            app.ISOSlider(2) = uicontrol('Parent',app.HBoxISOSlider,'Style','slider');
            addlistener(app.ISOSlider(1),'ContinuousValueChange',@app.RedrawISO);
            addlistener(app.ISOSlider(2),'ContinuousValueChange',@app.RedrawISO);
            set(app.HBoxISO,'Sizes',[-1 60 80]);
            
            % update font sizes based on version
            if ~verLessThan('matlab','8.4') % font sizes changed in R2014b, version 8.4
                set([app.PlayButton,app.AnimateButton,app.RecordButton],'FontSize',18)
                set([app.SampleText,app.TimeText],'FontSize',11)
                set([app.ISOLevelText,app.ISOAlphaText],'FontSize',10)
            end
            
            % disable the GUI components until data is imported, then link
            % axes rotation and set the viewing angle
            EnableDisableFig(app,'disable');
            
            % link axes rotation
            setappdata(app.SliceAxis,'graphics_linkprop',linkprop([app.SliceAxis app.ISOAxis],{'CameraPosition','CameraUpVector'}));
            view(app.SliceAxis,-30,50);
            
            % Factory defaults for Animation & Recording Options properties
            app.RecordFileName = 'animation';
            app.RecordFileFormat = '.avi';
            app.RecordFrameRate = 10;
            app.RecordMode = 'Manual';
            app.RecordArea = 'Whole App';
            app.RecordBy = 'Slice/ISO';
            app.RepeatOn = 0;
            app.RepeatTimes = 0;
            app.RecordSliceValuesExpr = {'linspace(minx,maxx,20)','maxy','minz'};
            app.RecordSliceValues = {1,1,1};
            app.RecordISOValuesExpr = {'linspace(mind,maxd,20)','mind+(maxd-mind)/2'};
            app.RecordISOValues = {1,1};
            app.RecordSampleValuesExpr = '1:slength';
            app.RecordSampleValues = 1;
            app.SliceCheckboxVal = [1,1,1];
            app.ISOCheckboxVal = [1,1];
            
            % turn off configuration options until data is loaded
            set(app.OptionsMenu(2),'Enable','off')
            set(app.OptionsMenu(9),'Enable','off')
                        
            % Save factory and default settings if file does not already
            % exist
            if exist('DE_confs.mat','file')==2
                load DE_confs
                if ~exist('conf_factory','var') || ~exist('conf_default','var')
                    conf_factory = getSettings(app);
                    conf_default = conf_factory; %#ok
                    save('DE_confs.mat','conf_factory','conf_default','-append');
                end
            elseif ~(exist('DE_confs.mat','file')==2)
                conf_factory = getSettings(app);
                conf_default = conf_factory; %#ok
                save('DE_confs.mat','conf_factory','conf_default');
            end
            
            % if functional syntax form is used, save into data variable in
            % same form as returned by uigetvariable so validation code can
            % be reused
            if numel(varargin) > 0
                names = {'','X','Y','Z'}; % put in same format returned by uigetvariables
            end
            if numel(varargin) == 1
                data{1} = varargin{1};
                data{2} = [];data{3} = [];data{4} = [];data{5} = [];
                ValidateData(app,data,names);
            elseif numel(varargin) == 4
                data{5} = [];data{1} = varargin{4};
                data{2} = varargin{1};data{3} = varargin{2};data{4} = varargin{3};
                ValidateData(app,data,names);
            elseif numel(varargin) == 5
                data{5} = varargin{5};data{1} = varargin{4};
                data{2} = varargin{1};data{3} = varargin{2};data{4} = varargin{3};
                ValidateData(app,data,names);
            elseif numel(varargin) > 0
                warning('Incorrect number of inputs');BeepOnWarn(app);
            end
        end
        
        function ResizeFcnCB(app,~,~)
            % ResizeFcnCB(app) is run after the DataExplorer figure window is
            % resized. This is needed to make sure the colorbar does not
            % overlap with the plot after resizing.
            if ishandle(app.SliceColorbar) % if colorbar exists
                delete(app.SliceColorbar) % delete it
                drawnow; % force graphics update, doesn't work without this!
                app.SliceColorbar = colorbar('peer',app.SliceAxis); % create new colorbar
            end
            if ishandle(app.ISOColorbar) % if colorbar exists
                delete(app.ISOColorbar) % delete it
                drawnow; % force graphics update, doesn't work without this!
                app.ISOColorbar = colorbar('peer',app.ISOAxis); % create new colorbar
            end
        end
        
        function CloseRequestFcnCB(app,~,~)
            % CloseRequestFcnCB(app) is run when the figure is closed. It closes
            % the recording file if it is actively being recorded to at the
            % time of closing.
            if get(app.RecordButton,'Value') == 1
                close(app.WriteObject);
            end
            delete(app.Figure);
        end
        
        % menu functions
        function SaveConfCB(app,~,~)
            % SaveConfCB(app) is the callback for the Options -> Save/Load
            % options -> Save current configuration menu option. It
            % will save the current set of options to a file in the app's
            % directory for later use.
            prompt = 'Enter a name for this set of preferences.';
            def = {'mypref'};
            answer = inputdlg(prompt,'Save Preferences',1,def);
            if ~(strcmp(answer{1},'default') || strcmp(answer{1},'factory'))
                try
                    eval(['conf_' answer{1} ' = getSettings(app);'])
                    save('DE_confs.mat',['conf_' answer{1}],'-append')
                catch
                    error('Preference name must be valid MATLAB variable name')
                end
            else
                warning('''default'' and ''factory'' are reserved preferences names, please resave under a different name');
                BeepOnWarn(app);
            end
        end
        
        function LoadConfCB(app,~,~)
            % LoadConfCB(app) is the callback for the Options -> Save/Load
            % options -> Load saved configuration menu option. It
            % will load a set of saved option settings from a mat file in
            % the app's directory.
            
            confs = load('DE_confs');
            names = fieldnames(confs);
            for ii = length(names):-1:1
                if ~isstruct(confs.(names{ii}))
                    names(ii) = [];
                end
            end
            names = cellfun(@(x) x(6:end),names,'UniformOutput',false);
            [answer,ok] = listdlg('ListString',names,'SelectionMode','single','Name','Select Preferences File','PromptString','Select a preferences file to load');
            if ok && ~strcmp(names{answer},'') && ~isempty(names{answer})
                setSettings(app,confs.(['conf_',names{answer}]));
                % need to redraw to use updated settings
                RedrawSlice(app)
                RedrawISO(app)
            end
        end
        
        function DeleteConfCB(~,~,~)
            % DeleteConfCB(app) is the callback for the Options ->
            % Save/Load options -> Deleted saved preferences menu option.
            % It allows the user to select a custom saved preferences file
            % and delete it.
            confs = load('DE_confs');
            names = fieldnames(confs);
            for ii = length(names):-1:1
                if ~isstruct(confs.(names{ii}))
                    names(ii) = [];
                end
            end
            idxrem = cellfun(@(x)(strcmp(x,'conf_factory')|strcmp(x,'conf_default')),names,'UniformOutput',false);
            idxrem = cell2mat(idxrem);
            names(idxrem) = [];
            names = cellfun(@(x) x(6:end),names,'UniformOutput',false);
            [answer,ok] = listdlg('ListString',names,'SelectionMode','single','Name','Select Preferences File','PromptString','Select a preferences file to delete');
            if ok && ~strcmp(names{answer},'') && ~isempty(names{answer})
                confs.(['conf_' names{answer}]) = [];
                save('DE_confs','-struct','confs')
            end
        end
        
        function DefaultConfCB(app,~,~)
            % DefaultConfCB(app) is the callback for the Options ->
            % Save/Load options -> Set Current As Default menu option. This
            % will set the current options configuration as both the
            % current and default set to use.
            conf = getSettings(app);
            confs = load('DE_confs');
            confs.conf_default = conf;
            save('DE_confs','-struct','confs');
        end
        
        function RestoreConfCB(app,~,~)
            % RestoreConfCB(app) is the callback for the Options ->
            % Save/Load options -> Restore factory defaults menu option.
            % This will reset both the current and default options
            % configuration to their original settings.
            confs = load('DE_confs');
            confs.conf_default = confs.conf_factory;
            save('DE_confs','-struct','confs')
            setSettings(app,confs.conf_factory)
            % need to redraw to use updated settings
            RedrawSlice(app)
            RedrawISO(app)
        end
        
        function conf = getSettings(app)
            % conf = getSettings(app) is a helper function that returns a
            % structure containing all the options information needed for
            % saving/loading options
            conf.LinkRotation = get(app.OptionsMenu(8),'Checked');
            conf.EqualAxes = get(app.OptionsMenu(9),'Checked');
            conf.LoopAnimation = get(app.OptionsMenu(10),'Checked');
            conf.RecordAnimation = get(app.OptionsMenu(11),'Checked');
            conf.Colorbar = get(app.OptionsMenu(12),'Checked');
            conf.Colormap = get(app.Figure,'Colormap');
            conf.FontSize = app.FontSize;
            conf.Beep = get(app.OptionsMenu(16),'Checked');
            conf.DynamicISO = get(app.OptionsMenu(17),'Checked');
            conf.RecordFileName = app.RecordFileName;
            conf.RecordFileFormat = app.RecordFileFormat;
            conf.RecordFrameRate = app.RecordFrameRate;
            conf.RecordMode = app.RecordMode;
            conf.RecordArea = app.RecordArea;
            conf.RecordBy = app.RecordBy;
            conf.RepeatOn = app.RepeatOn;
            conf.RepeatTimes = app.RepeatTimes;
            conf.RecordSliceValuesExpr = app.RecordSliceValuesExpr;
            conf.RecordISOValuesExpr = app.RecordISOValuesExpr;
            conf.RecordSampleValuesExpr = app.RecordSampleValuesExpr;
            conf.RecordSampleValues = app.RecordSampleValues;
            conf.SliceCheckboxVal = app.SliceCheckboxVal;
            conf.ISOCheckboxVal = app.ISOCheckboxVal;
        end
        
        function setSettings(app,conf)
            % setSettings(app,conf) is a helper function that accepts a
            % sturcture created by getSettings and will set all the
            % appropriate options to match the saved options configuration
            if strcmp(conf.LinkRotation,'on')
                % if the saved setting was 'on', set the current option
                % to 'off' and run callback so it gets turned on and
                % link axes
                set(app.OptionsMenu(8),'Checked','off');LinkRotationCB(app);
            else
                set(app.OptionsMenu(8),'Checked','on');LinkRotationCB(app);
            end
            if strcmp(conf.EqualAxes,'on')
                set(app.OptionsMenu(9),'Checked','off');EqualAxesCB(app);
            else
                set(app.OptionsMenu(9),'Checked','on');EqualAxesCB(app);
            end
            % no additional code to run for the LoopAnimation option, so we
            % can set it directly
            set(app.OptionsMenu(10),'Checked',conf.LoopAnimation);
            if strcmp(conf.Colorbar,'on')
                set(app.OptionsMenu(12),'Checked','off');ColorbarCB(app);
            else
                set(app.OptionsMenu(12),'Checked','on');ColorbarCB(app);
            end
            set(app.Figure,'Colormap',conf.Colormap);
            app.FontSize = conf.FontSize;
            set([app.SampleNum, app.TimeNum, app.DelayText, app.SliceValue,...
                app.ISOLevelText, app.ISOAlphaText],'FontSize',app.FontSize);
            set(app.OptionsMenu(16),'Checked',conf.Beep);
            if strcmp(conf.DynamicISO,'on')
                set(app.OptionsMenu(17),'Checked','off');DynamicISOCB(app);
            else
                set(app.OptionsMenu(17),'Checked','on');DynamicISOCB(app);
            end
            app.RecordFileName = conf.RecordFileName;
            app.RecordFileFormat = conf.RecordFileFormat;
            app.RecordFrameRate = conf.RecordFrameRate;
            app.RecordMode = conf.RecordMode;
            app.RecordArea = conf.RecordArea;
            app.RecordBy = conf.RecordBy;
            app.RepeatOn = conf.RepeatOn;
            app.RepeatTimes = conf.RepeatTimes;
            app.RecordSliceValuesExpr = conf.RecordSliceValuesExpr;
            app.RecordISOValuesExpr = conf.RecordISOValuesExpr;
            app.RecordSampleValuesExpr = conf.RecordSampleValuesExpr;
            app.RecordSampleValues = conf.RecordSampleValues;
            app.SliceCheckboxVal = conf.SliceCheckboxVal;
            app.ISOCheckboxVal = conf.ISOCheckboxVal;
        end
        
        function LoopAnimationCB(app,~,~)
            % LoopAnimationCB(app) is the callback for the Options -> Loop
            % Animation menu option. This toggles whether the animation
            % will loop upon finishing.
            if strcmp(get(app.OptionsMenu(10),'Checked'),'off')
                set(app.OptionsMenu(10),'Checked','on')
            else
                set(app.OptionsMenu(10),'Checked','off')
            end
        end
        
        function DynamicISOCB(app,~,~)
            % DynamicISOCB(app) is the callback  for the Options -> Dynamic
            % ISO Color menu option. Turning this on will cause the
            % isosurface colors to match the colors for the corresponding
            % value in the slice plot using the figures colormap. Turning
            % it off will set them to default values.
            if strcmp(get(app.OptionsMenu(17),'Checked'),'off')
                set(app.OptionsMenu(17),'Checked','on')
                if strcmp(get(app.OptionsMenu(12),'Checked'),'on')
                    % if colorbars are already on, create one for iso
                    app.ISOColorbar = colorbar('peer',app.ISOAxis);
                end
                ISOTextCB(app);
            else
                set(app.OptionsMenu(17),'Checked','off')
                if ~isempty(app.ISOColorbar)
                    if strcmp(get(app.OptionsMenu(12),'Checked'),'on') && ishandle(app.ISOColorbar)
                        % if colorbars are on, delete the iso one so as not to
                        % confuse the meaning of the default isosurface colors
                        delete(app.ISOColorbar)
                    end
                end
                ISOTextCB(app);
            end
        end
        
        function RecordAnimationCB(app,~,~)
            % RecordAnimationCB(app) is the callback for the Options ->
            % Animation & Recording Options menu option. It creates the
            % window for setting options for animating slices and
            % isosurfaces and creating recordings.
            EnableDisableFig(app,'disable');
            % create figure for UI
            screensize = get(0,'ScreenSize');
            app.RecordOptionFigure = figure('Position',[(screensize(3)-1100)/2 (screensize(4)-600)/2-20 ...
                800 650],'NumberTitle','off','Name','Animation & Recording Options',...
                'HandleVisibility','off','CloseRequestFcn',@app.RecordCancelCB,'Resize','off');
            delete(findall(app.RecordOptionFigure,'Type','uimenu'));
            delete(findall(app.RecordOptionFigure,'Type','uitoolbar'));
            app.RecordGrid = uiextras.Grid('Parent',app.RecordOptionFigure,'Padding',2,'Spacing',2);
            % Create the panels
            app.RecordPanel(1,1) = uiextras.BoxPanel('Parent',app.RecordGrid,'Title','Animate By'); % Animate By panel
            app.RecordPanel(2,1) = uiextras.BoxPanel('Parent',app.RecordGrid,'Title','Slice Values To Animate'); % Slice Values panel
            app.RecordPanel(3,1) = uiextras.BoxPanel('Parent',app.RecordGrid,'Title','Recording Mode'); % Recording Mode panel
            app.RecordPanel(4,1) = uiextras.BoxPanel('Parent',app.RecordGrid,'Title','File Options'); % File Options panel
            app.RecordPanel(1,2) = uiextras.BoxPanel('Parent',app.RecordGrid,'Title','Samples To Animate'); % Samples Values panel
            app.RecordPanel(2,2) = uiextras.BoxPanel('Parent',app.RecordGrid,'Title','ISO Values To Animate'); % ISO Values panel
            app.RecordPanel(3,2) = uiextras.BoxPanel('Parent',app.RecordGrid,'Title','Area To Record'); % Select Area panel
            app.RecordPanel(4,2) = uiextras.BoxPanel('Parent',app.RecordGrid,'Title','Close Window'); % Close Window panel
            set(app.RecordGrid,'RowSizes',[-1 -1 -1 -1],'ColumnSizes',[-1 -1]);
            % File Options panel
            app.FileOptionsVBox = uiextras.VBox('Parent',app.RecordPanel(4,1),'Padding',2,'Spacing',2);
            app.FileOptionsHBox(1) = uiextras.HBox('Parent',app.FileOptionsVBox,'Padding',2,'Spacing',2);
            app.FileOptionsHBox(2) = uiextras.HBox('Parent',app.FileOptionsVBox,'Padding',2,'Spacing',2);
            app.FileOptionsHBox(3) = uiextras.HBox('Parent',app.FileOptionsVBox,'Padding',2,'Spacing',2);
            app.FileOptionsText(1) = uicontrol('Style','text','Parent',app.FileOptionsHBox(1),'String','File Name:');
            app.FileOptionsEdit(1) = uicontrol('Style','edit','Parent',app.FileOptionsHBox(1),'String',app.RecordFileName,'BackgroundColor',[1 1 1]);%,'Callback',@app.FileNameCB);
            app.FileOptionsText(2) = uicontrol('Style','text','Parent',app.FileOptionsHBox(2),'String','Frame Rate (1/s):');
            app.FileOptionsEdit(2) = uicontrol('Style','edit','Parent',app.FileOptionsHBox(2),'String',num2str(app.RecordFrameRate),'BackgroundColor',[1 1 1],'Callback',@app.FrameRateCB);
            app.FileOptionsText(3) = uicontrol('Style','text','Parent',app.FileOptionsHBox(3),'String','File Format:');
            app.FileOptionsButtonGroup = uibuttongroup('Parent',app.FileOptionsHBox(3));
            app.FileOptionsRadio(1) = uicontrol('Style','radiobutton','String','.avi',...
                'Parent',app.FileOptionsButtonGroup,'pos',[0 2 45 30]);
            app.FileOptionsRadio(2) = uicontrol('Style','radiobutton','String','.mp4',...
                'Parent',app.FileOptionsButtonGroup,'pos',[60 2 55 30]);
            set(app.FileOptionsHBox,'Widths',[-1 -1.5])
            % set the radio option to be equal to the currently saved format
            set(findall(app.FileOptionsRadio,'String',app.RecordFileFormat),'Value',1)
            % Close Window panel
            app.CloseWindowVBox = uiextras.VBox('Parent',app.RecordPanel(4,2),'Padding',2,'Spacing',2);
            app.RecordEmpty(1) = uiextras.Empty('Parent',app.CloseWindowVBox);
            app.CloseWindowHBox = uiextras.HBox('Parent',app.CloseWindowVBox,'Padding',2,'Spacing',2);
            app.RecordEmpty(2) = uiextras.Empty('Parent',app.CloseWindowVBox);
            app.CloseWindowButton(1) = uicontrol('Style','pushbutton','Parent',app.CloseWindowHBox,'String','OK','Callback',@app.RecordOKCB);
            app.CloseWindowButton(2) = uicontrol('Style','pushbutton','Parent',app.CloseWindowHBox,'String','Cancel','Callback',@app.RecordCancelCB);
            app.CloseWindowButton(3) = uicontrol('Style','pushbutton','Parent',app.CloseWindowHBox,'String','Help','Callback',@app.RecordHelpCB);
            % Recording Mode panel
            app.StartStopHBox = uiextras.HBox('Parent',app.RecordPanel(3,1),'Padding',2,'Spacing',2);
            app.RecordButtonGroup(1) = uibuttongroup('Parent',app.StartStopHBox);%,'SelectionChangeFcn',@app.AnimateByCB);
            app.RecordStartRadio(1) = uicontrol('Style','radiobutton','String','Manual',...
                'parent',app.RecordButtonGroup(1),'pos',[30 60 80 30]);
            app.RecordStartRadio(2) = uicontrol('Style','radiobutton','String','Specified',...
                'parent',app.RecordButtonGroup(1),'pos',[30 20 80 30]);
            app.StartStopTextBox = uicontrol('Style','text','Parent',app.StartStopHBox,'String',...
                ['When recording is on, if you have selected ''Specified'', then the app will only record when using the ANIMATE button. ' ...
                'If you have selected ''Manual'', then every change to the plots will trigger the app to record another frame.']);
            set(app.StartStopHBox,'Sizes',[-1,-1.5])
            % set the radio option to be equal to the currently saved recording mode
            set(findall(app.RecordStartRadio,'String',app.RecordMode),'Value',1);
            % Select Area panel
            app.RecordButtonGroup(2) = uibuttongroup('Parent',app.RecordPanel(3,2));
            app.RecordAreaRadio(1) = uicontrol('Style','radiobutton','String','Whole App',...
                'parent',app.RecordButtonGroup(2),'pos',[20 60 120 30]);
            app.RecordAreaRadio(2) = uicontrol('Style','radiobutton','String','Plot Panels Only',...
                'parent',app.RecordButtonGroup(2),'pos',[150 60 120 30]);
            app.RecordAreaRadio(3) = uicontrol('Style','radiobutton','String','Slice Panel Only',...
                'parent',app.RecordButtonGroup(2),'pos',[20 20 120 30]);
            app.RecordAreaRadio(4) = uicontrol('Style','radiobutton','String','ISO Panel Only',...
                'parent',app.RecordButtonGroup(2),'pos',[150 20 120 30]);
            set(findall(app.RecordAreaRadio,'String',app.RecordArea),'Value',1);
            % Animate By panel
            app.AnimateVBox = uiextras.VBox('Parent',app.RecordPanel(1,1));
            app.AnimateHBox(1) = uiextras.HBox('Parent',app.AnimateVBox); 
            app.RecordButtonGroup(3) = uibuttongroup('Parent',app.AnimateHBox(1),'SelectionChangeFcn',@app.AnimateByCB);
            app.AnimateByRadio(1) = uicontrol('Style','radiobutton','String','Samples',...
                'parent',app.RecordButtonGroup(3),'pos',[40 40 80 20]);
            app.AnimateByRadio(2) = uicontrol('Style','radiobutton','String','Slice/ISO',... % modify here
                'parent',app.RecordButtonGroup(3),'pos',[40 10 80 20]);
            app.AnimateText = uicontrol('Style','text','Parent',app.AnimateHBox(1),...
                'String','You can choose to record a series of samples, similar to hitting PLAY button, or to animate the slices & isosurfaces.');
            app.AnimateHBox(2) = uiextras.HBox('Parent',app.AnimateVBox);
            app.RepeatOptions(1) = uicontrol('Style','checkbox','String','Repeat?','Parent',app.AnimateHBox(2),'Value',app.RepeatOn);
            app.RepeatOptions(2) = uicontrol('Style','text','String','Number of times:','Parent',app.AnimateHBox(2));
            app.RepeatOptions(3) = uicontrol('Style','edit','String',num2str(app.RepeatTimes),'Parent',app.AnimateHBox(2),'Callback',@app.RepeatTimesCB);
            set(app.AnimateVBox, 'Sizes', [-2 -1]);
            % set the radio option to be equal to the currently saved recording by option
            set(findall(app.RecordButtonGroup(3),'String',app.RecordBy),'Value',1);
            % Samples To Record panel
            app.SamplesVBox = uiextras.VBox('Parent',app.RecordPanel(1,2));
            app.RecordEmpty(3) = uiextras.Empty('Parent',app.SamplesVBox);
            app.SamplesText = uicontrol('Style','text','Parent',app.SamplesVBox,...
                'String','Sample numbers to use:');
            app.SamplesValues = uicontrol('Style','edit','Parent',app.SamplesVBox,...
                'String',app.RecordSampleValuesExpr,'BackgroundColor',[1 1 1],'Callback',@app.RecordSamplesCB);
            app.RecordSampleValuesExprTemp = app.RecordSampleValuesExpr;
            app.RecordSampleValuesTemp = app.RecordSampleValues;
            app.RecordEmpty(4) = uiextras.Empty('Parent',app.SamplesVBox);
            RecordSamplesCB(app) % run callback to execute expressions and populate values
            % Slice Values To Record panel
            app.SampleSliceGrid = uiextras.Grid('Parent',app.RecordPanel(2,1),'Spacing',2,'Padding',2);
            app.SliceCheckbox(1) = uicontrol('Style','checkbox','Parent',app.SampleSliceGrid,'String','X','Value',app.SliceCheckboxVal(1),'Callback',@app.RecordWhichSliceCB);
            app.SliceCheckbox(2) = uicontrol('Style','checkbox','Parent',app.SampleSliceGrid,'String','Y','Value',app.SliceCheckboxVal(2),'Callback',@app.RecordWhichSliceCB);
            app.SliceCheckbox(3) = uicontrol('Style','checkbox','Parent',app.SampleSliceGrid,'String','Z','Value',app.SliceCheckboxVal(3),'Callback',@app.RecordWhichSliceCB);
            app.SliceEdit(1) = uicontrol('Style','edit','Parent',app.SampleSliceGrid,'String',app.RecordSliceValuesExpr{1},'BackgroundColor',[1 1 1],'Callback',@app.RecordSliceValuesCB);
            app.SliceEdit(2) = uicontrol('Style','edit','Parent',app.SampleSliceGrid,'String',app.RecordSliceValuesExpr{2},'BackgroundColor',[1 1 1],'Callback',@app.RecordSliceValuesCB);
            app.SliceEdit(3) = uicontrol('Style','edit','Parent',app.SampleSliceGrid,'String',app.RecordSliceValuesExpr{3},'BackgroundColor',[1 1 1],'Callback',@app.RecordSliceValuesCB);
            app.RecordSliceValuesExprTemp = app.RecordSliceValuesExpr;
            app.RecordSliceValuesTemp = app.RecordSliceValues;
            set(app.SampleSliceGrid,'ColumnSizes',[-1 -4],'RowSizes',[-1 -1 -1]);
            RecordSliceValuesCB(app) % run callback to execute expressions and populate values
            % ISO Values To Record panel
            app.SampleISOGrid = uiextras.Grid('Parent',app.RecordPanel(2,2),'Spacing',2,'Padding',2);
            app.ISOCheckbox(1) = uicontrol('Style','checkbox','Parent',app.SampleISOGrid,'String','Level 1','Value',app.ISOCheckboxVal(1),'Callback',@app.RecordWhichISOCB);
            app.ISOCheckbox(2) = uicontrol('Style','checkbox','Parent',app.SampleISOGrid,'String','Level 2','Value',app.ISOCheckboxVal(2),'Callback',@app.RecordWhichISOCB);
            app.RecordEmpty(5) = uiextras.Empty('Parent',app.SampleISOGrid);
            app.ISOEdit(1) = uicontrol('Style','edit','Parent',app.SampleISOGrid,'String',app.RecordISOValuesExpr{1},'BackgroundColor',[1 1 1],'Callback',@app.RecordISOValuesCB);
            app.ISOEdit(2) = uicontrol('Style','edit','Parent',app.SampleISOGrid,'String',app.RecordISOValuesExpr{2},'BackgroundColor',[1 1 1],'Callback',@app.RecordISOValuesCB);
            app.RecordISOValuesTemp = app.RecordISOValues;
            app.RecordISOValuesExprTemp = app.RecordISOValuesExpr;
            app.RecordEmpty(6) = uiextras.Empty('Parent',app.SampleISOGrid);
            set(app.SampleISOGrid,'ColumnSizes',[-1 -4],'RowSizes',[-1 -1 -1]);
            RecordISOValuesCB(app) % run callback to execute expressions and populate values
            
            % if imported data cannot be animated by samples, turn off
            % those options
            if length(app.Data(1,1,1,:))==2 && all(all(all(app.Data(:,:,:,2))))==0
                set(findall(app.RecordButtonGroup(3),'String','Slice/ISO'),'Value',1)
                set(app.AnimateByRadio(1),'Enable','off')
                set(app.SamplesValues,'Enable','off')
            end
            
            % enable/disable components based on saved options
            AnimateByCB(app)
        end
        
        function FrameRateCB(app,~,~)
            % FrameRateCB(app) is the callback function for the frame rate
            % edit box in the Animation & Recording Options window opened
            % through Options -> Animation & Recording Options. It allows
            % the user to specify a frame rate to use, which must be a
            % positive integer.
            FR = str2double(get(app.FileOptionsEdit(2),'String'));
            if round(FR)~=FR || FR <= 0
                warning('Frame Rate must be positive integer, reseting to previously saved value');
                BeepOnWarn(app);
                set(app.FileOptionsEdit(2),'String',num2str(app.RecordFrameRate));
            end
        end
        
        function RecordOKCB(app,~,~)
            % RecordOKCB(app) is the callback for the OK button in the
            % Recording Options window opened through the menu Options ->
            % Recording Options. It confirms the current options to use for
            % recording.
            
            % check if entered values for slices and isos are good
            if CompareSliceISOSizes(app)
                % update options
                app.RecordFileName = get(app.FileOptionsEdit(1),'String');
                app.RecordFileFormat = get(get(app.FileOptionsButtonGroup,'SelectedObject'),'String');
                app.RecordFrameRate = str2double(get(app.FileOptionsEdit(2),'String'));
                app.RecordMode = get(get(app.RecordButtonGroup(1),'SelectedObject'),'String');
                app.RecordBy = get(get(app.RecordButtonGroup(3),'SelectedObject'),'String');
                app.RecordArea = get(get(app.RecordButtonGroup(2),'SelectedObject'),'String');
                app.RepeatOn = get(app.RepeatOptions(1),'Value');
                app.RepeatTimes = str2double(get(app.RepeatOptions(3),'String'));
                app.RecordSampleValues = app.RecordSampleValuesTemp;
                app.RecordSampleValuesExpr = app.RecordSampleValuesExprTemp;
                app.RecordSliceValues = app.RecordSliceValuesTemp;
                app.RecordSliceValuesExpr = app.RecordSliceValuesExprTemp;
                app.RecordISOValues = app.RecordISOValuesTemp;
                app.RecordISOValuesExpr = app.RecordISOValuesExprTemp;
                app.SliceCheckboxVal = [get(app.SliceCheckbox(1),'Value'),get(app.SliceCheckbox(2),'Value'),get(app.SliceCheckbox(3),'Value')];
                app.ISOCheckboxVal = [get(app.ISOCheckbox(1),'Value'),get(app.ISOCheckbox(2),'Value')];
                delete(app.RecordOptionFigure) % delete recording options window
                EnableDisableFig(app,'enable') % renable main app window
                % reset slices and isos
                for ii = 1:length(app.SliceHandle)
                    if app.SliceCheckboxVal(ii) == 0 && ishandle(app.SliceHandle(ii))
                        delete(app.SliceHandle(ii))
                    end
                end
                for ii = 1:length(app.PatchHandle)
                    if app.ISOCheckboxVal(ii) == 0 && ishandle(app.PatchHandle(ii))
                        delete(app.PatchHandle(ii))
                    end
                end
                RedrawSlice(app)
                RedrawISO(app)
            else
                warning('Nonscalar values for the slices and isosurfaces must have the same size.');
                BeepOnWarn(app);
            end
        end
        
        function RecordCancelCB(app,~,~)
            % RecordCancelCB(app) is the callback function for closing the
            % Animation & Recording Options window without hitting OK.
            delete(app.RecordOptionFigure)
            EnableDisableFig(app,'enable')
        end
        
        function RecordHelpCB(app,~,~) %#ok
            % RecordHelpCB(app) is the callback for the Help button in the
            % Animation & Recording Options window opened through the menu
            % Options -> Animation & Recording Options. It opens the help
            % documentation for more information on options for recording.
            web([pwd filesep 'html' filesep 'HelpDocFile.html']);
        end
        
        function AnimateByCB(app,~,~)
            % AnimateByCB(app) is the callback function for the button
            % group for Samples vs Sliders. Animating by Samples will be
            % similar to hitting the PLAY button. Animating by Sliders will
            % allow for recording moving sliders. Different options in the
            % Recording Options window will be disabled or enabled based on
            % the selection.
            if strcmp(get(get(app.RecordButtonGroup(3),'SelectedObject'),'String'),'Samples')
                set([app.SliceEdit(1),app.SliceEdit(2),app.SliceEdit(3),...
                    app.ISOEdit(1),app.ISOEdit(2)],'Enable','off');
                set(app.SamplesValues,'Enable','on');
            else
                set([app.SliceEdit(1),app.SliceEdit(2),app.SliceEdit(3),...
                    app.ISOEdit(1),app.ISOEdit(2)],'Enable','on');
                set(app.SamplesValues,'Enable','off');
                % when recording sliders, enable/disable options according
                % to checkboxes
                RecordWhichSliceCB(app);
                RecordWhichISOCB(app);
            end
        end
        
        function RepeatTimesCB(app,~,~)
            % RepeatTimesCB(app) is the callback function for the edit box
            % in the Animation & Recording Options window for specifying
            % the number of times to repeat a recording. The value must be
            % an integer of value zero or greater.
            RT = str2double(get(app.RepeatOptions(3),'String'));
            if round(RT)~=RT || RT < 0
                warning('Number of times to repeat must be integer greater or equal to zero');
                BeepOnWarn(app);
                set(app.RepeatOptions(3),'String',num2str(app.RepeatTimes));
            end
        end
        
        function RecordSamplesCB(app,~,~)
            % RecordSamplesCB(app) is the callback for the edit box in
            % Recording Options window for specifying which samples to use
            % in animation.
            slength = length(app.Data(1,1,1,:));
            try
                rs = eval(get(app.SamplesValues,'String'));
                if any(rs<1) || any(rs>slength) || any(round(rs)~=rs)
                    warning('x-slice values are out of range');BeepOnWarn(app);
                    set(app.SamplesValues,'String',app.RecordSampleValuesExprTemp);
                else
                    app.RecordSampleValuesTemp = rs;
                    app.RecordSampleValuesExprTemp = get(app.SamplesValues,'String');
                end
            catch
                warning('Invalid expression entered');BeepOnWarn(app);
                set(app.SamplesValues,'String',app.RecordSampleValuesExprTemp);
            end
        end
        
        function RecordWhichSliceCB(app,~,~)
            % RecordWhichSliceCB(app) is the callback function for the
            % checkboxes in the Animation & Recording Options window for
            % specifying which Slices to use in animating.
            usingsliders = strcmp(get(get(app.RecordButtonGroup(3),'SelectedObject'),'String'),'Slice/ISO');
            if get(app.SliceCheckbox(1),'Value')==0 || usingsliders==0
                set(app.SliceEdit(1),'Enable','off');
            else
                set(app.SliceEdit(1),'Enable','on');
            end
            if get(app.SliceCheckbox(2),'Value')==0 || usingsliders==0
                set(app.SliceEdit(2),'Enable','off');
            else
                set(app.SliceEdit(2),'Enable','on');
            end
            if get(app.SliceCheckbox(3),'Value')==0 || usingsliders==0
                set(app.SliceEdit(3),'Enable','off');
            else
                set(app.SliceEdit(3),'Enable','on');
            end
        end
        
        function RecordSliceValuesCB(app,~,~)
            % RecordSliceValuesCB(app) is the callback function for the
            % edit boxes in the Animation & Recording Options window for
            % specifying which slice values to use in the recording.
            ylength = length(app.Data(:,1,1,1));xlength = length(app.Data(1,:,1,1));zlength = length(app.Data(1,1,:,1)); %#ok
            minx = min(app.DimValues{1}(:));maxx = max(app.DimValues{1}(:));
            miny = min(app.DimValues{2}(:));maxy = max(app.DimValues{2}(:));
            minz = min(app.DimValues{3}(:));maxz = max(app.DimValues{3}(:));
            % if first slice checkbox is on
            if get(app.SliceCheckbox(1),'Value')
                % see if expression entered for slice values is valid
                try
                    rx = eval(get(app.SliceEdit(1),'String'));
                    if any(rx<minx) || any(rx>maxx)
                        warning('x-slice values are out of range');BeepOnWarn(app);
                        set(app.SliceEdit(1),'String',app.RecordSliceValuesExprTemp{1});
                    else
                        app.RecordSliceValuesTemp{1} = rx;
                        app.RecordSliceValuesExprTemp{1} = get(app.SliceEdit(1),'String');
                    end
                catch
                    warning('Invalid expression entered');BeepOnWarn(app);
                    set(app.SliceEdit(1),'String',app.RecordSliceValuesExprTemp{1});
                end
            end
            % if second slice checkbox is on
            if get(app.SliceCheckbox(2),'Value')
                % see if expression entered for slice values is valid
                try
                    ry = eval(get(app.SliceEdit(2),'String'));
                    if any(ry<miny) || any(ry>maxy)
                        warning('y-slice values are out of range');BeepOnWarn(app);
                        set(app.SliceEdit(2),'String',app.RecordSliceValuesExprTemp{2});
                    else
                        app.RecordSliceValuesTemp{2} = ry;
                        app.RecordSliceValuesExprTemp{2} = get(app.SliceEdit(2),'String');
                    end
                catch
                    warning('Invalid expression entered');BeepOnWarn(app);
                    set(app.SliceEdit(2),'String',app.RecordSliceValuesExprTemp{2});
                end
            end
            % if third slice checkbox is on
            if get(app.SliceCheckbox(3),'Value')
                % see if expression entered for slice values is valid
                try
                    rz = eval(get(app.SliceEdit(3),'String'));
                    if any(rz<minz) || any(rz>maxz)
                        warning('z-slice values are out of range');BeepOnWarn(app);
                        set(app.SliceEdit(3),'String',app.RecordSliceValuesExprTemp{3});
                    else
                        app.RecordSliceValuesTemp{3} = rz;
                        app.RecordSliceValuesExprTemp{3} = get(app.SliceEdit(3),'String');
                    end
                catch
                    warning('Invalid expression entered');BeepOnWarn(app);
                    set(app.SliceEdit(3),'String',app.RecordSliceValuesExprTemp{3});
                end
            end
        end
        
        function RecordWhichISOCB(app,~,~)
            % RecordWhichISOCB(app) is the callback function for the
            % checkboxes in the Animation & Recording Options window for
            % specifying which ISO Levels to use in animating.
            usingsliders = strcmp(get(get(app.RecordButtonGroup(3),'SelectedObject'),'String'),'Slice/ISO');
            if get(app.ISOCheckbox(1),'Value')==0 || usingsliders==0
                set(app.ISOEdit(1),'Enable','off');
            else
                set(app.ISOEdit(1),'Enable','on');
            end
            if get(app.ISOCheckbox(2),'Value')==0 || usingsliders==0
                set(app.ISOEdit(2),'Enable','off');
            else
                set(app.ISOEdit(2),'Enable','on');
            end
        end
        
        function RecordISOValuesCB(app,~,~)
            % RecordISOValuesCB(app) is the callback function for the edit
            % boxes in the Animation & Recording Options window for
            % specifying which ISO values to use in the recording.
            mind = min(app.Data(:));maxd = max(app.Data(:));
            % if the first ISO checkbox is on
            if get(app.ISOCheckbox(1),'Value')
                % see if expression entered for ISO values is valid
                try
                    rd = eval(get(app.ISOEdit(1),'String'));
                    if any(rd<mind) || any(rd>maxd)
                        warning('iso level 1 values are out of range');BeepOnWarn(app);
                        set(app.ISOEdit(1),'String',app.RecordISOValuesExprTemp{1});
                    else
                        app.RecordISOValuesTemp{1} = rd;
                        app.RecordISOValuesExprTemp{1} = get(app.ISOEdit(1),'String');
                    end
                catch
                    warning('Invalid expression entered');BeepOnWarn(app);
                    set(app.ISOEdit(1),'String',app.RecordISOValuesExprTemp{1});
                end
            end
            % if second ISO checkbox is on
            if get(app.ISOCheckbox(2),'Value')
                % see if expression entered for ISO values is valid
                try
                    rd2 = eval(get(app.ISOEdit(2),'String'));
                    if any(rd2<mind) || any(rd2>maxd)
                        warning('x-slice values are out of range');BeepOnWarn(app);
                        set(app.ISOEdit(2),'String',app.RecordISOValuesExprTemp{2});
                    else
                        app.RecordISOValuesTemp{2} = rd2;
                        app.RecordISOValuesExprTemp{2} = get(app.ISOEdit(2),'String');
                    end
                catch
                    warning('Invalid expression entered');BeepOnWarn(app);
                    set(app.ISOEdit(2),'String',app.RecordISOValuesExprTemp{2});
                end
            end
        end
        
        function allgood = CompareSliceISOSizes(app)
            % allgood = CompareSliceISOSizes(app) is a helper function that
            % checks if the non-scalar expressions entered in the Animation
            % & Recording Options window for the slice and iso values to
            % animate over are of equal size.
            allgood = 1;
            for ii = 1:3 % check slices
                if get(app.SliceCheckbox(ii),'Value')==1 % only if that slice is on
                    for jj = 1:2 % check isos
                        if get(app.ISOCheckbox(jj),'Value')==1 % only if that iso is on
                            % if any of them have a nonscalar length that is not
                            % the same as the other nonscalar entries, something is
                            % wrong
                            if length(app.RecordSliceValuesTemp{ii}) ~= length(app.RecordISOValuesTemp{jj}) && ...
                                    length(app.RecordSliceValuesTemp{ii}) ~= 1 && ...
                                    length(app.RecordISOValuesTemp{jj}) ~=1
                                allgood = 0;
                            end
                        end
                    end
                end
            end
        end
        
        function EqualAxesCB(app,~,~)
            % EqualAxesCB(app) is the callback for the Options -> Equalize
            % Axes menu option. This toggles whether the aspect ratio is
            % set so that the data units are the same in every direction
            if strcmp(get(app.OptionsMenu(9),'Checked'),'off')
                % if it's off when clicked, turn it on and set axis equal
                set(app.OptionsMenu(9),'Checked','on')
                axis([app.SliceAxis app.ISOAxis],'equal')
                % the bounds need to be reset for ISO axis
                minx = min(app.DimValues{1}(:));maxx = max(app.DimValues{1}(:));
                miny = min(app.DimValues{2}(:));maxy = max(app.DimValues{2}(:));
                minz = min(app.DimValues{3}(:));maxz = max(app.DimValues{3}(:));
                set(app.ISOAxis,'XLim',[minx maxx],'YLim',[miny maxy],'ZLim',[minz maxz]);
            else % if it's on, turn it off and return axis to normal
                set(app.OptionsMenu(9),'Checked','off')
                axis([app.SliceAxis app.ISOAxis],'normal')
            end
        end
        
        function HelpCB(~,~,~)
            % HelpCB() is the callback for the Help menu option. It opens
            % the MATLAB web browser to an html help page for this app.
            web([filesep 'html' filesep 'HelpDocFile.html']);
        end
        
        function ColorbarCB(app,~,~)
            % ColorbarCB(app) is the callback for the Options -> Colorbar
            % menu option. It toggles whether there is a colorbar for
            % SliceAxis.
            if strcmp(get(app.OptionsMenu(12),'Checked'),'off')
                % if it's not on when clicked, turn it on and create the
                % colorbar
                set(app.OptionsMenu(12),'Checked','on')
                app.SliceColorbar = colorbar('peer',app.SliceAxis);
                if strcmp(get(app.OptionsMenu(17),'Checked'),'on')
                    % only turn on colorbar for isosurfaces if dynamic iso
                    % coloring is on so as not to confusing meaning of
                    % default isosurface colors
                    app.ISOColorbar = colorbar('peer',app.ISOAxis);
                end
            else % if it's on, turn it off and delete colorbar
                set(app.OptionsMenu(12),'Checked','off');
                if ishandle(app.SliceColorbar)
                    delete(app.SliceColorbar)
                end
                if ishandle(app.ISOColorbar)
                    delete(app.ISOColorbar)
                end
            end
        end
        
        function BeepCB(app,~,~)
            % BeepCB(app) is the callback for the Options -> Beep on
            % Warning options. This toggles whether MATLAB will beep when
            % invalid data or settings are used.
            if strcmp(get(app.OptionsMenu(16),'Checked'),'on')
                set(app.OptionsMenu(16),'Checked','off')
            else
                set(app.OptionsMenu(16),'Checked','on')
            end
        end
        
        function LinkRotationCB(app,~,~)
            % LinkRotationCB(app) is the callback for the Options -> Link
            % Rotation menu option. It toggles whether the viewing azimuth
            % and elevation of both SliceAxis and ISOAxis are the same, and
            % allows the rotate tool to rotate both at the same time.
            if strcmp(get(app.OptionsMenu(8),'Checked'),'off')
                % if it's off, turn it on and link axes rotation
                set(app.OptionsMenu(8),'Checked','on');
                setappdata(app.SliceAxis,'graphics_linkprop',linkprop([app.SliceAxis app.ISOAxis],{'CameraPosition','CameraUpVector'}));
            else % if it's on, turn it off and decouple axes rotation
                set(app.OptionsMenu(8),'Checked','off');
                setappdata(app.SliceAxis,'graphics_linkprop',[]);
            end
        end
        
        function SetFontSizeCB(app,~,~)
            % SetFontSizeCB(app) is the callback for the Options -> Set
            % Font Size menu option. It creates a dialog box that asks the
            % user to enter a value for the font size used in the edit
            % boxes. If the entered value is invalid it warns and keeps the
            % original value.
            prompt = {'Enter font size for edit boxes'};
            answer = inputdlg(prompt,'Set Font Size',1,{num2str(app.FontSize,15)});
            fontans = str2double(answer{1});
            try
                validateattributes(fontans,{'double'},{'scalar','nonempty','>',0});
                app.FontSize = fontans;
                set([app.SampleNum, app.TimeNum, app.DelayText, app.SliceValue,...
                    app.ISOLevelText, app.ISOAlphaText],'FontSize',app.FontSize);
            catch
                warning('Invalid font size entered');BeepOnWarn(app);
            end
        end
        
        function SetAxesNamesCB(app,~,~)
            % SetAxesNamesCB(app) is the callback function for the Options
            % -> Set Axes Names option. It opens a dialog box that allows
            % the user to specify custom axes names.
            prompt = {'Enter x-axis label','Enter y-axis label','Enter z-axis label'};
            xl = get(app.SliceAxis,'XLabel');yl = get(app.SliceAxis,'YLabel');
            zl = get(app.SliceAxis,'ZLabel');
            def = {get(xl,'String'),get(yl,'String'),get(zl,'String')};
            answer = inputdlg(prompt,'Set Font Size',1,def);
            if ~isempty(answer)
                % set the axes names
                set(xl,'String',answer{1});set(yl,'String',answer{2});set(zl,'String',answer{3});
                xlabel(app.ISOAxis,answer{1});ylabel(app.ISOAxis,answer{2});zlabel(app.ISOAxis,answer{3});
                % store names in UserData for use in custom data tip
                set(app.SliceHandle,'UserData',answer);
                set(app.PatchHandle,'UserData',answer);
            end
        end
        
        function SetColormapCB(app,~,~)
            % SetColormapCB(app) is the callback for the Options -> Set
            % Colormap menu option. It opens a list box that lets the user
            % import a variable from the workspace, invert the current
            % colormap, or edit the current colormap. For importing a
            % colormap, another dialog box is opened that allows the user
            % to specify a variable in the base workspace that contains the
            % colormap.
            list = {'jet','hsv','hot','cool','spring','summer','autmn','winter','gray','bone','copper','pink','parula','Invert current colormap','Import from workspace','Open Colormap Editor'};
            [answer,ok] = listdlg('ListString',list,'SelectionMode','single','Name','Select Colormap','PromptString','Select a colormap');
            if ok
                switch answer
                    case 1; set(app.Figure,'ColorMap',jet); close(gcf);
                    case 2; set(app.Figure,'ColorMap',hsv); close(gcf);
                    case 3; set(app.Figure,'ColorMap',hot); close(gcf);
                    case 4; set(app.Figure,'ColorMap',cool); close(gcf);
                    case 5; set(app.Figure,'ColorMap',spring); close(gcf);
                    case 6; set(app.Figure,'ColorMap',summer); close(gcf);
                    case 7; set(app.Figure,'ColorMap',autumn); close(gcf);
                    case 8; set(app.Figure,'ColorMap',winter); close(gcf);
                    case 9; set(app.Figure,'ColorMap',gray); close(gcf);
                    case 10; set(app.Figure,'ColorMap',bone); close(gcf);
                    case 11; set(app.Figure,'ColorMap',copper); close(gcf);
                    case 12; set(app.Figure,'ColorMap',pink); close(gcf);
                    case 13; set(app.Figure,'ColorMap',parula); close(gcf);
                    case 14 % invert current colormap
                        set(app.Figure,'Colormap',flip(get(app.Figure,'Colormap'),1));
                    case 15 % custom colormap
                        text = {'Enter colormap variable from workspace'};
                        checkvars = @(in) ismatrix(in) && length(in(1,:))==3 ...
                            && max(in(:)) <= 1 && min(in(:)) >= 0;
                        customcm = uigetvariables(text,'ValidationFcn',checkvars);
                        if ~isempty(customcm)
                            try
                                set(app.Figure,'ColorMap',customcm{1});
                            catch
                                warning('Invalid colormap variable used');BeepOnWarn(app);
                            end
                        end
                    case 16 % open colormap editor
                        colormapeditor(app.Figure)
                end
                RedrawISO(app)
            end
        end
        % load data functions
        function ImportWorkspaceCB(app,~,~)
            % ImportWorkspaceCB(app) is the callback for the Import Data ->
            % Import from Workspace menu option. This opens up a dialog to
            % select variables from the workspace using the uigetvariables
            % command. It checks if selected variables are of proper size,
            % and issues warning and uses defaults if they are not.
            
            % strings for use in variable selection popup
            intro = ['Volumetric Data Explorer can accept data in a variety of formats. '...
                'For full details, see the Help menu in the app menubar.'];
            text = {'Select data values (1-D or 3-D or 4-D):';
                'X values (1-D or 3-D) (optional):';
                'Y values (1-D or 3-D) (optional):';
                'Z values (1-D or 3-D) (optional):';
                'Time values (1-D) (optional):'};
            % custom validation functions so user can only select
            % appropriately sized variables in popup dialog
            checkvars_data = @(in) (ndims(in)==4 || ndims(in)==3 || (isvector(in) && size(in,1)>0)) && isa(in,'double');
            checkvars_dim = @(in) ((isvector(in) && size(in,1)>0) || ndims(in)==3) && isa(in,'double');
            checkvars_time = @(in) (isvector(in) && size(in,1)>0) && isa(in,'double');
            checkvars = {checkvars_data;checkvars_dim;checkvars_dim;checkvars_dim;checkvars_time};
            % get data
            [data,names] = uigetvariables(text,'ValidationFcn',checkvars,'Introduction',intro);
            % check if selected variables are valid, and load if yes
            if ~isempty(data)
                ValidateData(app,data,names)
            end
        end
        
        function ValidateData(app,data,names)
            % ValidateData(app) checks if user selected inputs are valid,
            % and if so uses the valid ones. If any of them are invalid or
            % not selected by the user, it will warn and use defaults.
            % 'data' is a cell array returned by uigetvariables(), with the
            % plotted values first, followed by x, y, z and t if they were
            % imported. 'names' is a cell array with the names of the
            % variables selected from the workspace.
            
            % Data must be nonempty, and all its entries must be of double
            % data type.
            valinput = ~isempty(data) && all(cellfun(@(x)isa(x,'double'),data));
            % If selected v data is gridded, then it must be 3D or 4D.
            valmat =  (ndims(data{1}) == 3 || ndims(data{1}) == 4);
            % If select v is a vector, then x, y, and z data must also be
            % entered and be vectors of the same length.
            valvect = (isvector(data{1}) && isvector(data{2}) && isvector(data{3}) && isvector(data{4}))...
                && (length(data{1}) == length(data{2}) && ...
                length(data{1}) == length(data{3}) && ...
                length(data{1}) == length(data{4}));
            % Combine the three
            valdata = valinput && (valmat || valvect);
            if ~valdata
                % if data to be visualized is not of right dimensions and
                % type, set it to empty
                data{1} = [];
            end
            if ~isempty(data{1}) && ~isvector(data{1})
                % if data dimensions are valid and pre-gridded
                app.Data = data{1}; % store data
                SetupGriddedData(app,data,names)
            elseif ~isempty(data{1}) && isvector(data{1}) % if data dimensions were not pre-gridded
                SetupUngriddedData(app,data,names)
            else % otherwise selected data was invalid, warn and do not setup app
                warning('Invalid or missing input data. Check sizes and data types. See app''s Help Documentation for more information.');BeepOnWarn(app);
            end
        end
        
        function SetupGriddedData(app,data,names)
            % SetupGriddedData(app) is a helper function that prepares
            % imported data that is already arranged in a grid. 'data' is a
            % cell array returned by uigetvariables(), with the plotted
            % values first, followed by x, y, z and t if they were
            % imported. 'names' is a cell array with the names of the
            % variables selected from the workspace.
            try
                % This try statement checks if the X, Y, and Z data
                % selected were all 3D matricies of proper size, or if
                % all were left blank. If any were 3D matricies of
                % improper size, it will issue a warning and use
                % default axis values. If any were not 3D matricies it
                % will error and move to the catch block.
                validmatricies = ((size(data{2}) == size(app.Data(:,:,:,1)))...
                    & (size(data{3}) == size(app.Data(:,:,:,1)))...
                    & (size(data{4}) == size(app.Data(:,:,:,1))))...
                    | (isempty(data{2}) && isempty(data{3}) && isempty(data{4}));
                if validmatricies
                    app.DimValues{1} = data{2};xlabel(app.SliceAxis,names{2});xlabel(app.ISOAxis,names{2});
                    app.DimValues{2} = data{3};ylabel(app.SliceAxis,names{3});ylabel(app.ISOAxis,names{3});
                    app.DimValues{3} = data{4};zlabel(app.SliceAxis,names{4});zlabel(app.ISOAxis,names{4});
                else
                    warning('Invalid X, Y and/or Z data selected, using defaults');BeepOnWarn(app);
                    [app.DimValues{1},app.DimValues{2},app.DimValues{3}]=...
                        meshgrid(1:length(app.Data(1,:,1,1)),1:length(app.Data(:,1,1,1)),1:length(app.Data(1,1,:,1)));
                    xlabel(app.SliceAxis,'X');ylabel(app.SliceAxis,'Y');zlabel(app.SliceAxis,'Z');
                    xlabel(app.ISOAxis,'X');ylabel(app.ISOAxis,'Y');zlabel(app.ISOAxis,'Z');
                end
            catch
                % the catch block checks if the X, Y, and Z variables
                % selected were all 1D vectors of proper sizes. If not,
                % it will use default values. If the variables were not
                % empty and are invalid, it will issue a warning.
                validvectors = (isvector(data{2}) && length(data{2})==length(app.Data(1,:,1,1)))...
                    && (isvector(data{3}) && length(data{3})==length(app.Data(:,1,1,1)))...
                    && (isvector(data{4}) && length(data{4})==length(app.Data(1,1,:,1)));
                if validvectors
                    app.DimValues{1} = data{2};app.DimValues{2} = data{3};app.DimValues{3} = data{4};
                    xlabel(app.SliceAxis,names{2});ylabel(app.SliceAxis,names{3});zlabel(app.SliceAxis,names{4});
                    xlabel(app.ISOAxis,names{2});ylabel(app.ISOAxis,names{3});zlabel(app.ISOAxis,names{4});
                else
                    if ~isempty(data{2}) || ~isempty(data{3}) || ~isempty(data{4})
                        warning('Invalid X, Y and/or Z data, using defaults');BeepOnWarn(app);
                    end
                    [app.DimValues{1},app.DimValues{2},app.DimValues{3}]=...
                        meshgrid(1:length(app.Data(1,:,1,1)),1:length(app.Data(:,1,1,1)),1:length(app.Data(1,1,:,1)));
                    xlabel(app.SliceAxis,'X');ylabel(app.SliceAxis,'Y');zlabel(app.SliceAxis,'Z');
                    xlabel(app.ISOAxis,'X');ylabel(app.ISOAxis,'Y');zlabel(app.ISOAxis,'Z');
                end
            end
            % check if the selected time vector is of the right size
            validtime = isvector(data{5}) && length(data{5})==length(app.Data(1,1,1,:));
            if validtime
                app.TimeValues = data{5};
            else
                if ~isempty(data{5})
                    warning('Invalid time data, using defaults');BeepOnWarn(app);
                end
                app.TimeValues = [];
            end
            % setup ISO Levels
            SetupISOLevels(app)
            % setup app with imported data
            LoadDataSetup(app)
        end
        
        function SetupUngriddedData(app,data,names)
            % SetupUngriddedData(app,data,names) is used when importing
            % un-gridded data. It transforms it into whichever n-D shape is
            % appropriate. 'data' is a cell array returned by
            % uigetvariables(), with the plotted values first, followed by
            % x, y, z and t if they were imported. 'names' is a cell array
            % with the names of the variables selected from the workspace.
            
            % Here we use the full data set for reference, so the
            % axes limits don't change during animation
            numpoints = 30;
            [app.DimValues{1},app.DimValues{2},app.DimValues{3}] = meshgrid(linspace(min(data{2}),max(data{2}),numpoints),...
                linspace(min(data{3}),max(data{3}),numpoints),linspace(min(data{4}),max(data{4}),numpoints));
            if length(data{5}) == length(data{1})
                % if time data selected is of appropriate size
                app.TimeValues = sort(unique(data{5}));
                app.Data(length(app.DimValues{1}),length(app.DimValues{2}),length(app.DimValues{3}),length(app.TimeValues)) = 0;
                for ii = 1:length(app.TimeValues)
                    idx = data{5} == app.TimeValues(ii);
                    x = data{2}(idx); y = data{3}(idx); z = data{4}(idx);
                    v = data{1}(idx);
                    if verLessThan('matlab','8.1')
                        % scatteredInterpolant was introduced in R2013a
                        % (MATLAB 8.1) so on earlier versions must use
                        % TriScatteredInterp
                        F = TriScatteredInterp(x,y,z,v); %#ok
                    else
                        F = scatteredInterpolant(x,y,z,v);
                    end
                    app.Data(:,:,:,ii) = F(app.DimValues{1},app.DimValues{2},app.DimValues{3});
                end
                xlabel(app.SliceAxis,names{2});ylabel(app.SliceAxis,names{3});zlabel(app.SliceAxis,names{4});
                xlabel(app.ISOAxis,names{2});ylabel(app.ISOAxis,names{3});zlabel(app.ISOAxis,names{4});
                % setup ISO Levels
                SetupISOLevels(app)
                % setup app with imported data
                LoadDataSetup(app)
            elseif isempty(data{5}) % if no time data was selected
                app.Data(length(app.DimValues{1}),length(app.DimValues{2}),length(app.DimValues{3})) = 0;
                if verLessThan('matlab','8.1')
                    % scatteredInterpolant was introduced in R2013a
                    % (MATLAB 8.1) so on earlier versions must use
                    % TriScatteredInterp
                    F = TriScatteredInterp(data{2},data{3},data{4},data{1}); %#ok
                else
                    F = scatteredInterpolant(data{2},data{3},data{4},data{1});
                end
                app.Data = F(app.DimValues{1},app.DimValues{2},app.DimValues{3});
                xlabel(app.SliceAxis,names{2});ylabel(app.SliceAxis,names{3});zlabel(app.SliceAxis,names{4});
                xlabel(app.ISOAxis,names{2});ylabel(app.ISOAxis,names{3});zlabel(app.ISOAxis,names{4});
                % setup ISO Levels
                SetupISOLevels(app)
                % setup app with imported data
                LoadDataSetup(app)
            else
                % time of improper size, need to warn
                warning(['Selected variable for time values must be the same size as variables for x, y, z, and v.' ...
                    'Please reselect data to import.'])
            end
        end
        
        function SetupISOLevels(app)
            % SetupISOLevels(app) is used to determine the initial ISO
            % levels to use when importing a new set of data.
            
            % take guess at ISO levels (this is before we expand 3D
            % only input data to 4D in LoadDataSetup, so we don't have
            % to worry about filtering out the added zeros)
            mind = min(app.Data(:));maxd = max(app.Data(:));
            diff = maxd-mind;
            % set the first iso level at 1/3 from min towards max
            set(app.ISOLevelText(1),'String',num2str(min(app.Data(:))+diff/3,15));
            set(app.ISOSlider(1),'Min',mind,'Max',maxd,'Value',min(app.Data(:))+diff/3);
            % set the second iso level at 2/3 from min towards max
            set(app.ISOLevelText(2),'String',num2str(min(app.Data(:))+2*diff/3,15));
            set(app.ISOSlider(2),'Min',mind,'Max',maxd,'Value',min(app.Data(:))+2*diff/3);
        end
        
        function LoadEllipsoidCB(app,~,~)
            % LoadEllipsoidCB(app) is the callback for the Import Data ->
            % Example Data -> Oscillating Ellipsoid menu option. It creates
            % sample data for an ellipsoid whose semi-principal axes in the
            % X and Z directions vary with time according to
            % x^2/(5-4*cos(t)) + y^2 + z^2/(5+4*cos(t)).
            x = -10:1:10;t = 0:0.2:40;
            [xx,yy,zz] = meshgrid(x,x,x);
            data(length(x),length(x),length(x),length(t)) = 0;
            for ii=1:length(t)
                data(:,:,:,ii) = xx.^2/(5-4*cos(t(ii))) + yy.^2 + zz.^2/(5+4*cos(t(ii)));
            end
            app.Data = data;
            app.DimValues{1} = xx;
            app.DimValues{2} = yy;
            app.DimValues{3} = zz;
            app.TimeValues = t;
            set(app.ISOLevelText(1),'String','20');
            set(app.ISOSlider(1),'Min',min(app.Data(:)),'Max',max(app.Data(:)),'Value',20);
            set(app.ISOLevelText(2),'String','6');
            set(app.ISOSlider(2),'Min',min(app.Data(:)),'Max',max(app.Data(:)),'Value',6);
            set([app.SliceAxis,app.ISOAxis],'CLim',[min(app.Data(:)),max(app.Data(:))])
            LoadDataSetup(app);
        end
        
        function LoadFluidFlowCB(app,~,~)
            % LoadFluidFlowCB(app) is the callback for the Import Data ->
            % Example Data -> Fluid Flow menu option. It uses a modified
            % version of MATLAB's built-in FLOW command to vary the flow
            % parameters A and nu over time.
            x = 0:0.5:10;y = -3:0.2:3;z = y;t = 0:0.1:10;
            [app.DimValues{1},app.DimValues{2},app.DimValues{3},app.Data] = flowmod(x,y,z,t);
            app.TimeValues = t;
            set(app.ISOLevelText(1),'String','4');
            set(app.ISOSlider(1),'Min',min(app.Data(:)),'Max',max(app.Data(:)),'Value',4);
            set(app.ISOLevelText(2),'String','5');
            set(app.ISOSlider(2),'Min',min(app.Data(:)),'Max',max(app.Data(:)),'Value',5);
            set(app.ISOAlphaText(2),'String','0.7');
            set([app.SliceAxis,app.ISOAxis],'CLim',[min(app.Data(:)),max(app.Data(:))])
            LoadDataSetup(app);
        end
        
        function LoadDataSetup(app)
            % LoadDataSetup(app) is run after importing a new data set. It
            % sets up the limits and values for the various sliders and
            % axes in the app, and creates visualizations of the first
            % sample.
            
            app.CurrentSample = 1;
            confs = load('DE_confs');
            setSettings(app,confs.conf_default); % set default settings
            EnableDisableFig(app,'enable'); % enable entire figure
            if ndims(app.Data) == 3 % if 3-D data (no animation dimension)
                app.Data(:,:,:,2) = 0; % add a 4th dimension so we can reuse same code (costs extra memory, planning to fix that in future)
                EnableDisableFig(app,'enable') % enables entire fig then turns off right parts
                set(app.OptionsMenu(10),'Checked','off','Enable','off'); % no animation, so disabling looping option
            end
            ylength = length(app.Data(:,1,1,1));
            xlength = length(app.Data(1,:,1,1));
            zlength = length(app.Data(1,1,:,1));
            tlength = length(app.Data(1,1,1,:));
            minx = min(app.DimValues{1}(:));maxx = max(app.DimValues{1}(:));
            miny = min(app.DimValues{2}(:));maxy = max(app.DimValues{2}(:));
            minz = min(app.DimValues{3}(:));maxz = max(app.DimValues{3}(:));
            mind = min(app.Data(:));maxd = max(app.Data(:));
            % set slider values
            set(app.GridComp{3,2},'Min',minx,'Max',maxx,'Value',maxx,'SliderStep',[.1/(xlength-1) 1/(xlength-1)]);
            set(app.GridComp{1,2},'Min',miny,'Max',maxy,'Value',maxy,'SliderStep',[.1/(ylength-1) 1/(ylength-1)]);
            set(app.GridComp{2,1},'Min',minz,'Max',maxz,'Value',minz,'SliderStep',[.1/(zlength-1) 1/(zlength-1)]);
            set(app.ISOAxis,'CLim',[mind maxd],'XLim',[minx maxx],'YLim',[miny maxy],'ZLim',[minz maxz]);
            set(app.SliceAxis,'CLim',[mind maxd],'XLim',[minx maxx],'YLim',...
                [miny maxy],'ZLim',[minz maxz]);
            set(app.PlaySlider,'Min',1,'Max',tlength,'Value',app.CurrentSample,'SliderStep',[1/(tlength-1) 10/(tlength-1)]);
            set(app.SampleNum,'String',num2str(app.CurrentSample,15));
            if isempty(app.TimeValues)
                % if valid time variable WAS NOT selected
                set(app.TimeNum,'Style','text'); % remove ability to edit time value
                app.CurrentTime = 0;
            else % if valid time variable WAS selected
                app.CurrentTime = app.TimeValues(app.CurrentSample);
                set(app.TimeNum,'Style','edit','String',num2str(app.CurrentTime,15)); % set current time and make it editable
            end
            % create Slice and Isosurface visualizations
            if ishandle(app.SliceHandle)
                % if slices exist, clear them for recreation in
                % RedrawSlice, otherwise their sizes won't be right
                delete(app.SliceHandle)
                app.SliceHandle = [];
            end
            view(app.SliceAxis,-30,50) % set viewing angle
            RedrawSlice(app)
            ISOTextCB(app)
            set(app.OptionsMenu(2),'Enable','on') % enable save/load configurations
            set(app.OptionsMenu(9),'Enable','on')
            % determine values for animating
            for ii = 1:3
                if app.SliceCheckboxVal(ii) == 1 && strcmp(app.RecordBy,'Slice/ISO')
                    app.RecordSliceValues{ii} = eval(app.RecordSliceValuesExpr{ii});
                end
            end
            for jj = 1:2
                if app.ISOCheckboxVal(jj) == 1 && strcmp(app.RecordBy,'Slice/ISO')
                    app.RecordISOValues{jj} = eval(app.RecordISOValuesExpr{jj});
                end
            end
            if strcmp(app.RecordBy,'Samples')
                app.RecordSampleValues = eval(app.RecordSampleValuesExpr);
            end
        end
        
        function BeepOnWarn(app)
            % BeepCB(app) is used to determine whether an audible beep is
            % issued when a custom warning occurs.
            if strcmp(get(app.OptionsMenu(16),'Checked'),'on');
                beep
            end
        end
        % control panel functions
        function DelayPlusCB(app,~,~)
            % DelayPlusCB(app) is the callback for the + push button in the
            % Delay panel. It adds 0.1 to the delay value.
            app.Delay = app.Delay + 0.1;
            set(app.DelayText,'String',num2str(app.Delay,15));
        end
        
        function DelayTextCB(app,~,~)
            % DelayTextCB(app) is the callback for the edit box in the
            % Delay panel. It displays the current delay time. If manually
            % set to a value less than zero, it will warn and reset the
            % delay value to the previous value.
            Delaynum = str2double(get(app.DelayText,'String'));
            if Delaynum >= 0
                app.Delay = Delaynum;
            else
                set(app.DelayText,'String',num2str(app.Delay,15));
                warning('Delay time must be positive');BeepOnWarn(app);
            end
        end
        
        function DelayMinusCB(app,~,~)
            % DelayMinusCB(app) is the callback for the - button in the
            % Delay panel. It will decrease the delay value by 0.1, unless
            % it would decrease it to less than 0.01 in which case it
            % rounds the delay value to zero.
            if (app.Delay - 0.1) < 0.01
                app.Delay = round(app.Delay);
                set(app.DelayText,'String',num2str(app.Delay,15));
            else
                app.Delay = app.Delay - 0.1;
                set(app.DelayText,'String',num2str(app.Delay,15));
            end
        end
        
        function SampleNumCB(app,~,~)
            % SampleNumCB(app) is the callback for the top edit box in the
            % Display panel. It allows the user to enter a sample number to
            % visualize, and is updated during animation. It must be an
            % integer in the range of samples taken. If it is outside the
            % range, it will be reset to the original value. If it is a
            % non-integer in range, it will be rounded to the nearest
            % integer.
            StNum = str2double(get(app.SampleNum,'String'));
            if (round(StNum)-StNum) ~= 0
                warning('Starting sample must be an integer value...entry rounded');BeepOnWarn(app);
            end
            if StNum >= get(app.PlaySlider,'Min') && StNum <= get(app.PlaySlider,'Max')
                app.CurrentSample = round(StNum);
                set(app.SampleNum,'String',app.CurrentSample);
                set(app.PlaySlider,'Value',app.CurrentSample);
                if ~isempty(app.TimeValues)
                    dispt = app.TimeValues;
                    set(app.TimeNum,'String',num2str(dispt(app.CurrentSample),15));
                end
                RedrawSlice(app);
                ISOTextCB(app);
            else
                warning('Starting Sample must be within number of samples');BeepOnWarn(app);
                set(app.SampleNum,'String',num2str(app.CurrentSample,15));
            end
        end
        
        function TimeNumCB(app,~,~)
            % TimeNumCB(app) is the callback for the bottom edit box in the
            % Display panel. If no valid time variable was selected when
            % importing data, it is disabled. Otherwise it allows the user
            % to select the time sample to display, and is updated during
            % animation. If the entered time value is not one that was
            % sampled at, then it will warn and round to the nearest one.
            CT = str2double(get(app.TimeNum,'String'));
            dispt = app.TimeValues;
            if ~ismember(CT,dispt)
                warning('Must be a time value sampled at, setting to nearest valid value');BeepOnWarn(app);
                [minl,lowt] = min(abs(dispt-CT));
                [minh,hight] = min(abs(dispt+CT));
                if minl < minh
                    app.CurrentSample = lowt;
                    app.CurrentTime = dispt(lowt);
                else
                    app.CurrentSample = hight;
                    app.CurrentTime = dispt(hight);
                end
            else
                app.CurrentTime = CT;
                app.CurrentSample = find(dispt == app.CurrentTime);
            end
            set(app.TimeNum,'String',num2str(app.CurrentTime,15));
            set(app.SampleNum,'String',num2str(app.CurrentSample,15));
            set(app.PlaySlider,'Value',app.CurrentSample);
            RedrawSlice(app);
            ISOTextCB(app);
        end
        
        function PlaySliderCB(app,~,~)
            % PlaySliderCB(app) is the callback for the Play Slider in the
            % Controls panel. It is disabled during animation. Outside of
            % animation it allows the user to scroll through visualizations
            % of time samples.
            app.CurrentSample = round(get(app.PlaySlider,'Value'));
            set(app.SampleNum,'String',num2str(app.CurrentSample,15));
            if ~isempty(app.TimeValues)
                dispt = app.TimeValues;
                set(app.TimeNum,'String',num2str(dispt(app.CurrentSample),15));
            end
            RedrawSlice(app);
            ISOTextCB(app);
        end
        
        function PlayButtonCB(app,~,~)
            % PlayButtonCB(app) is the callback for the PLAY toggle button
            % in the Controls panel. When it displays PLAY, clicking on it
            % will begin animation at the current sample, and change it to
            % a STOP button. When the STOP button is clicked, it halts the
            % animation and changes back to PLAY. While animation is in
            % progress, the Display panel and Play Slider are disabled.
            if get(app.PlayButton,'Value') == 1
                % disable animate button while PLAY button is on
                set(app.AnimateButton,'Enable','off')
                % if recording is not already on, disable it
                if get(app.RecordButton,'Value')==0
                    set(app.RecordButton,'Enable','off')
                end
                % only run if button was in PLAY state when clicked
                set(app.PlayButton,'String','STOP')
                % disable certain controls during animation
                set(app.SampleNum,'Enable','off')
                set(app.PlaySlider,'Enable','off')
                set(app.TimeNum,'Enable','off')
                set(app.ImportMenu(1),'Enable','off')
                % reset to beginning if starting animation from the last
                % sample
                if app.CurrentSample == length(app.Data(1,1,1,:))
                    app.CurrentSample = 1;
                end
                % animation
                PlayBySamples(app)
                % reenable controls after animation is complete
                set(app.SampleNum,'Enable','on')
                set(app.PlaySlider,'Enable','on')
                if ~isempty(app.TimeValues)
                    % only enable if valid time values exist
                    set(app.TimeNum,'Enable','on')
                end
                set(app.ImportMenu(1),'Enable','on') % import new data
                set(app.PlayButton,'String','PLAY','Value',0)
                % enable animate button back on
                set(app.AnimateButton,'Enable','on')
                % if recording had been disabled, reenable it
                if strcmp(get(app.RecordButton,'Enable'),'off')
                    set(app.RecordButton,'Enable','on')
                end
            end
        end
        
        function PlayBySamples(app)
            % PlayBySamples(app) animates the figures sample by samples, as
            % opposed to animating by changing slice location or iso
            % values.
            
            % while the play button is clicked and the time is less
            % than or equal to the max time, keep running
            t = app.CurrentSample;
            while t <= length(app.Data(1,1,1,:)) && get(app.PlayButton,'Value') == 1
                app.CurrentSample = t;
                set(app.SampleNum,'String',num2str(app.CurrentSample,15))
                set(app.PlaySlider,'Value',app.CurrentSample)
                if ~isempty(app.TimeValues)
                    set(app.TimeNum,'String',num2str(app.TimeValues(app.CurrentSample),15))
                end
                % update visualziations
                RedrawSlice(app)
                ISOTextCB(app)
                drawnow
                if strcmp(app.RecordMode,'Manual') % only record PLAY button if in Manual recording mode
                    RecordToFile(app) % save screencapture to file
                end
                pause(app.Delay) % introduce additional pause specified in app
                % if it is the last sample and loop animation option is
                % on, reset to beginning
                if t == length(app.Data(1,1,1,:)) && strcmp(get(app.OptionsMenu(10),'Checked'),'on')
                    t = 1;
                else
                    t = t + 1;
                end
            end
        end
        
        function AnimateButtonCB(app,~,~)
            % AnimationButtonCB(app) is called when using the ANIMATE
            % button in the Control panel. It will use the animation
            % options set in the Options menu -> Animation & Recording
            % Options window. When clicked it will begin the animation and
            % turn to a STOP ANIMATION button, which can be clicked to stop
            % animation early.
            EnableDisableFig(app,'disable')
            set(app.AnimateButton,'Enable','on','String','STOP ANIMATION','FontSize',12)
            if ~verLessThan('matlab','8.4') % font sizes changed in R2014b, version 8.4
                set(app.AnimateButton,'FontSize',10)
            end
            LoopCount = 1;
            if app.RepeatOn==1 % if repeating is on, get how many times
                LoopCount = app.RepeatTimes+1;
            end
            while LoopCount > 0 && get(app.AnimateButton,'Value')==1% if we need to repeat
                if strcmp(app.RecordBy,'Samples')
                    % if we are animating by samples
                    ii = 1;
                    while ii <= length(app.RecordSampleValues) && get(app.AnimateButton,'Value')==1
                        app.CurrentSample = app.RecordSampleValues(ii);
                        set(app.SampleNum,'String',num2str(app.CurrentSample,15));
                        set(app.PlaySlider,'Value',app.CurrentSample);
                        if ~isempty(app.TimeValues)
                            set(app.TimeNum,'String',num2str(app.TimeValues(app.CurrentSample),15));
                        end
                        % update visualziations
                        RedrawSlice(app);
                        RedrawISO(app);
                        RecordToFile(app); % save screencapture to file
                        pause(app.Delay); % introduce additional pause specified in app
                        ii = ii+1;
                    end
                elseif strcmp(app.RecordBy,'Slice/ISO')
                    ii = 1; 
                    numsamples = 1; 
                    if any(app.SliceCheckboxVal) % if slices are checked
                        SliceSliders = [app.GridComp{3,2},app.GridComp{1,2},app.GridComp{2,1}]; % [x,y,z] sliders
                        for jj = 1:3 % find length of slice values
                            if app.SliceCheckboxVal(jj) % only search if that slice is being animated
                                samplength = length(app.RecordSliceValues{jj});
                                if samplength > numsamples
                                    numsamples = samplength;
                                end
                                if samplength == 1 % if single value is used for constant slice, set slider value here
                                    set(SliceSliders(jj),'Value',app.RecordSliceValues{jj});
                                end
                            end
                        end
                    end
                    if any(app.ISOCheckboxVal) % if isosurfaces are checked
                        for jj = 1:2 % find length of iso values
                            if app.ISOCheckboxVal(jj) % only search if that isosurface is being animated
                                samplength = length(app.RecordISOValues{jj});
                                if samplength > numsamples
                                    numsamples = samplength;
                                end
                                if samplength == 1 % if single value is used for constant isosurface, set slider value here
                                    set(app.ISOSlider(jj),'Value',app.RecordISOValues{jj});
                                end
                            end
                        end
                    end
                    while ii <= numsamples && get(app.AnimateButton,'Value')==1
                        for kk = 1:3 % loop through x,y,z and set new slider values if appropriate
                            if app.SliceCheckboxVal(kk)==1 && length(app.RecordSliceValues{kk})>1
                                % if that slice box is checked and the values are non constant
                                set(SliceSliders(kk),'Value',app.RecordSliceValues{kk}(ii));
                            end
                        end
                        for kk = 1:2
                            if app.ISOCheckboxVal(kk)==1 && length(app.RecordISOValues{kk})>1
                                set(app.ISOSlider(kk),'Value',app.RecordISOValues{kk}(ii))
                            end
                        end
                        % update visualziations
                        RedrawSlice(app);
                        RedrawISO(app); % need to call this instead of RedrawISO so that new slider values get picked up
                        RecordToFile(app); % save screencapture to file
                        pause(app.Delay); % introduce additional pause specified in app
                        ii = ii+1;
                    end

                end
                LoopCount = LoopCount - 1;
            end
            EnableDisableFig(app,'enable')
            set(app.AnimateButton,'String','ANIMATE','FontSize',24,'Value',0)
            if ~verLessThan('matlab','8.4') % font sizes changed in R2014b, version 8.4
                set(app.AnimateButton,'FontSize',18)
            end
        end
        
        function RecordButtonCB(app,~,~)
            % RecordButtonCB(app) is used by the RECORD button in the
            % Control panel. Clicking record will open a connection to a
            % file and will then record any changes to the plots to the
            % file. The file used and other recording options are set in
            % the Options menu -> Animation & Recording Options window. Once
            % clicked it will turn into a STOP ANIMATION button, which when
            % clicked will close the connection to the file.
            if get(app.RecordButton,'Value')==1
                set(app.Figure,'Resize','off');
                set(app.RecordButton,'String','STOP RECORDING','FontSize',12);
                switch app.RecordFileFormat
                    case '.avi'
                        app.WriteObject = VideoWriter([app.RecordFileName app.RecordFileFormat]);
                    case '.mp4'
                        app.WriteObject = VideoWriter([app.RecordFileName app.RecordFileFormat],'MPEG-4');
                end
                app.WriteObject.FrameRate = app.RecordFrameRate;
                open(app.WriteObject);
            else
                close(app.WriteObject);
                set(app.Figure,'Resize','on');
                set(app.RecordButton,'String','RECORD','FontSize',24);
            end
        end
        % slice panel functions
        function RedrawSlice(app,~,~)
            % RedrawSlice(app) is the callback for the Slice panel sliders
            % and also creates the slice visualization. It uses a modified
            % version of MATLAB's built-in SLICE command.
            SliceLoc{1} = get(app.GridComp{3,2},'Value');
            set(app.SliceValue(1),'String',num2str(SliceLoc{1},15));
            if app.SliceCheckboxVal(1) == 0 % if this slice is disabled in Record & Animation Options window, turn it off
                SliceLoc{1} = [];
            end
            SliceLoc{2} = get(app.GridComp{1,2},'Value');
            set(app.SliceValue(2),'String',num2str(SliceLoc{2},15));
            if app.SliceCheckboxVal(2) == 0 % if this slice is disabled in Record & Animation Options window, turn it off
                SliceLoc{2} = [];
            end
            SliceLoc{3} = get(app.GridComp{2,1},'Value');
            set(app.SliceValue(3),'String',num2str(SliceLoc{3},15));
            if app.SliceCheckboxVal(3) == 0 % if this slice is disabled in Record & Animation Options window, turn it off
                SliceLoc{3} = [];
            end
            createslice = 0; % should we recreate slices
            if isempty(app.SliceHandle) % if empty then this is first data set so yes
                createslice = 1;
            else % otherwise check if one was deleted and needs to be recreated
                for ii = 1:length(app.SliceHandle) % if any of them have been deleted and turned back on
                    if ~ishandle(app.SliceHandle(ii)) && ~isempty(SliceLoc{ii})
                        createslice = 1; % recreate slices
                    end
                end
                if createslice % if we need to recreate slices
                    for ii = 1:length(app.SliceHandle)
                        if ishandle(app.SliceHandle(ii))
                            delete(app.SliceHandle(ii)) % delete existing ones so there won't be duplicates
                        end
                    end
                end
            end
            if createslice % check if they need to be recreated
                app.SliceHandle = slicemod(app.SliceAxis,app.DimValues{1},app.DimValues{2},app.DimValues{3},app.Data(:,:,:,app.CurrentSample),SliceLoc{1},SliceLoc{2},SliceLoc{3});
                set(app.SliceHandle,'FaceColor','interp','EdgeAlpha',0.2);
                xl = get(app.SliceAxis,'XLabel');yl = get(app.SliceAxis,'YLabel');zl = get(app.SliceAxis,'ZLabel');
                set(app.SliceHandle,'UserData',{get(xl,'String');get(yl,'String');get(zl,'String')});
            else % otherwise, update underlying data rather than recreate them
                UpdateSlice(app,SliceLoc{1},SliceLoc{2},SliceLoc{3});
            end
            % if in Manual recording mode, so if not in PLAY mode AND not in ANIMATE mode AND are in RECORD mode
            if strcmp(app.RecordMode,'Manual')
                RecordToFile(app); % save screencapture to file
            end
        end
        
        function UpdateSlice(app,xh,yh,zh)
            % UpdateSlice(app,xh,yh,zh) updates the slices without
            % recreating them. This avoids creating and deleting graphics
            % objects unnecessarily. xh, yh and zh are new slice locations.
            % It uses the same method for computing the color of slices as
            % the built-in |slice| command.
            % x-slice
            ii = 1;
            if ~isempty(xh) && app.SliceCheckboxVal(1)==1
                xs = size(get(app.SliceHandle(ii),'XData'));
                xd = xh*ones(xs(1),xs(2)); % new x values for x-slice
                yd = get(app.SliceHandle(ii),'YData'); % keep y values for x-slice
                zd = get(app.SliceHandle(ii),'ZData'); % keep z values for x-slice
                % new color data for x-slice, found using interp3 which is
                % used under the hood of the slice command
                vi = interp3(app.DimValues{1},app.DimValues{2},app.DimValues{3},app.Data(:,:,:,app.CurrentSample),xd,yd,zd);
                set(app.SliceHandle(ii),'XData',xd,'CData',vi);
                ii=ii+1;
            end
            % y-slice
            if ~isempty(yh) && app.SliceCheckboxVal(2) == 1
                xd = get(app.SliceHandle(ii),'XData');
                ys = size(get(app.SliceHandle(ii),'XData'));
                yd = yh*ones(ys(1),ys(2));
                zd = get(app.SliceHandle(ii),'ZData');
                vi = interp3(app.DimValues{1},app.DimValues{2},app.DimValues{3},app.Data(:,:,:,app.CurrentSample),xd,yd,zd);
                set(app.SliceHandle(ii),'YData',yd,'CData',vi);
                ii=ii+1;
            end
            % z-slice
            if ~isempty(zh) && app.SliceCheckboxVal(3) == 1
                xd = get(app.SliceHandle(ii),'XData');
                yd = get(app.SliceHandle(ii),'YData');
                zs = size(get(app.SliceHandle(ii),'ZData'));
                zd = zh*ones(zs(1),zs(2));
                vi = interp3(app.DimValues{1},app.DimValues{2},app.DimValues{3},app.Data(:,:,:,app.CurrentSample),xd,yd,zd);
                set(app.SliceHandle(ii),'ZData',zd,'CData',vi);
            end
        end
        
        function SliceValueCB(app,~,~)
            % SliceValueCB(app) is the callback for the edit boxes in the
            % Slice panel. It checks if the entered values are within the
            % axes ranges. If they are, it will update the values of the
            % Slice panel sliders and execute RedrawSlice. If not, it will
            % warn and return to the previous values.
            SV(3) = str2double(get(app.SliceValue(3),'String'));
            SV(2) = str2double(get(app.SliceValue(2),'String'));
            SV(1) = str2double(get(app.SliceValue(1),'String'));
            validrange = SV(3) >= get(app.GridComp{2,1},'Min') && SV(3) <= get(app.GridComp{2,1},'Max')...
                && SV(2) >= get(app.GridComp{1,2},'Min') && SV(2) <= get(app.GridComp{1,2},'Max')...
                && SV(1) >= get(app.GridComp{3,2},'Min') && SV(1) <= get(app.GridComp{3,2},'Max');
            if validrange
                set(app.GridComp{2,1},'Value',SV(3));
                set(app.GridComp{1,2},'Value',SV(2));
                set(app.GridComp{3,2},'Value',SV(1));
                RedrawSlice(app);
            else
                warning('Slice value must be within range of data');BeepOnWarn(app);
                set(app.SliceValue(3),'String',num2str(get(app.GridComp{2,1},'Value'),15));
                set(app.SliceValue(2),'String',num2str(get(app.GridComp{1,2},'Value'),15));
                set(app.SliceValue(1),'String',num2str(get(app.GridComp{3,2},'Value'),15));
            end
        end
        % ISO panel functions
        function ISOTextCB(app,~,~)
            % ISOTextCB is used when ISO level values are updated by
            % changing the values in the text boxes.
            ISOLevel(1) = str2double(get(app.ISOLevelText(1),'String'));
            ISOLevel(2) = str2double(get(app.ISOLevelText(2),'String'));
            for ii = 1:2
                if ISOLevel(ii) < min(app.Data(:)) || ISOLevel(ii) > max(app.Data(:)) % if iso level is not in right range
                    warning('ISO Level must be between the min and max values of the imported data');BeepOnWarn(app);
                    set(app.ISOLevelText(ii),'String',num2str(get(app.ISOSlider(ii),'Value'),15));
                else
                    set(app.ISOSlider(ii),'Value',ISOLevel(ii))
                end
            end
            RedrawISO(app)
        end
        
        function RedrawISO(app,~,~)
            % RedrawISO(app) is used when updating the ISO plots.
            ISOLevel(1) = get(app.ISOSlider(1),'Value');
            set(app.ISOLevelText(1),'String',num2str(ISOLevel(1),15));
            ISOAlpha(1) = str2double(get(app.ISOAlphaText(1),'String'));
            ISOLevel(2) = get(app.ISOSlider(2),'Value');
            set(app.ISOLevelText(2),'String',num2str(ISOLevel(2),15));
            ISOAlpha(2) = str2double(get(app.ISOAlphaText(2),'String'));
            if app.ISOCheckboxVal(1) == 0
                ISOAlpha(1) = -1; % if level 1 is checked off in recording options window, force it to disable
            end
            if app.ISOCheckboxVal(2) == 0
                ISOAlpha(2) = -1; % if level 2 is checked off in recording options window, force it to disable
            end
            for ii = 1:2
                if ISOAlpha(ii) > 0 && ISOAlpha(ii) <= 1 % if alpha in right range
                    surf1data = isosurface(app.DimValues{1},app.DimValues{2},app.DimValues{3},app.Data(:,:,:,app.CurrentSample),ISOLevel(ii));
                    if strcmp(get(app.OptionsMenu(17),'Checked'),'on') && (length(app.PatchHandle)>=ii)
                        cm=get(app.Figure,'Colormap');
                        cl=get(app.SliceAxis,'CLim');
                        pos=(ISOLevel(ii)-cl(1))/(cl(2)-cl(1))*length(cm);
                        newcolor=interp1(linspace(0,length(cm),length(cm)),cm,pos);
                    else
                        if ii==1
                            newcolor = [0 1 0];
                        else
                            newcolor = [0 0 1];
                        end
                    end
                    if length(app.PatchHandle)<ii || ~ishandle(app.PatchHandle(ii));% if patch doesn't exist, create it
                        app.PatchHandle(ii) = patch(surf1data,'FaceColor',newcolor,'EdgeColor','none','FaceAlpha',ISOAlpha(ii),'Parent',app.ISOAxis);
                        xl = get(app.ISOAxis,'XLabel');yl = get(app.ISOAxis,'YLabel');zl = get(app.ISOAxis,'ZLabel');
                        set(app.PatchHandle(ii),'UserData',{get(xl,'String');get(yl,'String');get(zl,'String')});
                    else % otherwise don't create new patch, just update underlying data
                        set(app.PatchHandle(ii), 'Faces', surf1data.faces, 'Vertices', surf1data.vertices,'FaceAlpha',ISOAlpha(ii),'FaceColor',newcolor);
                    end
                elseif (ISOAlpha(ii) <= 0 || ISOAlpha(ii) > 1 ) && ishandle(app.PatchHandle(ii))
                    % if the alpha is invalid AND the patch is currently
                    % rendered, delete the patch rather than setting it's
                    % Visible property to 'off' so animation speeds up
                    delete(app.PatchHandle(ii));
                end
            end
            % if in Manual recording mode, write to file
            if strcmp(app.RecordMode,'Manual')
                RecordToFile(app); % save screencapture to file
            end
        end
        
        function EnableDisableFig(app,EDF)
            % EnableDisableFig(app,EDF) is a helper function to enable or
            % disable the UI components in the app. The EDF flag is either
            % 'enable' or 'disable'.
            if strcmp(EDF,'enable')
                set([app.PlaySlider,app.GridComp{2,1},app.GridComp{1,2},...
                    app.GridComp{3,2},app.ISOSlider,app.DelayText,...
                    app.SampleNum,app.TimeNum,app.SliceValue,...
                    app.ISOLevelText,app.ISOAlphaText,app.PlayButton,...
                    app.DelayPlus,app.DelayMinus,app.OptionsMenu,...
                    app.RecordButton,app.AnimateButton],'Enable','on');
                if length(app.Data(1,1,1,:))==2 && all(all(all(app.Data(:,:,:,2)==0)))
                    % in this case there is no time dimension to animate
                    % on, so turn those components back off
                    set([app.DelayPlus,app.DelayMinus,app.DelayText,app.PlayButton,app.PlaySlider,...
                        app.SampleNum,app.TimeNum,app.OptionsMenu(10)],'Enable','off');
                end
            elseif strcmp(EDF,'disable')
                set([app.PlaySlider,app.GridComp{2,1},app.GridComp{1,2},...
                    app.GridComp{3,2},app.ISOSlider,app.DelayText,...
                    app.SampleNum,app.TimeNum,app.SliceValue,...
                    app.ISOLevelText,app.ISOAlphaText,app.PlayButton,...
                    app.DelayPlus,app.DelayMinus,app.OptionsMenu,app.AnimateButton,...
                    app.RecordButton],'Enable','off');
            else
                warning('Acceptable options for EnableDisableFig are ''enable'' and ''disable''');
                BeepOnWarn(app);
            end
        end
        
        function RecordToFile(app)
            % RecordToFile(app) is used to write a screenshot of the
            % desired app area to an open file.
            figure(app.Figure); % ensure figure is in front
            if get(app.RecordButton,'Value')==1
                % if recording, capture the figure and write to the file
                pos = get(app.Figure,'Position');
                switch app.RecordArea
                    case 'Whole App'
                        im = screencapture(app.Figure);
                    case 'Plot Panels Only'
                        im = screencapture(app.Figure,[10,160,pos(3)-5,pos(4)-160]);
                    case 'Slice Panel Only'
                        im = screencapture(app.Figure,[10,160,pos(3)/2-5,pos(4)-160]);
                    case 'ISO Panel Only'
                        im = screencapture(app.Figure,[pos(3)/2+7,160,pos(3)/2-5,pos(4)-160]);
                end
                writeVideo(app.WriteObject,im);
            end
        end
    end
    
    methods(Static)
        function output_txt = DataTipCB(~,event_obj)
            % output_txt = DataTipCB(obj,event_obj)
            % Custom data tip callback function
            % obj          Currently not used (empty)
            % event_obj    Handle to event object (i.e. slice)
            % output_txt   Data cursor text string (string or cell array of strings).
            pos = get(event_obj,'Position');
            tg = get(event_obj,'Target');
            idx = get(event_obj,'dataindex');
            names = get(tg,'UserData');
            SliceColor = get(tg,'CData'); % doesn't exist for patch
            output_txt = {[names{1},': ',num2str(pos(1),4)],...
                [names{2},': ',num2str(pos(2),4)],...
                [names{3},': ',num2str(pos(3),4)]};
            if ~isempty(SliceColor) % if it's not the patch
                output_txt{end+1} = ['Value: ',num2str(SliceColor(idx),4)];
            end
        end
    end
    
end