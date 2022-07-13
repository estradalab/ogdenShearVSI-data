%% Measure analog input channels on the Measurement(DT987x or DT887x) instrument
%   
%  Copyright (C) 2009 DataTranslation Inc.


%% Open connection to instrument.
% RsrcName is an IVI logical name or an instrument specific string that 
% identifies the address of the instrument, such as a VISA resource descriptor
% string. For USB Measurement instruments, specify USB::InstrumentName,
% where InstrumentName is the name of the Measurement instrument in the 
% Open Layers Control Panel. USB::DT987x(00) is an example of a ResourceName
% for the DT987x.

RsrcName = 'TCPIP::192.43.218.113::SOCKETS';
dev = icdevice('DtxMeasurement_DtxMeasurement.mdd', RsrcName);

try
    connect(dev);
    
    % Enable protected commands to function
    comobj = get(dev,'System');
    invoke(comobj, 'EnableProtectedCommands', 'admin');
    
    % Get the instrument identity
    comobj = get(dev, 'Identity');
    propertyValue = get(comobj, 'InstrumentModel');
    str = strcat ('InstrumentModel= ',propertyValue);  
    disp(str);
    propertyValue = get(comobj, 'InstrumentManufacturer');
    str = strcat ('InstrumentManufacturer= ',propertyValue);  
    disp(str);
    propertyValue = get(comobj, 'InstrumentFirmwareRevision');
    str = strcat ('InstrumentFirmwareRevision= ',propertyValue);  
    disp(str);
    propertyValue = get(comobj, 'Description');
    str = strcat ('Description= ',propertyValue);  
    disp(str);
    propertyValue = get(comobj, 'Identifier');
    str = strcat ('Identifier= ',propertyValue);  
    disp(str);
    propertyValue = get(comobj, 'Vendor');
    str = strcat ('Vendor= ',propertyValue);  
    disp(str);
    propertyValue = get(comobj, 'Revision');
    str = strcat ('Revision= ',propertyValue);  
    disp(str);
    
    % If the first channel type is thermocouple, set its thermocouple type
    % to K
    chanType = get(dev.Channel(1), 'ChannelType');
    if strcmp(chanType,'DtxMeasurementChannelTypeThermocouple')
        set(dev.Channel(1), 'ThermocoupleType', 'DtxMeasurementThermocoupleTypeK');
    end
    
    % If the first channel type is RTD, set its RTD type to PT100
    if strcmp(chanType,'DtxMeasurementChannelTypeRtd')
        set(dev.Channel(1), 'RtdType', 'DtxMeasurementRtdTypeAmericanPT100');
    end
    
    % If the first channel type is MultiRange, set the RangeType type to
    % Bipolar 10 Volts
    if strcmp(chanType,'DtxMeasurementChannelTypeMultiRange')
        set(dev.Channel(1), 'RangeType', 'DtxMeasurementRangeTypeBip10Volts');
    end
        
    % Enable channel 0,1 and 2 for scanning
    comobj = get(dev, 'Channels');
    comobj = comobj(1);
    channels = int32([0;1;2]);
    invoke(comobj, 'Configure', channels);
  
    % Specify the rate at which to scan the list of enabled channels (in seconds).
    comobj = get(dev, 'Trigger');
    set(comobj, 'TimerInterval', .2);

    % Initiate the measurements
    comobj = get(dev, 'Acquisition');
    invoke(comobj, 'Initiate');
    pause(2);
    isRunning = 1;
    
    % Read 10 scans every 2 seconds for 1 minute
    RequestedScansToRead = 10;
    RequestedScansIndex = 0;
    
    for timeInc = 1:30 
        [ScansIndex, State] = invoke(comobj, 'GetStatus', 0, 0);
  
        if (State ==1)
            [ActualScansIndex, ActualScansRead, StartTimeInSeconds, StartTimeInMilliSeconds, Samples] = invoke(comobj, 'Fetch', int32(RequestedScansIndex), int32(RequestedScansToRead), int32(0), int32(0), int32([0;0]), int32([0;0]), double([0;0]));
            RequestedScansIndex = ActualScansIndex+ActualScansRead;
            disp(['Actual Scans Index: ',num2str(ActualScansIndex) ,'  Actual Scans Read: ', num2str(ActualScansRead)]);
            disp(Samples);
            pause(2);
        else
            break;
        end
    end
    
    invoke(comobj, 'Abort');
    
    % Disable protected commands 
    comobj = get(dev,'System');
    invoke(comobj, 'DisableProtectedCommands', 'admin');
 
catch DtxMeasurementerror
    disp(['Error id: ', DtxMeasurementerror.identifier]);
    disp(['Error Message: ',DtxMeasurementerror.message]);
end

% Close the connection with instrument
disconnect(dev);
delete(dev);   
    