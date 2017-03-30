function varargout = EphysSpikes(varargin)
% EPHYSSPIKES MATLAB code for EphysSpikes.fig
%      EPHYSSPIKES, by itself, creates a new EPHYSSPIKES or raises the existing
%      singleton*.
%
%      H = EPHYSSPIKES returns the handle to a new EPHYSSPIKES or the handle to
%      the existing singleton*.
%
%      EPHYSSPIKES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EPHYSSPIKES.M with the given input arguments.
%
%      EPHYSSPIKES('Property','Value',...) creates a new EPHYSSPIKES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EphysSpikes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EphysSpikes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EphysSpikes

% Last Modified by GUIDE v2.5 21-Feb-2016 15:10:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EphysSpikes_OpeningFcn, ...
                   'gui_OutputFcn',  @EphysSpikes_OutputFcn, ...
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


% --- Executes just before EphysSpikes is made visible.
function EphysSpikes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EphysSpikes (see VARARGIN)

% Choose default command line output for EphysSpikes
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EphysSpikes wait for user response (see UIRESUME)
% uiwait(handles.figure1);
if exist(getSettingsPath,'file')
    load(getSettingsPath);
    set(handles.textStatus,'string','Saved settings loaded');
else
    user = lower(input('Cannot find settings file. What is your name?: ','s'));
    useSynologyRaw = input('Do you want to use the Synology drive for getting raw data? (1=yes, 0=no): ');
    disp('Creating required folder structures and setting files. Please check, reconfirm and manually reconfigure folders by editing the settings.mat file. Also, please create outShare folder structure manually.');
    settings = defaultSettings(user,useSynologyRaw);
    
    if ~useSynologyRaw
        setupFolderStructure('',1,0);
    end
    setupFolderStructure('',0,1);
    save(getSettingsPath,'settings');
    set(handles.textStatus,'string','Default settings loaded');
end

set(handles.analyzerPath,'string',settings.rawAnalyzerPath);
set(handles.path,'string',settings.rawDataPath);
set(handles.animal,'string',settings.animal);
set(handles.unit,'string',settings.unit);
set(handles.experiment,'string',settings.experiment);


% --- Outputs from this function are returned to the command line.
function varargout = EphysSpikes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function animal_Callback(hObject, eventdata, handles)
% hObject    handle to animal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of animal as text
%        str2double(get(hObject,'String')) returns contents of animal as a double
load(getSettingsPath);
settings.animal = get(hObject,'String');
save(getSettingsPath,'settings');


% --- Executes during object creation, after setting all properties.
function animal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to animal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unit_Callback(hObject, eventdata, handles)
% hObject    handle to unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unit as text
%        str2double(get(hObject,'String')) returns contents of unit as a double
load(getSettingsPath);
settings.unit = get(hObject,'String');
save(getSettingsPath,'settings');

% --- Executes during object creation, after setting all properties.
function unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function experiment_Callback(hObject, eventdata, handles)
% hObject    handle to experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of experiment as text
%        str2double(get(hObject,'String')) returns contents of experiment as a double
load(getSettingsPath);
settings.experiment = get(hObject,'String');
save(getSettingsPath,'settings');

% --- Executes during object creation, after setting all properties.
function experiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonLoadData.
function buttonLoadData_Callback(hObject, eventdata, handles)
clear Thresholds d ChId NS3 NS6 NEV
global NS6 Thresholds d ChID Analyzer NS3 NEV

load(getSettingsPath);
path = get(handles.path,'string');
analyzerpath = get(handles.analyzerPath,'string');
animal = lower(lower(get(handles.animal,'string')));
unit = get(handles.unit,'string');
experiment = get(handles.experiment,'string');

fullFileName = [animal '_' unit '_' experiment];

analyzerFullPath = [analyzerpath settings.filepathSlash animal settings.filepathSlash fullFileName '.analyzer'];
rawDataFullPath = [path settings.filepathSlash animal settings.filepathSlash unit '_' experiment settings.filepathSlash fullFileName];

load(analyzerFullPath,'-mat');
openNSx([rawDataFullPath '.ns3']);
openNSx([rawDataFullPath '.ns6']);
NEV = openNEV([settings.nevFilePrefix rawDataFullPath '.nev']);

% load(['Z:\Ephys\AnalyzerFiles\',get(handles.animal,'string'), '\', get(handles.animal,'string'),'_',get(handles.unit,'string'),'_',get(handles.experiment,'string'),'.analyzer'],'-mat');
% openNSx(['Z:\Ephys\' get(handles.animal,'string') '\' get(handles.unit,'string') '_' get(handles.experiment,'string') '\' get(handles.animal,'string') '_' get(handles.unit,'string') '_' get(handles.experiment,'string') '.ns6']);
% openNSx(['Z:\Ephys\' get(handles.animal,'string') '\' get(handles.unit,'string') '_' get(handles.experiment,'string') '\' get(handles.animal,'string') '_' get(handles.unit,'string') '_' get(handles.experiment,'string') '.ns3']);
% NEV = openNEV(['Z:\Ephys\' get(handles.animal,'string') '\' get(handles.unit,'string') '_' get(handles.experiment,'string') '\' get(handles.animal,'string') '_' get(handles.unit,'string') '_' get(handles.experiment,'string') '.NEV']);

disp('Filtering Data'); set(handles.textStatus,'string','Filtering Data'); drawnow
[b,a] = butter(3,[250 5000]/15000,'bandpass');
if length(NS6.Data) == 2
NS6.Data = [NS6.Data{2}];
end
NS6.Data = double(NS6.Data);
Thresholds = zeros(length(NS6.Data(:,1)),1);
for i = 1:length(NS6.Data(:,1));
    disp(['Filtering Ch ' num2str(i)]);
    set(handles.textStatus,'string',['Filtering Ch ' num2str(i)]); drawnow
    NS6.Data(i,:) = filter(b,a,NS6.Data(i,:));
    ChStr{i} = i;
    Mapping{i,1} = NS6.ElectrodesInfo(1,i).ElectrodeID;
    ChID(i) = NS6.ElectrodesInfo(1,i).ElectrodeID;
    Mapping{i,2} = 1;
    Thresholds(i) = mean(NS6.Data(i,:)) - std(NS6.Data(i,:))*3;
end
disp('Done filtering'); set(handles.textStatus,'string','Done filtering'); drawnow
set(handles.Chs,'String',ChStr);
set(handles.Mapping,'Data',Mapping);
d = figure;
figure(d)
subplot(3,1,1)
hold off
plot(NS6.Data(1,1:360000));
hold on
title('First minute')
plot([1 360000],[Thresholds(1) Thresholds(1)], 'color', [1 0 0])
subplot(3,1,2)
Mddle = length(NS6.Data(1,:))/2;
hold off
plot(NS6.Data(1,Mddle:Mddle+360000-1));
hold on
title('Middle minute')
plot([1 360000],[Thresholds(1) Thresholds(1)], 'color', [1 0 0])
subplot(3,1,3)
hold off
plot(NS6.Data(1,end-360000+1:end));
hold on
title('Last minute')
plot([1 360000],[Thresholds(1) Thresholds(1)], 'color', [1 0 0])

% hObject    handle to buttonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit4_Callback(hObject, eventdata, handles)
global d Thresholds NS6
channel = get(handles.Chs,'Value');
Thresholds(channel) = eval(get(handles.edit4,'String'));
figure(d)
subplot(3,1,1)
hold off
plot(NS6.Data(channel,1:360000));
hold on
title('First minute')
plot([1 360000],[Thresholds(channel) Thresholds(channel)], 'color', [1 0 0])
subplot(3,1,2)
Mddle = length(NS6.Data(1,:))/2;
hold off
plot(NS6.Data(channel,Mddle:Mddle+360000-1));
hold on
title('Middle minute')
plot([1 360000],[Thresholds(channel) Thresholds(channel)], 'color', [1 0 0])
subplot(3,1,3)
hold off
plot(NS6.Data(channel,end-360000+1:end));
hold on
title('Last minute')
plot([1 360000],[Thresholds(channel) Thresholds(channel)], 'color', [1 0 0])


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonSaveData.
function buttonSaveData_Callback(hObject, eventdata, handles)
    global Thresholds NS6 ChID Analyzer NS3 NEV Pulses
    Times = [];

    %Find spikes and waveforms
    Mapping = get(handles.Mapping,'Data');
    for i = 1:length(Mapping(1,:))
        for j = 1:length(Mapping(:,1))
            Mappingtmp(i,j) = double(Mapping{j,i});
        end
    end
    
    Mapping = Mappingtmp;
    sites = unique(Mapping(2,:));
    for site = 1:length(sites)
        channels = Mapping(1,(Mapping(2,:) == site));
        for Ch = 1:length(channels)
            channel = find(ChID == channels(Ch));
            channels(Ch) = find(ChID == channels(Ch));
            disp(['Extracting channel ' int2str(channel)]);
            set(handles.textStatus,'string',['Extracting channel ' int2str(channel)]); drawnow
            BelowTh = NS6.Data(channel,:) < Thresholds(channel);
            BelowTh = BelowTh(2:end) - BelowTh(1:end-1);
            BelowTh(1:24) = 0;
            BelowTh(end-24:end) = 0;
            LT = length(Times);
            if Thresholds(channel) > 0
                Times = [Times find(BelowTh == -1)];
            else
                Times = [Times find(BelowTh == 1)];
            end
            i = LT+2;
            while i <= length(Times)
                if Times(i) - Times(i-1) < 60
                    Times = [Times(1:i-1) Times(i+1:end)];
                end
                i = i+1;
            end
        end
        Times = sort(Times);
        dT = Times(2:end)-Times(1:end-1);
        I = dT<25;
        Times(I)=[];
        I = Times<25;
        Times(I)=[];
        I = Times>length(NS6.Data(1,:))-40;
        Times(I)=[];
        Spikes{site}.TimeStamp = Times';
        Spikes{site}.Unit = ones(length(Times),1);
        disp(['Getting wvfs of site ' int2str(site)]);
        set(handles.textStatus,'string',['Getting wvfs of site ' int2str(site)]); drawnow
        fprintf('\n');
        for i = 1:length(Times)
            Spikes{site}.Waveform(i,:,:) = double(NS6.Data(channels,Times(i)-14:Times(i)+25))-repmat(mean(NS6.Data(channels,Times(i)-14:Times(i)+25),2),1,40);
            if mod(i,floor(length(Times)/100)) == 0
                fprintf('.');
                set(handles.textStatus,'string',['Getting wvfs of site ' int2str(site) ': ' num2str(floor(i/floor(length(Times)/100))) '% complete']); drawnow
            end
        end
    end
    for i = 1:length(Spikes)
        UnitType{i} = [3];
    end

    %find trial events
    Events.Type{1} = [];
    if length(NS3.Data(:)) == 2
        NS3.Data = NS3.Data{2};
    end
    
    samplingFreq = 30000;
    analogChSamplingFreq = 2000;
    StimT = getparam('stim_time',Analyzer)*samplingFreq; % Analyzer.P.param{1,3}{1,3}*30000;
    PostD = getparam('postdelay',Analyzer)*samplingFreq; % Analyzer.P.param{1,2}{1,3}*30000;
    PreD = getparam('predelay',Analyzer)*samplingFreq; % Analyzer.P.param{1,1}{1,3}*30000;

    Reps = getnorepeats(1,Analyzer);
    [~,domval,blankid]=getdomainvalue(Analyzer);
    try
    BReps = getnorepeats(blankid,Analyzer);
    end
    trialList = getcondtrial(Analyzer);
    
    Times = zeros(length(trialList)*4,1);
    
    %with pulses
    if get(handles.Pulses,'Value')
        % use the analyzer utils apis instead of doing this yourself.
        % already done outside the if block
%         Reps = length(Analyzer.loops.conds{1,1}.repeats);
%         BReps =length(Analyzer.loops.conds{1,end}.repeats);
%         Reps = getnorepeats(1,Analyzer);
%         [~,domval,blankid]=getdomainvalue(Analyzer);
%         BReps = getnorepeats(blankid,Analyzer);
        
%         Times = zeros(size(domval,1)*Reps*4 + BReps*4,1);
        Pulses = [];
        PulseChs = eval(get(handles.PulsesCh,'String'));
        for j = 1:length(PulseChs)
            CH = 1;
            while NS3.ElectrodesInfo(1,CH).ElectrodeID ~= PulseChs(j)
                CH = CH+1;
                if CH > length(NS3.ElectrodesInfo)
                    return
                end
            end
            AinpData = abs(NS3.Data(CH,2:end)-NS3.Data(CH,1:end-1));
            Thres = max(AinpData)/5;
            Pulses = [Pulses find(AinpData > Thres)];
        end
        
        Pulses = sort(Pulses);
        Pulses = Pulses*15;

        % Identifing real pulses based on Delay expectation
%         StimT = Analyzer.P.param{1,3}{1,3}*30000;
%         PostD = Analyzer.P.param{1,2}{1,3}*30000;
%         PreD = Analyzer.P.param{1,1}{1,3}*30000;
        ind = 1;
        PInd = 0;
        Pulses = [Pulses Pulses(end)+PostD+analogChSamplingFreq];
        while ind<length(Times)
            error = StimT+PostD+PreD;
            while error >(StimT+PostD+PreD)*0.1
                PInd = PInd +1;
                error = 0;
                TriBeg = Pulses(PInd);
                StimBeg = TriBeg+PreD;
                StimEnd = StimBeg + StimT;
                TriEnd = StimEnd+PostD;
                SB = 1;
                while StimBeg-Pulses(PInd + SB)>0
                    SB = SB+1;
                end
                if abs(StimBeg-Pulses(PInd + SB)) > abs(StimBeg-Pulses(PInd + SB-1))
                    SB = SB-1;
                end
                error = error + abs(StimBeg-Pulses(PInd + SB));
                SE = SB+1;
                while StimEnd-Pulses(PInd + SE)>0
                    SE = SE+1;
                end
                if abs(StimEnd-Pulses(PInd + SE)) > abs(StimEnd-Pulses(PInd + SE-1))
                    SE = SE-1;
                end
                error = error + abs(StimEnd-Pulses(PInd + SE));
                    TE = SE+1;
                while TriEnd-Pulses(PInd + TE)>0
                    TE = TE+1;
                end
                if abs(TriEnd-Pulses(PInd + TE)) > abs(TriEnd-Pulses(PInd + TE-1))
                    TE = TE-1;
                end
                error = error + abs(TriEnd-Pulses(PInd + TE));
                if error >40000
                    disp('Following are the error events based on pulses: ');
                    disp(error)
                    disp(TE)
                    disp(SE)
                    disp(SB)
                    disp(TriBeg)
                end
            end
            Times(ind:ind+3) = [Pulses(PInd) Pulses(PInd+SB) Pulses(PInd+SE) Pulses(PInd+TE)];
            PInd = PInd+TE;
            ind = ind+4;
        end
        Events.Type{1} = 'Pulses';
        Times = reshape(Times,4,length(Times)/4);
        Times = Times';
        Events.Timestamp{1} = Times;
    end

    % with photodiode
    if get(handles.Photodiode,'Value')
        PulseCh = eval(get(handles.PhotodCh,'String'));
        AinpData = (NS3.Data(PulseCh,:));
        p = zeros(length(AinpData),1);
        for ii=1:length(AinpData)-1
            if (AinpData(ii) <= 100 && AinpData(ii+1) > 100) || ...
                (AinpData(ii) >= 100 && AinpData(ii+1) < 100)
                p(ii) = 1;
            end
        end
        Pulses = find(p);
%         Pulses = find(AinpData >= 100);
        Pulses = sort(Pulses);
        Postd = getparam('postdelay',Analyzer)*analogChSamplingFreq;
        Pulses = [Pulses;Pulses(end)+Postd];
        Pulses = Pulses*15;

        % Identifing real pulses based on Delay expectation
        ind = 1;
        PInd = 0;
        Pulses = [Pulses;Pulses(end)+PostD+analogChSamplingFreq];
        while ind<length(Times)
            error = StimT+PostD+PreD;
            while error >(StimT+PostD+PreD)*0.1
                PInd = PInd +1;
                error = 0;
                TriBeg = Pulses(PInd);
                StimBeg = TriBeg+PreD;
                StimEnd = StimBeg + StimT;
                TriEnd = StimEnd+PostD;
                SB = 1;
                while StimBeg-Pulses(PInd + SB)>0
                    SB = SB+1;
                end
                if abs(StimBeg-Pulses(PInd + SB)) > abs(StimBeg-Pulses(PInd + SB-1))
                    SB = SB-1;
                end
                error = error + abs(StimBeg-Pulses(PInd + SB));
                SE = SB+1;
                while StimEnd-Pulses(PInd + SE)>0
                    SE = SE+1;
                end
                if abs(StimEnd-Pulses(PInd + SE)) > abs(StimEnd-Pulses(PInd + SE-1))
                    SE = SE-1;
                end
                error = error + abs(StimEnd-Pulses(PInd + SE));
                    TE = SE+1;
                while TriEnd-Pulses(PInd + TE)>0
                    TE = TE+1;
                end
                if abs(TriEnd-Pulses(PInd + TE)) > abs(TriEnd-Pulses(PInd + TE-1))
                    TE = TE-1;
                end
                error = error + abs(TriEnd-Pulses(PInd + TE));
                if error > 40000
                    disp('Following are the error events based on photodiode: ');
                    disp(error)
                    disp(TE)
                    disp(SE)
                    disp(SB)
                    disp(TriBeg)
                end
            end
            Times(ind:ind+3) = [Pulses(PInd) Pulses(PInd+SB) Pulses(PInd+SE) Pulses(PInd+TE)];
            PInd = PInd+TE;
            ind = ind+4;
        end
        if isempty(Events.Type{length(Events.Type)})
            Events.Type{length(Events.Type)} = 'Photodiode';
        else
            Events.Type{length(Events.Type)+1} = 'Photodiode';
        end
        Times = reshape(Times,4,length(Times)/4);
        Times = Times';
        Events.Timestamp{length(Events.Type)} = Times;
    end
    
    % with digital input
    Times = zeros(length(trialList)*4,1);
    Pulses = double(NEV.Data.SerialDigitalIO.TimeStamp);
    Pulses = sort(Pulses);
    UnpDta = abs(diff(double(NEV.Data.SerialDigitalIO.UnparsedData')));
    if ~isempty(Pulses)
        % already done outside the if block
%         StimT = Analyzer.P.param{1,3}{1,3}*30000;
%         PostD = Analyzer.P.param{1,2}{1,3}*30000;
%         PreD = Analyzer.P.param{1,1}{1,3}*30000;

        % use the analyzer utils apis instead of doing this yourself.
        % already done outside the if block
%         Reps = length(Analyzer.loops.conds{1,1}.repeats);
%         BReps =length(Analyzer.loops.conds{1,end}.repeats);
%         Reps = getnorepeats(1,Analyzer);
%         [~,~,blankid]=getdomainvalue(Analyzer);
%         BReps = getnorepeats(blankid,Analyzer);

        dPulses = diff(Pulses);
        Pulses([find(dPulses<2) find(dPulses<5)+1]) = [];
        UnpDta([find(dPulses<2) find(dPulses<5)-1]) = [];
        Pulses(find((UnpDta ~= 8)&(UnpDta ~= 4)&(UnpDta ~= 24))+1) = [];
%         Times = zeros(((length(Analyzer.loops.conds)-1)*Reps*4+BReps*4),1);
        if length(Pulses) ~= length(Times)
            warndlg('Digital pulses are more than 4 X NunberOfTrials using expected times')
            ind = 1;
            PInd = 0;
            Pulses = [Pulses Pulses(end)+PostD Pulses(end)+PostD+analogChSamplingFreq];
            while ind<length(Times)
                error = StimT+PostD+PreD;
                while error >(StimT+PostD+PreD)*0.1
                    PInd = PInd +1;
                    error = 0;
                    TriBeg = Pulses(PInd);
                    StimBeg = TriBeg+PreD;
                    StimEnd = StimBeg + StimT;
                    TriEnd = StimEnd+PostD;
                    SB = 1;
                    while StimBeg-Pulses(PInd + SB)>0
                        SB = SB+1;
                    end
                    if abs(StimBeg-Pulses(PInd + SB)) > abs(StimBeg-Pulses(PInd + SB-1))
                        SB = SB-1;
                    end
                    error = error + abs(StimBeg-Pulses(PInd + SB));
                    SE = SB+1;
                    while StimEnd-Pulses(PInd + SE)>0
                        SE = SE+1;
                    end
                    if abs(StimEnd-Pulses(PInd + SE)) > abs(StimEnd-Pulses(PInd + SE-1))
                        SE = SE-1;
                    end
                    error = error + abs(StimEnd-Pulses(PInd + SE));
                        TE = SE+1;
                    while TriEnd-Pulses(PInd + TE)>0
                        TE = TE+1;
                    end
                    if abs(TriEnd-Pulses(PInd + TE)) > abs(TriEnd-Pulses(PInd + TE-1))
                        TE = TE-1;
                    end
                    error = error + abs(TriEnd-Pulses(PInd + TE));
                    if error >(StimT+PostD+PreD)*0.1
                        disp('Following are the error events based on digital pulses: ');
                        disp(error)
                        disp(TE)
                        disp(SE)
                        disp(SB)
                        disp(TriBeg)
                    end
                end
                Times(ind:ind+3) = [Pulses(PInd) Pulses(PInd+SB) Pulses(PInd+SE) Pulses(PInd+TE)];
                PInd = PInd+TE;
                ind = ind+4;
            end
        else
            disp('Using straight Digi pulses')
            Times = Pulses;
        end
        
        Times = reshape(Times,4,length(Times)/4);
        Times = Times';
        if isempty(Events.Type{length(Events.Type)})
        Events.Type{length(Events.Type)} = 'Digital';
        else
        Events.Type{length(Events.Type)+1} = 'Digital';
        end
        Events.Timestamp{length(Events.Type)} = Times;
    end

    %save
    load(getSettingsPath);
    save([settings.outSpikeFilePath settings.filepathSlash lower(get(handles.animal,'string')) '_' get(handles.unit,'String') '_' get(handles.experiment,'String') '_spikes'], 'Spikes','UnitType','Mapping','Events');
    save([settings.outShareSpikeFilePath settings.filepathSlash lower(get(handles.animal,'string')) '_' get(handles.unit,'String') '_' get(handles.experiment,'String') '_spikes'], 'Spikes','UnitType','Mapping','Events');
    
    choice = questdlg('Do you want to save this record in the summary file?', ...
        'Summary file save', ...
        'Yes','No','Yes');

    switch choice
        case 'Yes'
            set(handles.textStatus,'string','Saving record in summary file. Please go to command window to enter details.'); drawnow
        case 'No'
            set(handles.textStatus,'string','Finished processing. Spike file saved. Proceed to sorting.'); drawnow
            return;
    end
    
    % open table data
    [num,~,raw] = xlsread([settings.outSummaryFilePath settings.filepathSlash settings.summaryFileName],1);
    AnimalAge = raw([0;num(:,6)] ==  1,[1 2]);
    [~,id,~] = unique(AnimalAge(:,1));
    AnimalAge = AnimalAge(id,:);
    ID = 0;
    for i = 1:length(AnimalAge(:,1))
        if isequal(AnimalAge{i,1},lower(get(handles.animal,'string')))
            ID = i;
        end
    end

    %Ask for area and age and experiment
    Area = input('Please provide recording area: ','s');
    Experiment = input('Please provide experiment class: ','s');
    if ID == 0
        Age = input('Please provide animal Age: ');
    else
        Age = AnimalAge(ID,2);
    end
    
    
    %mod table
    Rows = [];
    for i = 1:length(raw(:,1))
        if isequal(raw(i,3),{[lower(get(handles.animal,'string')) '_' get(handles.unit,'String') '_' get(handles.experiment,'String')]})
            Rows = [Rows i];
        end
    end
    raw(Rows,:) = [];
    raw(end+1,1) = {lower(get(handles.animal,'string'))};
    raw(end,2) = {Age};
    raw(end,3) = {[lower(get(handles.animal,'string')) '_' get(handles.unit,'String') '_' get(handles.experiment,'String')]};
    raw(end,5) = {Experiment};
    raw(end,9) = {'Not sorted'};
    raw(end,18) = {Area};
    xlswrite([settings.outSummaryFilePath settings.filepathSlash settings.summaryFileName],raw,1)
    xlswrite([settings.outShareSummaryFilePath settings.filepathSlash settings.summaryFileName],raw,1)
    set(handles.textStatus,'string','Finished processing. Spike/summary file saved. Proceed to sorting.');

% --- Executes on selection change in Chs.
function Chs_Callback(hObject, eventdata, handles)
global d Thresholds NS6
channel = get(handles.Chs,'Value');
figure(d)
subplot(3,1,1)
hold off
plot(NS6.Data(channel,1:360000));
hold on
title('First minute')
plot([1 360000],[Thresholds(channel) Thresholds(channel)], 'color', [1 0 0])
subplot(3,1,2)
Mddle = length(NS6.Data(1,:))/2;
hold off
plot(NS6.Data(channel,Mddle:Mddle+360000-1));
hold on
title('Middle minute')
plot([1 360000],[Thresholds(channel) Thresholds(channel)], 'color', [1 0 0])
subplot(3,1,3)
hold off
plot(NS6.Data(channel,end-360000+1:end));
hold on
title('Last minute')
plot([1 360000],[Thresholds(channel) Thresholds(channel)], 'color', [1 0 0])
% hObject    handle to Chs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Chs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Chs


% --- Executes during object creation, after setting all properties.
function Chs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Chs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TetrodeBool.
function TetrodeBool_Callback(hObject, eventdata, handles)
% hObject    handle to TetrodeBool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TetrodeBool


% --- Executes on button press in Digital.
function Digital_Callback(hObject, eventdata, handles)
% hObject    handle to Digital (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Digital


% --- Executes on button press in AllSngle.
function AllSngle_Callback(hObject, eventdata, handles)
Mapping = get(handles.Mapping,'Data');
for i = 1:length(Mapping(:,1))
Mapping{i,2} = i;
end
set(handles.Mapping,'Data',Mapping);
% hObject    handle to AllSngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AllTtd.
function AllTtd_Callback(hObject, eventdata, handles)
Mapping = get(handles.Mapping,'Data');
for i = 1:length(Mapping(:,1))
Mapping{i,2} = ceil(i/4);
end
set(handles.Mapping,'Data',Mapping);
% hObject    handle to AllTtd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Photodiode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Photodiode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function PhotodCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhotodCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over buttonLoadData.
function buttonLoadData_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to buttonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Pulses.
function Pulses_Callback(hObject, eventdata, handles)
% hObject    handle to Pulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Pulses


% --- Executes during object creation, after setting all properties.
function Pulses_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function PulsesCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PulsesCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function path_Callback(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path as text
%        str2double(get(hObject,'String')) returns contents of path as a double
load(getSettingsPath);
settings.path = get(hObject,'String');
save(getSettingsPath,'settings');

% --- Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function analyzerPath_Callback(hObject, eventdata, handles)
% hObject    handle to analyzerPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of analyzerPath as text
%        str2double(get(hObject,'String')) returns contents of analyzerPath as a double
load(getSettingsPath);
settings.analyzerPath = get(hObject,'String');
save(getSettingsPath,'settings');

% --- Executes during object creation, after setting all properties.
function analyzerPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analyzerPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function settingsPath = getSettingsPath
if ispc
    settingsPath = 'settings\currSettings.mat';
else
    settingsPath = 'settings/currSettings.mat';
end
