function varargout = CSDgui(varargin)
%CSDGUI M-file for CSDgui.fig
%      CSDGUI, by itself, creates a new CSDGUI or raises the existing
%      singleton*.
%
%      H = CSDGUI returns the handle to a new CSDGUI or the handle to
%      the existing singleton*.
%
%      CSDGUI('Property','Value',...) creates a new CSDGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CSDgui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CSDGUI('CALLBACK') and CSDGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CSDGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CSDgui

% Last Modified by GUIDE v2.5 27-Mar-2017 13:29:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CSDgui_OpeningFcn, ...
                   'gui_OutputFcn',  @CSDgui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before CSDgui is made visible.
function CSDgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for CSDgui
handles.output = hObject;
% Set file path
handles.filepath = 'Z:\EphysNew\';
handles.data.VEPdur = 0.3;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CSDgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CSDgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in SetPath.
function SetPath_Callback(hObject, ~, handles)
[filePath] = uigetdir('Z:\','Choose a folder');
handles.filePath = filePath;
set(handles.FilePath,'String',handles.filePath);
% Update handles structure
guidata(hObject,handles);
% --- 
function FilePath_CreateFcn(~, ~, ~)
function FilePath_Callback(hObject, ~, handles)
% Get file path
filePath = get(hObject,'String');
if ~strcmp('\',filePath(end))
    filePath = [filePath '\'];
end
% If filepath exists and is different from current file path
if exist(filePath,'dir') && ~strcmp(filePath,handles.filePath)
    % set new file path
    handles.filePath = filePath;
    % reset the gui
    guiReset(handles);
elseif ~exist(filePath,'dir')
    % If no such filepath send a warning
    warndlg('Can not find the chosen directory!','Invalid Path');
    set(handles.FilePath,'String',handles.filePath);
end
guidata(hObject,handles);

% --- 
function Animal_CreateFcn(~, ~, ~)
function Animal_Callback(hObject, eventdata, handles)
animalID = get(hObject,'String');
handles.animalID = animalID;
guidata(hObject,handles);
% --- 
function Experiment_CreateFcn(hObject, ~, ~)
function Experiment_Callback(hObject, eventdata, handles)
experimentID = get(hObject,'String');
handles.experimentID = experimentID;
guidata(hObject,handles);
% --- 
function VEPdur_CreateFcn(hObject, eventdata, handles)
function VEPdur_Callback(hObject, eventdata, handles)
VEPdur = get(hObject,'String');
handles.data.VEPdur = str2num(VEPdur);
guidata(hObject,handles);

% --- Executes on button press in LoadData.
function LoadData_Callback(hObject, ~, handles)
handles.NS3 = openNSx(strcat(handles.filepath,'data/',handles.animalID,'/',handles.experimentID,'/',strcat(handles.animalID,'_',handles.experimentID),'.ns3'));
handles.NEV = openNEV(strcat(handles.filepath,'data/',handles.animalID,'/',handles.experimentID,'/',strcat(handles.animalID,'_',handles.experimentID),'.nev'));
load([handles.filepath 'analyzer\',handles.animalID,'\',strcat(handles.animalID,'_',handles.experimentID),'.analyzer'],'-mat');
handles.Analyzer = Analyzer;
for j = 1:length(handles.Analyzer.P.param);
    handles.experimentParams.(eval(strcat('handles.Analyzer.P.param{',num2str(j),'}{1}'))) = ...
        eval(strcat('handles.Analyzer.P.param{',num2str(j),'}{3}'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.screenRate = 60; %Hz
        handles.Nelectrodes = 16;
        handles.elocs = [9 8 10 7 13 4 12 5 15 2 16 1 14 3 11 6];
        handles.elDist = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpEvents=de2bi(handles.NEV.Data.SerialDigitalIO.UnparsedData);
tmpEvents=double(tmpEvents);
stimOnI=find(diff(tmpEvents(:,4))==1 & tmpEvents(1:end-1,3)==1);
stimOn=handles.NEV.Data.SerialDigitalIO.TimeStampSec(stimOnI+1); 
Nstim=length(stimOn);
%compare number of events to analyzer file
try
if Nstim ~= getnotrials(handles.Analyzer)
    disp('not all events were detected!')
end
end
%get the reverals
reversal_rate = handles.experimentParams.t_period/handles.screenRate;
t_reversals = 0:reversal_rate:handles.experimentParams.stim_time;
t_reversals = t_reversals(2:end-1);
t_reversals = repmat(t_reversals,length(stimOn),1) + repmat(stimOn',1,length(t_reversals));
Nreversals = Nstim*size(t_reversals,2);
% compute VEP across all stimuli
% compute the number of samples per VEP
fsSlow=handles.NS3.MetaTags.SamplingFreq;
NVEP=fsSlow*handles.data.VEPdur; 
handles.data.NVEP = NVEP;
handles.data.fsSlow = fsSlow;
%filter the LFP - butterworth filter
LFPout = LFPfilt(double(handles.NS3.Data(1:handles.Nelectrodes,:)),0,fsSlow,50,1); 
LFPout = zscore(LFPout')';
LFPTimes = [0:1/fsSlow:(1/fsSlow)*(size(LFPout,2)-1)];
% get the spike-triggered-average from the VEPs
sta=0;
reversalIdx = zeros(1,Nreversals);
for k=1:Nreversals %for each stimulus presentation 
[~,reversalIdx(k)] = min(abs(LFPTimes - t_reversals(k)));
    pc = LFPout(:,reversalIdx(k):reversalIdx(k)+NVEP-1); %get the VEP 
    sta = sta+pc/Nreversals; 
end
% reorganize according to depth
staOrg = zeros(size(sta));
for n = 1:handles.Nelectrodes
   staOrg(n,:) = sta(handles.elocs(n),:); 
end
handles.data.staOrg = staOrg;
handles.data.smoothed = 0;
handles.data.normalized = 0;
guidata(hObject,handles);
disp('File loaded successfully')

% --- Executes on button press in Smoothed.
function Smoothed_Callback(hObject, ~, handles)
handles.data.smoothed = get(hObject,'Value');
guidata(hObject,handles);
    
% --- Executes on button press in Normalized.
function Normalized_Callback(hObject, ~, handles)
handles.data.normalized = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on selection change in PlotType.
function PlotType_CreateFcn(~, ~, ~)
function PlotType_Callback(hObject, ~, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns PlotType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotType
plots = cellstr(get(hObject,'String'));
plotType = plots{get(hObject,'Value')};
handles.plotType = plotType;
guidata(hObject,handles);

function CSDType_CreateFcn(hObject, eventdata, handles)
function CSDType_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
CSDType = contents{get(hObject,'Value')};
handles.CSDType = CSDType;
guidata(hObject,handles);

function update_Callback(hObject, eventdata, handles)
% hObject    handle to update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % compute the number of samples per VEP
    fsSlow=handles.NS3.MetaTags.SamplingFreq;
    NVEP=fsSlow*handles.data.VEPdur;
    % convert samples to ms
    timevec=(0:NVEP-1)/fsSlow*1000;
    staOrg = handles.data.staOrg;
    
    if handles.data.smoothed
        W = 1; %smooth time with a hann window of 1
        hh = [hann(W)' zeros(1,length(staOrg(1,:))-W)];
        hh = ones(length(staOrg(:,1)),1)*hh;
        staOrg = ifft(fft(staOrg,[],2).*abs(fft(hh,[],2)),[],2);
        staOrg = staOrg(1:1:end,:);
        W = 4; %smooth space with a hann window of 4
        hh = [hann(W); zeros(length(staOrg(:,1))-W,1)];
        hh = hh*ones(1,length(staOrg(1,:)));
        staOrg = ifft(fft(staOrg,[],1).*abs(fft(hh,[],1)),[],1);
    end
    
    if handles.data.normalized
        staOrg = zscore(staOrg');
        staOrg = staOrg';
    end
    
    if strcmp(handles.plotType,'LFP')
        axes(handles.axes2); cla; hold on; box on;
        sdom = (0:length(staOrg(:,1))-1)*handles.elDist;
        imagesc(timevec,sdom,staOrg), ylabel('depth (um)'); xlabel('time after stimulus');title('Stimulus triggered LFP')
        axis tight; set(gca,'YDir','reverse')
    elseif strcmp(handles.plotType,'CSD')
        if strcmp(handles.CSDType,'Convolution')
            k = normpdf([-1.5:0.24:1.5],0,1);
            csdfilt = [0.2.*k; 0.5.*k; k; 0.5.*k; 0.2.*k];
            handles.data.CSD = conv2(diff(conv2(diff(conv2(staOrg,csdfilt,'same'),1,1),csdfilt,'same'),1,1),csdfilt,'same');
            sdom = (0:length(staOrg(:,1))-1)*handles.elDist;
            handles.data.sdom = sdom(1):sdom(end)/(size(handles.data.CSD,1)-1):sdom(end);
            handles.data.CSDflag = 'conv';
        elseif strcmp(handles.CSDType,'Laplacian')
            handles.data.sdom = (0:length(staOrg(:,1))-1)*handles.elDist;
            handles.data.sdom = handles.data.sdom(2:end-1);
            [~, CSD] = gradient(staOrg(2:end-1,:));
            [~, handles.data.CSD] = gradient(CSD);
            handles.data.CSDflag = 'grad';
        end
        axes(handles.axes2); cla; hold on; box on; axis tight; set(gca,'YDir','reverse')
        imagesc(timevec,handles.data.sdom,handles.data.CSD), ylabel('depth (um)'); title('Current source density')
    elseif strcmp(handles.plotType,'VEPs')
        axes(handles.axes2); cla; hold on; box on; set(gca,'YDir','normal');
        for m=1:handles.Nelectrodes
        hold on
        plot(timevec,staOrg(m,:)*handles.elDist/2-handles.elDist*(m-1))
        end
        title('VEPs, organized by depth');
    end
    guidata(hObject,handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% support functions %

function [LFPout, H] = LFPfilt(LFP,lnflag,fs,hl,hh)

%This one has the line noise cancelation as optional (lnflag).  Also, you
%can put in 0 for hh if you don't wan't a high pass filter, and/or inf for
%hl if you don't want a low-pass filter.

%This one takes in the channels (i.e. not organized by space) and uses a
%band-pass Butterworth... It also uses a notch filter.

%hh - high-pass cutoff in Hz
%hl - low-pass cutoff in Hz

N = length(LFP(1,:));

%wfund = pepparam('refresh')/pepparam('h_period');  %Fundamental harmonic
wfund = 34/60;

%wfund = 120.3115/4;  %for most cat randori

w = linspace(-pi,pi,N);
ws = 2*pi*wfund*(1:4)/fs;  %This gets rid of harmonics as well
%ws = 2*pi*wfund/fs;
ws = [-fliplr(ws) ws];
%%%
%lnflag = 1;
%%%
if lnflag == 1
    rs = ones(1,length(ws));
    comb = ones(1,N);
    for i = 1:length(ws)
        comb = comb.*(exp(1i*w) - rs(i)*exp(1i*ws(i)));
    end
    rs = .9*ones(1,length(ws));
    for i = 1:length(ws)
        comb = comb./(exp(1i*w) - rs(i)*exp(1i*ws(i)));
    end
    comb = comb-min(comb);
end

ordL = 8;
ordH = 4;
wl = 2*pi*hl/fs;
wh = 2*pi*hh/fs;
Butter_LP = 1./(1+(1i*w/wl).^(2*ordL));
Butter_HP = 1./(1+(wh./(1i*w)).^(2*ordH));


if hh ~= 0 & hl == inf
    Butter = Butter_HP;
elseif hh == 0 & hl ~= inf
    Butter = Butter_LP;
elseif hh ~= 0 & hl ~= inf
    Butter = Butter_LP.*Butter_HP;
else 
    Butter = ones(size(Butter_LP));
end


if lnflag == 1
    H = Butter.*abs(comb);
else
    H = Butter;
end

H = fftshift(H);  %make it go 0 to 2pi

LFPout = NaN*zeros(size(LFP));


for i = 1:length(LFP(:,1))
    s = LFP(i,:);
    s = real(ifft(abs(H).*fft(s)));
    LFPout(i,:) = s;
end
