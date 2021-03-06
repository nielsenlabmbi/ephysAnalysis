function varargout = Sorter(varargin)
% SORTER MATLAB code for Sorter.fig
%      SORTER, by itself, creates a new SORTER or raises the existing
%      singleton*.
%
%      H = SORTER returns the handle to a new SORTER or the handle to
%      the existing singleton*.
%
%      SORTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SORTER.M with the given input arguments.
%
%      SORTER('Property','Value',...) creates a new SORTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sorter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sorter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Sorter

% Last Modified by GUIDE v2.5 06-Oct-2015 12:14:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Sorter_OpeningFcn, ...
                   'gui_OutputFcn',  @Sorter_OutputFcn, ...
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


% --- Executes just before Sorter is made visible.
function Sorter_OpeningFcn(hObject, eventdata, handles, varargin)
global Str UnitIndx PCs PC TetrodeBool
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Sorter (see VARARGIN)

% Choose default command line output for Sorter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes Sorter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Sorter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ChannelMenu.
function ChannelMenu_Callback(hObject, eventdata, handles)
global Spikes UnitType PC
site = get(handles.SiteMenu,'Value')-1;
unit = get(handles.UnitMenu,'Value')-1;
In = get(handles.ChannelMenu,'Value')-1;
if In == 0
    disp(In)
    return 
end
axes(handles.axes1)
hold off
if get(handles.allUnits,'Value')
Spks = find(Spikes{site}.Unit>0);
UnitsPloted = unique(Spikes{site}.Unit(Spks));
else
if get(handles.CompPair,'Value')
unit2 = get(handles.UCompare,'Value')-1;
if unit2 >= unit
    unit2 = unit2 + 1;
end
Spks = find(Spikes{site}.Unit==unit);
Spks2 = find(Spikes{site}.Unit==unit2);
else
Spks = find(Spikes{site}.Unit==unit);
end
end
if get(handles.allUnits,'Value')
    for i = 1:length(UnitsPloted)
    Spks = find(Spikes{site}.Unit==UnitsPloted(i));
    Tmp = length(UnitsPloted)-1;
    if Tmp == 0
        Tmp = 1;
    end
    color = [normpdf((i-1)/Tmp,0,0.3)+normpdf((i-1)/Tmp,1,0.3) normpdf((i-1)/Tmp,0.33,0.3) normpdf((i-1)/Tmp,0.66,0.3)];
    color = color./max(color);
    if ~(get(handles.All,'Value')) && length(Spks) > eval(get(handles.MaxSp,'String')) 
    Spks = randsample(Spks,eval(get(handles.MaxSp,'String')));
    end
    plot(squeeze(Spikes{site}.Waveform(Spks,In,:))','color',color)
    hold on
    end
    hold off
    
else
if get(handles.CompPair,'Value')
if get(handles.All,'Value')
    plot(squeeze(Spikes{site}.Waveform(Spks,In,:))','Color',[0 1 0])
    hold on
    plot(squeeze(Spikes{site}.Waveform(Spks2,In,:))','Color',[1 0 1])
    hold off
else
    if length(Spks) > eval(get(handles.MaxSp,'String')) && length(Spks2) > eval(get(handles.MaxSp,'String'))
    plot(squeeze(Spikes{site}.Waveform(randsample(Spks,eval(get(handles.MaxSp,'String'))),In,:))','Color',[0 1 0])
    hold on
    plot(squeeze(Spikes{site}.Waveform(randsample(Spks2,eval(get(handles.MaxSp,'String'))),In,:))','Color',[1 0 1])
    hold off
    else
    plot(squeeze(Spikes{site}.Waveform(Spks,In,:))','Color',[0 1 0])
    hold on
    plot(squeeze(Spikes{site}.Waveform(Spks2,In,:))','Color',[1 0 1])
    hold off
    end
end
else
if get(handles.All,'Value')
    plot(squeeze(Spikes{site}.Waveform(Spks,In,:))')
else
    if length(Spks) > eval(get(handles.MaxSp,'String')) 
    plot(squeeze(Spikes{site}.Waveform(randsample(Spks,eval(get(handles.MaxSp,'String'))),In,:))')
    else
    plot(squeeze(Spikes{site}.Waveform(Spks,In,:))')
    end
end
end
end
PC = 0;
set(handles.pushbutton6,'Visible','off')
set(handles.points,'Visible','off')
Spks = find(Spikes{site}.Unit==unit);
PCtoD = zeros(length(Spks),1);

% hObject    handle to ChannelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChannelMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChannelMenu


% --- Executes during object creation, after setting all properties.
function ChannelMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChannelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global Spikes UnitType
In = get(handles.UnitMenu,'Value')-1;
Site = get(handles.SiteMenu,'Value')-1;
Spikes{Site}.Unit(Spikes{Site}.Unit==In) = 0;
for i = In+1:length(UnitType{Site})
    Spikes{Site}.Unit(Spikes{Site}.Unit==i)  = i-1;
end
UnitType{Site}(In:end-1) = UnitType{Site}(In+1:end);
UnitType{Site} = UnitType{Site}(1:end-1);
set(handles.UnitMenu,'Value',1);
set(handles.UnitChg,'Value',1);
set(handles.ChannelMenu,'Value',1);
set(handles.UType,'Value',1);
set(handles.UCompare,'Value',1);
UnitMenuStr{1} = 'unit';
for i = 1:length(UnitType{Site})
    UnitMenuStr{i+1} = int2str(i);
end
set(handles.UnitMenu,'String',UnitMenuStr);
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in UnitMenu.
function UnitMenu_Callback(hObject, eventdata, handles)
global UnitType Spikes PCs PC Mapping PCtoD
In = get(handles.UnitMenu,'Value')-1;
if In == 0
    return 
end
Site = get(handles.SiteMenu,'Value')-1;
set(handles.ChannelMenu,'Value',1);
set(handles.UCompare,'Value',1);
set(handles.UType,'Value',UnitType{Site}(In));

SpkPerSec = sum(Spikes{Site}.Unit == In)/(max(Spikes{Site}.TimeStamp(Spikes{Site}.Unit == In))/30000);
set(handles.SpkPerSec,'String',num2str(SpkPerSec));
UnitChgStr{1} = 'Unit to merge';
g = 1;
for i = 1:length(UnitType{Site})
    if i ~= In
    UnitChgStr{g+1} = int2str(i);
    g = g+1;
    end
end
set(handles.UnitChg,'Value',1);
set(handles.UnitChg,'String',UnitChgStr);
UnitChgStr{1} = 'Unit to Comp';
g = 1;
for i = 1:length(UnitType{Site})
    if i ~= In
    UnitChgStr{g+1} = int2str(i);
    g = g+1;
    end
end
set(handles.UCompare,'String',UnitChgStr);
Spks = find(Spikes{Site}.Unit == In);
PCtoD = zeros(length(Spks),1);

% hObject    handle to UnitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns UnitMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UnitMenu


% --- Executes during object creation, after setting all properties.
function UnitMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UnitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global UnitType Spikes PC PCtoD
unit = get(handles.UnitMenu,'Value')-1;
Site = get(handles.SiteMenu,'Value')-1;
newUnit = length(UnitType{Site})+1;
if ~PC
In = get(handles.ChannelMenu,'Value')-1;
if In == 0
    return 
end
Values = ginput(2);
Values(:,1) = round(Values(:,1));
Spks = find(Spikes{Site}.Unit== unit);
ConvertedSpks = 0;
for i =1:length(Spks)
    if not(isempty(find(Spikes{Site}.Waveform(Spks(i),In,min(Values(:,1)):max(Values(:,1))) > min(Values(:,2)) & Spikes{Site}.Waveform(Spks(i),In,min(Values(:,1)):max(Values(:,1))) < max(Values(:,2)))))
        Spikes{Site}.Unit(Spks(i)) = newUnit;
        ConvertedSpks = 1;
    end
end
else
    
Spks = find(Spikes{Site}.Unit== unit);
ConvertedSpks = 0;

if sum(PCtoD)>0
    Spikes{Site}.Unit(Spks(PCtoD)) = newUnit;
    ConvertedSpks = 1;
end

end
Units = length(unique(Spikes{Site}.Unit));

if length(UnitType{Site}) ~= Units
UnitType{Site} = [UnitType{Site};3];
end

UnitMenuStr{1} = 'unit';
for i = 1:length(UnitType{Site})
    UnitMenuStr{i+1} = int2str(i);
end

set(handles.UnitMenu,'String',UnitMenuStr);
set(handles.UnitMenu,'Value',1);
set(handles.UnitChg,'Value',1);
set(handles.ChannelMenu,'Value',1);
set(handles.UType,'Value',1);

        
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
global  Spikes UnitType

Site = get(handles.SiteMenu,'Value')-1;
In = get(handles.UnitMenu,'Value')-1;
if In < get(handles.UnitChg,'Value')-1
Spikes{Site}.Unit(find(Spikes{Site}.Unit==In)) = get(handles.UnitChg,'Value');
else
Spikes{Site}.Unit(find(Spikes{Site}.Unit==In)) = get(handles.UnitChg,'Value')-1;
end
for i = In+1:length(UnitType{Site})
   Spikes{Site}.Unit(find(Spikes{Site}.Unit==i)) = i-1;
end
UnitType{Site}(In:end-1) = UnitType{Site}(In+1:end);
UnitType{Site} = UnitType{Site}(1:end-1);
if In < get(handles.UnitChg,'Value')-1
UnitType{Site}(get(handles.UnitChg,'Value')) = 3;
else
UnitType{Site}(get(handles.UnitChg,'Value')-1) = 3;
end
set(handles.UnitMenu,'Value',1);
set(handles.ChannelMenu,'Value',1);
set(handles.UType,'Value',1);
set(handles.UnitChg,'Value',1);
UnitMenuStr{1} = 'unit';
for i = 1:length(UnitType{Site})
    UnitMenuStr{i+1} = int2str(i);
end
set(handles.UnitMenu,'String',UnitMenuStr);


% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in UnitChg.
function UnitChg_Callback(hObject, eventdata, handles)
% hObject    handle to UnitChg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns UnitChg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UnitChg


% --- Executes during object creation, after setting all properties.
function UnitChg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UnitChg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxSp_Callback(hObject, eventdata, handles)
% hObject    handle to MaxSp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxSp as text
%        str2double(get(hObject,'String')) returns contents of MaxSp as a double


% --- Executes during object creation, after setting all properties.
function MaxSp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxSp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in All.
function All_Callback(hObject, eventdata, handles)
% hObject    handle to All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of All


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
global UnitType Spikes PCs PC Mapping PCtoD PCVs
PCVs = []; 

Site = get(handles.SiteMenu,'Value') -1;
set(handles.ChannelMenu,'Value',1)
unit = get(handles.UnitMenu,'Value')-1;
if get(handles.allUnits,'Value')
    allSpks = find(Spikes{Site}.Unit>0);
else
if get(handles.CompPair,'Value')
    unit2 = get(handles.UCompare,'Value')-1;
    if unit2 >= unit
        unit2=unit2+1;
    end
    allSpks = find(Spikes{Site}.Unit==unit | Spikes{Site}.Unit==unit2);
else
PC = 1;
set(handles.pushbutton6,'Visible','on')
set(handles.points,'Visible','on')
allSpks = find(Spikes{Site}.Unit== unit);
end
end
UnitsPloted = unique(Spikes{Site}.Unit(allSpks));
chs = [get(handles.CH1,'Value') get(handles.CH2,'Value')];
PCVs = [];
for i = 1:length(UnitsPloted)
if UnitsPloted(i)>0
    
    if PCs(1)<4 || PCs(2)<4
    Waveforms = squeeze(Spikes{Site}.Waveform(allSpks,chs(1),:))';
    G{1}.PC = princomp(Waveforms.');
    Waveforms = squeeze(Spikes{Site}.Waveform(allSpks,chs(2),:))';
    G{2}.PC = princomp(Waveforms.');
    end
    PCVs = [];
    for j =1:2
    Spks = find(Spikes{Site}.Unit== UnitsPloted(i));
    Waveforms = squeeze(Spikes{Site}.Waveform(Spks,chs(j),:))';
    
    if PCs(j) < 4
    PCVs(j,:) = G{j}.PC(PCs(j),:)*Waveforms;
    end
    
    if PCs(j) == 4
    Waveforms = squeeze(Spikes{Site}.Waveform(Spks,chs(j),:))';
    for h = 1:length(Waveforms(1,:))
        Min = min(Waveforms(:,h));
        Minp = find(Waveforms(:,h) == Min);
        Max = max(Waveforms(1:Minp,h));
        PCVs(j,h) = Max;
    end
    end

    if PCs(j) == 5
    Waveforms = squeeze(Spikes{Site}.Waveform(Spks,chs(j),:))';
    for h = 1:length(Waveforms(1,:))
        Energy = sqrt(sum((diff(Waveforms(:,h))).^2));
        PCVs(j,h) = Energy;
    end
    end
    
    if PCs(j) == 6
    Waveforms = Spikes{Site}.Waveform(Spks,:,:);
    for h = 1:length(Waveforms(:,1,1))
    Wf = squeeze(Waveforms(h,:,:));
    Min = min(min(Wf));
    [Y MinP] = find(Wf == Min,1);
    Min = min(Wf(chs(j),MinP));
    PCVs(j,h) = abs(Min);
    end
    end
    if PCs(j) == 7
    Waveforms = squeeze(Spikes{Site}.Waveform(Spks,chs(j),:))';
    for h = 1:length(Waveforms(1,:))
        Min = min(Waveforms(:,h));
        MinP = find(Waveforms(:,h) == Min,1);
        Max = max(Waveforms(MinP:end,h));
        MaxP = find(Waveforms(:,h) == Max,1);
        PCVs(j,h) = (Max-Min)/abs(MaxP-MinP);
    end
    end
    
    if PCs(j) == 8
    Waveforms = squeeze(Spikes{Site}.Waveform(Spks,chs(j),:))';
    for h = 1:length(Waveforms(1,:))
    Min = min(Waveforms(:,h));
    MinP = find(Waveforms(:,h) == Min,1);
    PCVs(j,h) = MinP;
    end
    end
    if PCs(j) == 9
    Waveforms = squeeze(Spikes{Site}.Waveform(Spks,chs(j),:))';
    for h = 1:length(Waveforms(1,:))
        Min = min(Waveforms(:,h));
        Minp = find(Waveforms(:,h) == Min);
        Max = max(Waveforms(Minp:end,h));
        PCVs(j,h) = Max;
    end
    end
    
    end
    
    
Tmp = length(UnitsPloted)-1;
if Tmp == 0
   Tmp = 1;
end

color = [normpdf((i-1)/Tmp,0,0.3)+normpdf((i-1)/Tmp,1,0.3) normpdf((i-1)/Tmp,0.33,0.3) normpdf((i-1)/Tmp,0.66,0.3)];
color = color./max(color);
if not(get(handles.All,'Value')) && length(PCVs(1,:)) > eval(get(handles.MaxSp,'String'))
ToP = randsample([1:length(PCVs(1,:))],eval(get(handles.MaxSp,'String')) );
PCVsTOp = PCVs(:,ToP);
else
PCVsTOp = PCVs;
end
if i == 1
scatter(PCVsTOp(1,:),PCVsTOp(2,:),'CData',[0 0 0],'Marker','.','SizeData',eval(get(handles.DotSize,'String')))
hold on
else
scatter(PCVsTOp(1,:),PCVsTOp(2,:),'CData',color,'Marker','.','SizeData',eval(get(handles.DotSize,'String')))
hold on
end
end
end

hold off
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
global Analyzer  Spikes UnitType Events Mapping

% Save Sort File
save(['out\SpikesEphys\' get(handles.Animal,'string') '_' get(handles.Unit,'string') '-' get(handles.Exp,'string') ' ' 'Spikes'],'Spikes','UnitType','-append')

% Analyze responses

Params = Analyzer.L.param{1,1};
Reps = length(Analyzer.loops.conds{1,1}.repeats);
BReps =length(Analyzer.loops.conds{1,end}.repeats);
Conds = length(Analyzer.loops.conds);
if isequal(Analyzer.loops.conds{1,end}.symbol{1,1}, 'blank')
    Conds = length(Analyzer.loops.conds)-1;
    for r = 1:BReps
        Trial = Analyzer.loops.conds{1,end}.repeats{1,r}.trialno;
        TrialInfo(Trial,:) = [(zeros(length(Analyzer.loops.conds{1,end-1}.val(:,:)),1)-1).'  double(Events.Timestamp{1}(Trial,:))];
    end
end
for i = 1:Conds
    for r = 1:Reps
        Trial = Analyzer.loops.conds{1,i}.repeats{1,r}.trialno;
        TrialInfo(Trial,:) = [vertcat(Analyzer.loops.conds{1,i}.val{:,:}).' double(Events.Timestamp{1}(Trial,:))];
    end
    CondInfo(i,:)= vertcat(Analyzer.loops.conds{1,i}.val{:,:}).';
end
Parmass = Analyzer.L.param;
VarValDims = 1;
for i=1:length(Parmass(1,:))
    VarValDims = VarValDims*(length(eval(Parmass{1,i}{2})));
end

for site = 1:length(Spikes)
Units = unique(Spikes{site}.Unit);
Units = sort(Units);
if Units(1) == 0
    Units = Units(2:end);
end
Units = length(Units);
Data{site}.Spiking=cell(VarValDims, Units, Reps);
Data{site}.BSpiking=cell(Units, BReps);
Repetition = ones(VarValDims,1);
for T = 1:length(TrialInfo(:,1))
    if not(isequal(TrialInfo(T,1:length(Parmass(1,:))), (zeros(1,length(Parmass(1,:))) -1)))
    for Unit=1:Units
    VariableIndex = VarIndx(TrialInfo(T,1:length(Parmass(1,:))));
    Data{site}.Spiking{VariableIndex,Unit,Repetition(VariableIndex)} = double(Spikes{site}.TimeStamp(find(Spikes{site}.Unit == Unit & Spikes{site}.TimeStamp > TrialInfo(T,length(Parmass(1,:))+1) & Spikes{site}.TimeStamp < TrialInfo(T,length(Parmass(1,:))+4))))-TrialInfo(T,length(Parmass(1,:))+2);
    end
    Repetition(VariableIndex) = Repetition(VariableIndex)+1;
    end
end
Repetition = 1;
for T = 1:length(TrialInfo(:,1))
    if TrialInfo(T,1:length(Parmass(1,:))) == zeros(1,length(Parmass(1,:))) -1
    for Unit=1:Units
    Data{site}.BSpiking{Unit,Repetition} = double(Spikes{site}.TimeStamp(find(Spikes{site}.Unit == Unit & Spikes{site}.TimeStamp > TrialInfo(T,length(Parmass(1,:))+1) & Spikes{site}.TimeStamp < TrialInfo(T,length(Parmass(1,:))+4))))-TrialInfo(T,length(Parmass(1,:))+2);
    end
    Repetition = Repetition+1;
    end
end

Paramss = Analyzer.L.param;
FixIndx = zeros(length(Paramss(1,:)),1)-1;
for Unit = 1:length(UnitType{site})
%Calc Blank Responses
RepVal = [];
    for rep = 1:length(Data{site}.BSpiking(1,:))
        Spks = [Data{site}.BSpiking{Unit,rep}];
        RepVal(rep) =((sum(Spks>0 & Spks < 30000*Analyzer.P.param{1,3}{3})/(Analyzer.P.param{1,3}{3}))-((sum(Spks<0))/Analyzer.P.param{1,1}{3}));
    end
Data{site}.BRespMean(Unit) = mean(RepVal);
Data{site}.BRespVar(Unit) = std(RepVal)/sqrt(length(RepVal));
Data{site}.AllBResp(:,Unit) = RepVal;

%Calc responses for all conditions
for i = 1:length(Data{site}.Spiking(:,1,1))
RepVal = [];
    for rep = 1:length(Data{site}.Spiking(1,1,:))
        Spks = [Data{site}.Spiking{i,Unit,rep}];
        RepVal(rep) =((sum(Spks>0 & Spks < 30000*Analyzer.P.param{1,3}{3})/(Analyzer.P.param{1,3}{3}))-((sum(Spks<0))/Analyzer.P.param{1,1}{3}));
    end
Data{site}.RespVar(Unit,i) = std(RepVal)/sqrt(length(RepVal));
Data{site}.RespMean(Unit,i) = mean(RepVal);
Data{site}.AllResp(Unit,:,i) = RepVal;
end
end
end
RespFunc = [0 30000*Analyzer.P.param{1,3}{3} 0 -30000*Analyzer.P.param{1,1}{3}];
Values = cell(length(Paramss(1,:)),1);
for i = 1:length(Paramss(1,:))
    Params = Paramss{1,i};
    Variables{i} =Params{1,1};
    Values{i} = eval(Params{1,2});
end
save(['out/AnalyzedEphys/' get(handles.Animal,'string') '_' get(handles.Unit,'string') '-' get(handles.Exp,'string') ' ' 'Data'], 'Data', 'UnitType','Variables','Values','RespFunc','CondInfo')



% Write Table

if length(UnitType)>1
    warndlg('Table writting not programmed for more than one site')
    return
end
FileName = 'Ephys Data.xlsx';
Path = 'C:\Users\Nielsen Lab\Documents\';
[num,txt,raw] = xlsread([Path FileName],1);
Rows = [];
for i = 1:length(raw(:,1))
if isequal(raw(i,3),{[get(handles.Animal,'String') '_' get(handles.Unit,'String') '_' get(handles.Exp,'String')]})
Rows = [Rows i];
end
end
if isempty(Rows)
    errordlg('Experiment has not been added to table yet')
    return
end
if length(Rows) > 1
    warndlg('Table showing previous units, deleting old units')
    raw(Rows(2:end),:) = [];
    Rows = Rows(1);
end
raw = [raw(1:Rows(1),:);cell(length(UnitType{1})-1,length(raw(1,:)));raw(Rows(1)+1:end,:)];
raw(Rows(1):Rows(1)+length(UnitType{1})-1,1) = {get(handles.Animal,'String')};
raw(Rows(1):Rows(1)+length(UnitType{1})-1,3) = {[get(handles.Animal,'String') '_' get(handles.Unit,'String') '_' get(handles.Exp,'String')]};
raw(Rows(1):Rows(1)+length(UnitType{1})-1,10) = cell(length(UnitType{1}),1);
raw(Rows(1):Rows(1)+length(UnitType{1})-1,18) = raw(Rows(1),18);
raw(Rows(1):Rows(1)+length(UnitType{1})-1,2) = raw(Rows(1),2);
raw(Rows(1):Rows(1)+length(UnitType{1})-1,4) = raw(Rows(1),4);
raw(Rows(1):Rows(1)+length(UnitType{1})-1,5) = raw(Rows(1),5);
raw(Rows(1):Rows(1)+length(UnitType{1})-1,6) = raw(Rows(1),6);
raw(Rows(1):Rows(1)+length(UnitType{1})-1,11:14) = repmat([{3} {date} {3} {date}],length(UnitType{1}),1);
for i = 1:length(UnitType{1})
    raw(Rows(1)+i-1,7) = {i};
    raw(Rows(1)+i-1,8) = {UnitType{1}(i)};
    Spks = Spikes{1}.Unit == i;
    TStmps = double(Spikes{1}.TimeStamp(Spks));
    ISI = TStmps(2:end) - TStmps(1:end-1);
    Cont = sum(ISI<45)/length(ISI);
    raw(Rows(1)+i-1,9) = {Cont};
end
xlswrite([Path FileName],raw,1)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
global PCVs PCtoD
if get(handles.points,'value') == 1
PCLs = ginput(3);
else
PCLs = ginput(get(handles.points,'value')+1);
end
if sum(PCtoD) == 0
PCtoD = inpolygon(PCVs(1,:).',PCVs(2,:).',PCLs(:,1),PCLs(:,2));
else
PCtoD = PCtoD.*(inpolygon(PCVs(1,:).',PCVs(2,:).',PCLs(:,1),PCLs(:,2)));
end
if sum(PCtoD)>0
axes(handles.axes1)
scatter(PCVs(1,find(PCtoD == 0)),PCVs(2,find(PCtoD == 0)),'CData',[0 0 1],'Marker','.','SizeData',eval(get(handles.DotSize,'String')))
hold on
scatter(PCVs(1,find(PCtoD == 1)),PCVs(2,find(PCtoD == 1)),'CData',[1 0 0],'Marker','.','SizeData',eval(get(handles.DotSize,'String')))
hold off
end
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in points.
function points_Callback(hObject, eventdata, handles)
% hObject    handle to points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns points contents as cell array
%        contents{get(hObject,'Value')} returns selected item from points


% --- Executes during object creation, after setting all properties.
function points_CreateFcn(hObject, eventdata, handles)
% hObject    handle to points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CP12.
function CP12_Callback(hObject, eventdata, handles)
global PCs
if ~(get(handles.CP12,'Value'));
set(handles.CP12,'Value',1);
end
set(handles.CP13,'Value',0);
set(handles.CP23,'Value',0);
PCs = [1 2];

% hObject    handle to CP12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CP12


% --- Executes on button press in CP23.
function CP23_Callback(hObject, eventdata, handles)
global PCs
if ~(get(handles.CP23,'Value'));
set(handles.CP23,'Value',1);
end
set(handles.CP13,'Value',0);
set(handles.CP12,'Value',0);
PCs = [1 3];
% hObject    handle to CP23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CP23


% --- Executes on button press in CP13.
function CP13_Callback(hObject, eventdata, handles)
global PCs
if ~(get(handles.CP13,'Value'));
set(handles.CP13,'Value',1);
end
set(handles.CP12,'Value',0);
set(handles.CP23,'Value',0);
PCs = [2 3];
% hObject    handle to CP13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CP13


% --- Executes on button press in allUnits.
function allUnits_Callback(hObject, eventdata, handles)
global TetrodeBool
if get(handles.allUnits,'Value')
    if get(handles.CompPair,'Value')
        set(handles.CompPair,'Value',0)
        set(handles.UCompare,'Visible','off')
    end
    set(handles.UnitMenu,'Value',1)
    set(handles.UnitMenu,'Visible','off')
    ChannelMenuStr{1} = 'Channel';
    if TetrodeBool
    Chs = [1:4];
    else
    Chs = [1:16];
    end
    ChannelMenuStr{2} = Chs(1);
    for i = 2:length(Chs)
    ChannelMenuStr{i+1} = int2str(Chs(i));
    end
    set(handles.ChannelMenu,'Value',1);
    set(handles.UType,'Value',1);
    set(handles.ChannelMenu,'String',ChannelMenuStr);
else
    set(handles.UnitMenu,'Visible','on')
    set(handles.ChannelMenu,'Value',1);
    set(handles.UType,'Value',1);
    set(handles.ChannelMenu,'String','Channel');
end
    
% hObject    handle to allUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allUnits


% --- Executes on selection change in UType.
function UType_Callback(hObject, eventdata, handles)
global UnitType
In = get(handles.UnitMenu,'Value')-1;
Site = get(handles.SiteMenu,'Value')-1;
UnitType{Site}(In)=get(handles.UType,'Value');

% hObject    handle to UType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns UType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UType


% --- Executes during object creation, after setting all properties.
function UType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CompPair.
function CompPair_Callback(hObject, eventdata, handles)
if get(handles.CompPair,'Value')
    set(handles.UCompare,'Visible','On')
    if get(handles.allUnits,'Value')
       set(handles.allUnits,'Value',0)
       set(handles.UnitMenu,'Visible','on')
    end
else
    set(handles.UCompare,'Visible','Off')
end
% hObject    handle to CompPair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CompPair


% --- Executes on selection change in UCompare.
function UCompare_Callback(hObject, eventdata, handles)
% hObject    handle to UCompare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns UCompare contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UCompare


% --- Executes during object creation, after setting all properties.
function UCompare_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UCompare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Var2.
function Var2_Callback(hObject, eventdata, handles)
global PCs
PCs(2) = get(handles.Var2,'Value');
% hObject    handle to Var2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Var2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Var2


% --- Executes during object creation, after setting all properties.
function Var2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Var2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Var1.
function Var1_Callback(hObject, eventdata, handles)
global PCs
PCs(1) = get(handles.Var1,'Value');
% hObject    handle to Var1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Var1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Var1


% --- Executes during object creation, after setting all properties.
function Var1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Var1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CH1.
function CH1_Callback(hObject, eventdata, handles)
% hObject    handle to CH1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CH1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CH1


% --- Executes during object creation, after setting all properties.
function CH1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CH1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CH2.
function CH2_Callback(hObject, eventdata, handles)
% hObject    handle to CH2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CH2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CH2


% --- Executes during object creation, after setting all properties.
function CH2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CH2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Animal_Callback(hObject, eventdata, handles)
% hObject    handle to Animal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Animal as text
%        str2double(get(hObject,'String')) returns contents of Animal as a double


% --- Executes during object creation, after setting all properties.
function Animal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Animal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Unit_Callback(hObject, eventdata, handles)
% hObject    handle to Unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Unit as text
%        str2double(get(hObject,'String')) returns contents of Unit as a double


% --- Executes during object creation, after setting all properties.
function Unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Exp_Callback(hObject, eventdata, handles)
% hObject    handle to Exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Exp as text
%        str2double(get(hObject,'String')) returns contents of Exp as a double


% --- Executes during object creation, after setting all properties.
function Exp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadSpks.
function LoadSpks_Callback(hObject, eventdata, handles)
clear UnitType Spikes PCs PC Mapping
global UnitType Spikes PCs PC Mapping Events Analyzer
eval(strcat('load(''','out\SpikesEphys\',get(handles.Animal,'string'),'_',get(handles.Unit,'string'),'-',get(handles.Exp,'string'),' Spikes.mat'',''-mat'')'))
eval(strcat('load(''Z:\Ephys\AnalyzerFiles\',get(handles.Animal,'string'), '\', get(handles.Animal,'string'),'_',get(handles.Unit,'string'),'_',get(handles.Exp,'string'),'.analyzer'',''-mat'')'))

SiteMenuStr{1} = 'site';
for i = 1:length(Spikes)
    SiteMenuStr{i+1} = int2str(i);
end
set(handles.SiteMenu,'String',SiteMenuStr);
set(handles.SiteMenu,'Value',1);
set(handles.UnitMenu,'Value',1);

PCs = [1 1];
PC = 0;
% hObject    handle to LoadSpks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in SiteMenu.
function SiteMenu_Callback(hObject, eventdata, handles)
global UnitType Spikes PCs PC Mapping
input = get(handles.SiteMenu,'Value')-1;
if input<1
    return
end
UnitMenuStr{1} = 'unit';
for i = 1:length(UnitType{input})
    UnitMenuStr{i+1} = int2str(i);
end
Chs = Mapping(1,(find(Mapping(2,:) == input)));
Channels =[];
Str = [];
for i =1:length(Chs)
Channels = strcat(Channels, '[',int2str(Chs(i)),'] ');
set(handles.Chs,'String',Channels);
Str{i} = int2str(Chs(i));
end
set(handles.CH1,'String',Str);
set(handles.CH2,'String',Str);
Str = [];
Str{1} = 'Channel';
for i =1:length(Chs)
Str{i+1} = int2str(Chs(i));
end
set(handles.ChannelMenu,'String',Str);
set(handles.CH1,'Visible','on');
set(handles.CH2,'Visible','on');
set(handles.UnitMenu,'String',UnitMenuStr);
set(handles.UnitMenu,'Value',1);
% hObject    handle to SiteMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SiteMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SiteMenu


% --- Executes during object creation, after setting all properties.
function SiteMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SiteMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DotSize_Callback(hObject, eventdata, handles)
% hObject    handle to DotSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DotSize as text
%        str2double(get(hObject,'String')) returns contents of DotSize as a double


% --- Executes during object creation, after setting all properties.
function DotSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DotSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
