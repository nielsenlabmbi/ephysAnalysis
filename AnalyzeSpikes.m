function varargout = AnalyzeSpikes(varargin)
% ANALYZESPIKES MATLAB code for AnalyzeSpikes.fig
%      ANALYZESPIKES, by itself, creates a new ANALYZESPIKES or raises the existing
%      singleton*.
%
%      H = ANALYZESPIKES returns the handle to a new ANALYZESPIKES or the handle to
%      the existing singleton*.
%
%      ANALYZESPIKES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZESPIKES.M with the given input arguments.
%
%      ANALYZESPIKES('Property','Value',...) creates a new ANALYZESPIKES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AnalyzeSpikes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AnalyzeSpikes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AnalyzeSpikes

% Last Modified by GUIDE v2.5 17-Sep-2015 16:55:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AnalyzeSpikes_OpeningFcn, ...
                   'gui_OutputFcn',  @AnalyzeSpikes_OutputFcn, ...
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


% --- Executes just before AnalyzeSpikes is made visible.
function AnalyzeSpikes_OpeningFcn(hObject, eventdata, handles, varargin)
global UnitIndx Analyzer FixIndx BSpiking BRespMean BRespVar RespVar RespMean Spiking AllResp AllBResp
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AnalyzeSpikes (see VARARGIN)

% Choose default command line output for AnalyzeSpikes
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
axes(handles.axes3)
axis off

% UIWAIT makes AnalyzeSpikes wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AnalyzeSpikes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in unit.
function unit_Callback(hObject, eventdata, handles)
global Data UnitType CondInfo Variables Values ValsMax tmpvalsvar
set(handles.popupmenu5,'Value',1)
set(handles.FixedVars,'String','')
site = get(handles.SiteMenu,'Value')-1;
if get(handles.unit,'Value')==1
    return
end
set(handles.Channel,'Value',1);

In = get(handles.unit,'Value')-1;
tmp = get(handles.uitable2,'Data');

set(handles.uitable2,'Data',repmat({''},length(Variables),4));

for i = 1:length(Variables)
    tmp = get(handles.uitable2,'Data');
    tmp(i,1) =Variables(i);
    set(handles.uitable2,'Data',tmp)
end
set(handles.Var1,'String',tmp(:,1))
set(handles.Var2,'String',tmp(:,1))
set(handles.popupmenu5,'String',tmp(:,1))
tmpvalsvar = Values{1};
tmp = cell(length(tmpvalsvar),1);
for i = 1:length(tmpvalsvar)
    tmp{i} = num2str(tmpvalsvar(i));
end
set(handles.popupmenu6,'String',tmp)
Resp = Data{site}.RespMean(In,:);
RespVr = Data{site}.RespVar(In,:);
ValsMax =ReverseVarIndx(find(Resp == max(Resp),1));
Vals = repmat({''},length(ValsMax),1);
for i = 1:length(ValsMax)
    Vals{i} = mat2str(ValsMax(i));
end
tmp = get(handles.uitable2,'Data');
tmp(:,2) =Vals;
set(handles.uitable2,'Data',tmp);
TunCurs = cell(length(ValsMax),1);
TunCursVar = cell(length(ValsMax),1);
for i = 1:length(ValsMax)
    tmpvalsvar = Values{i};
        tmpTunCur = [];
        tmpTunCurVar = [];
    for Val = 1:length(tmpvalsvar)
        tmpVals = ValsMax;
        tmpVals(i) = tmpvalsvar(Val);
        tmpTunCur(Val) = Resp(VarIndx(tmpVals));
        tmpTunCurVar(Val) = RespVr(VarIndx(tmpVals));
    end
    MxMnMx{i} = (max(tmpTunCur)-min(tmpTunCur))/max(tmpTunCur);
    HalfPoint = abs(tmpTunCur -(max(tmpTunCur)-min(tmpTunCur))/2);
    Width{i} = abs(find(tmpTunCur == max(tmpTunCur)) - find(HalfPoint == min(HalfPoint)));
end
tmp = get(handles.uitable2,'Data');
tmp(:,4) =MxMnMx;
tmp(:,3) =Width;
set(handles.uitable2,'Data',tmp);
% hObject    handle to unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns unit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from unit


% --- Executes during object creation, after setting all properties.
function unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Var1.
function Var1_Callback(hObject, eventdata, handles)
global FixIndx Analyzer
Paramss = Analyzer.L.param;
FixIndx = zeros(length(Paramss(1,:)),1)-1;
set(handles.FixedVars,'String','')
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


% --- Executes on selection change in Var2.
function Var2_Callback(hObject, eventdata, handles)
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


% --- Executes on button press in TD.
function TD_Callback(hObject, eventdata, handles)
if get(handles.TD,'value')
    set(handles.Var2,'Visible','On')
else
    set(handles.Var2,'Visible','Off')
end
% hObject    handle to TD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TD


% --- Executes on selection change in Channel.
function Channel_Callback(hObject, eventdata, handles)
global   Spikes UnitType
unit = get(handles.unit,'Value')-1;
In = get(handles.Channel,'Value')-1;
site = get(handles.SiteMenu,'Value')-1;
if In == 0 || unit == 0 || site == 0
    return 
end
axes(handles.axes4)
Spks = find(Spikes{site}.Unit==unit);
Mean = mean(double(squeeze(Spikes{site}.Waveform(Spks,In,:))),1);
Upper = Mean + std(double(squeeze(Spikes{site}.Waveform(Spks,In,:))))/2;
Lower = Mean - std(double(squeeze(Spikes{site}.Waveform(Spks,In,:))))/2;
plot(Mean);
hold on
plot(Upper,'color',[1 0 0])
plot(Lower,'color',[1 0 0])
hold off
% hObject    handle to Channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Channel


% --- Executes during object creation, after setting all properties.
function Channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Rasterplot.
function Rasterplot_Callback(hObject, eventdata, handles)
% hObject    handle to Rasterplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Rasterplot


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in range.
function range_Callback(hObject, eventdata, handles)
if get(handles.range,'Value')
    set(handles.from,'Visible','on')
    set(handles.to,'Visible','on')
else
    set(handles.from,'Visible','off')
    set(handles.to,'Visible','off')
    set(handles.autocorr,'Value',0)
    set(handles.ISI,'Value',1)
end
% hObject    handle to range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of range



function from_Callback(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of from as text
%        str2double(get(hObject,'String')) returns contents of from as a double


% --- Executes during object creation, after setting all properties.
function from_CreateFcn(hObject, eventdata, handles)
% hObject    handle to from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function to_Callback(hObject, eventdata, handles)
% hObject    handle to to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of to as text
%        str2double(get(hObject,'String')) returns contents of to as a double


% --- Executes during object creation, after setting all properties.
function to_CreateFcn(hObject, eventdata, handles)
% hObject    handle to to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autocorr.
function autocorr_Callback(hObject, eventdata, handles)
if get(handles.autocorr,'Value')
    set(handles.ISI,'Value',0)
    set(handles.range,'Value',1)
    set(handles.from,'Visible','on')
    set(handles.to,'Visible','on')
else
    set(handles.ISI,'Value',1)
end
% hObject    handle to autocorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autocorr


% --- Executes on button press in ISI.
function ISI_Callback(hObject, eventdata, handles)
if get(handles.ISI,'Value')
    set(handles.autocorr,'Value',0)
else
    set(handles.autocorr,'Value',1)
end
% hObject    handle to ISI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ISI


% --- Executes on button press in Log.
function Log_Callback(hObject, eventdata, handles)
if get(handles.Log,'Value')
    set(handles.axes1,'XScale','log');
else
    set(handles.axes1,'XScale','linear');
end
% hObject    handle to Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Log


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global Data UnitType Spikes
In = get(handles.unit,'Value')-1;
site = get(handles.SiteMenu,'Value')-1;
if In == 0 || site == 0
    return
end
axes(handles.axes1)
Spks = find(Spikes{site}.Unit == In);
TStmps = double(Spikes{site}.TimeStamp(Spks));
ISI = TStmps(2:end) - TStmps(1:end-1);
if get(handles.autocorr,'Value')
binsize = floor((str2double(get(handles.to,'String'))-str2double(get(handles.from,'String')))/15);
if binsize == 0
    binsize = 1;
end
ACF = zeros(length(str2double(get(handles.from,'String')):binsize:str2double(get(handles.to,'String'))),1);
indx = 1;
for i = [str2double(get(handles.from,'String')):binsize:str2double(get(handles.to,'String'))]
    for L = 1:binsize*30
    ACF(indx) = ACF(indx)+length(intersect(TStmps,(TStmps-i*30-L)));
    end
    indx = indx + 1;
end
plot(1000./(str2double(get(handles.from,'String')):binsize:str2double(get(handles.to,'String'))),((double(ACF)*double(Spikes{site}.TimeStamp(end)))/length(Spks)).')
else
if get(handles.range,'Value')
[V C] = hist(ISI(find(ISI>str2double(get(handles.from,'String'))*30 & ISI<str2double(get(handles.to,'String'))*30)),50);
else
[V C] = (hist(ISI,50));
end
bar(C./30,V);
end
if get(handles.Log,'Value')
    set(handles.axes1,'XScale','log');
end
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global Analyzer TunCur TunCurVar Data FixIndx ValsMax tmpvalsvar
TunCur = [];
TunCurVar = [];
Parmass = Analyzer.L.param;
tmpvalsvar = eval(Parmass{1,get(handles.Var1,'Value')}{2});
tmpvalsvar = sort(tmpvalsvar);
Vals = FixIndx;
Vals(FixIndx<0) = ValsMax(FixIndx<0);
In = get(handles.unit,'Value')-1;
site = get(handles.SiteMenu,'Value')-1;
Resp = Data{site}.RespMean(In,:);
RespVr = Data{site}.RespVar(In,:);
if get(handles.TD,'Value')
   set(handles.axes3,'Visible','on')
   set(handles.ymin,'Visible','on')
   set(handles.ymax,'Visible','on')
   set(handles.xmin,'Visible','on')
   set(handles.xmax,'Visible','on')
   set(handles.scalemin,'Visible','on')
   set(handles.scalemax,'Visible','on')
   tmpvalsvar2 = eval(Parmass{1,get(handles.Var2,'Value')}{2});
   tmpvalsvar2 = sort(tmpvalsvar2);
   TunCur = zeros(length(tmpvalsvar),length(tmpvalsvar2));
   for Val2 = 1:length(tmpvalsvar2)
   tmpVals = Vals;
   tmpVals(get(handles.Var2,'Value')) = tmpvalsvar2(Val2);
   for Val = 1:length(tmpvalsvar)
   tmpVals(get(handles.Var1,'Value')) = tmpvalsvar(Val);
   TunCur(Val,Val2) = Resp(VarIndx(tmpVals));
   end
   end
   Mn = min(min(TunCur));
   Mx = max(max(TunCur));
   TunCur(TunCur<0) = TunCur(TunCur<0)*(5./abs(Mn));
   TunCur(TunCur>0) = TunCur(TunCur>0)*(5./abs(Mx));
   TunCur = round(TunCur);
   TunCur = TunCur+6;
   ColorCode = redbluecmap;
   for i = 1:length(TunCur(1,:))
       for h = 1:length(TunCur(:,1))
           HtMap(h,i,:) = ColorCode(TunCur(end-h+1,i),:);
       end
   end
   axes(handles.axes2)
   hold off
   image(HtMap);
   axis off
   axes(handles.axes3)
   ColorCode = imrotate(reshape(ColorCode,11,1,3),180);
   image(ColorCode);
   axis off
   set(handles.ymin,'String',num2str(min(tmpvalsvar)))
   set(handles.ymax,'String',num2str(max(tmpvalsvar)))
   set(handles.xmin,'String',num2str(min(tmpvalsvar2)))
   set(handles.xmax,'String',num2str(max(tmpvalsvar2)))
   set(handles.scalemin,'String',Mn)
   set(handles.scalemax,'String',Mx)
else
   set(handles.ymin,'Visible','off')
   set(handles.ymax,'Visible','off')
   set(handles.xmin,'Visible','off')
   set(handles.xmax,'Visible','off')
   set(handles.scalemin,'Visible','off')
   set(handles.scalemax,'Visible','off')
   axes(handles.axes3)
   cla reset
   set(handles.axes3,'Visible','off')
    for Val = 1:length(tmpvalsvar)
    tmpVals = Vals;
    tmpVals(get(handles.Var1,'Value')) = tmpvalsvar(Val);
    TunCur(Val) = Resp(VarIndx(tmpVals));
    TunCurVar(Val) = RespVr(VarIndx(tmpVals));
    end
    axes(handles.axes2)
    hold
    hold off
    plot(tmpvalsvar,TunCur)
    hold all
    set(handles.axes2,'ButtonDownFcn',@(hObject, eventdata)click(hObject, eventdata, guidata(hObject)))
    plot(tmpvalsvar,TunCur + TunCurVar,'color',[1 0 0])
    plot(tmpvalsvar,TunCur - TunCurVar,'color',[1 0 0])
end
    
    
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
global FixIndx S
tmp = get(handles.popupmenu6,'String');
tmp2 =  get(handles.popupmenu5,'String');
FixIndx(get(handles.popupmenu5,'Value')) = str2double(tmp{get(handles.popupmenu6,'Value')});
fixvalues = find(FixIndx >=0);
S={};
for i =1:length(fixvalues)
String = strcat(tmp2{fixvalues(i)},'=',num2str(FixIndx(fixvalues(i))));
S{end+1} = String;
end
set(handles.FixedVars,'String',S);
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
global Analyzer
Parmass = Analyzer.L.param;
tmpvalsvar = eval(Parmass{1,get(handles.popupmenu5,'Value')}{2});
tmp = cell(length(tmpvalsvar),1);
for i = 1:length(tmpvalsvar)
    tmp{i} = num2str(tmpvalsvar(i));
end
set(handles.popupmenu6,'Value',1)
set(handles.popupmenu6,'String',tmp)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on mouse press over axes background.
function click(hObject, eventdata, handles)
global tmpvalsvar Analyzer FixIndx ValsMax Data condition
site = get(handles.SiteMenu,'Value')-1; 
reps = length(Data{site}.Spiking(1,1,:));
point = get(handles.axes2,'CurrentPoint');
tmp = abs(tmpvalsvar-point(1,1));
point = (find(tmp==min(tmp)));
figure(1)
subplot(2,1,1)
hold off
plot([0 0], [0 2+reps], 'color','red')
hold on
StimDur = Analyzer.P.param{1,3}{3};
plot([StimDur StimDur], [0 2+reps], 'color','green')
Vals = FixIndx;
Vals(FixIndx<0) = ValsMax(FixIndx<0);
Vals(get(handles.Var1,'Value')) = tmpvalsvar(point);
condition = VarIndx(Vals); 
Unit = get(handles.unit,'Value')-1;
allSpks = [];
Spikes = [];
for i = 1:reps
    Spikes = Data{site}.Spiking{condition,Unit,i}./30000;
    allSpks = [allSpks;Spikes(:)];
    for S = 1:length(Spikes)
        plot([Spikes(S) Spikes(S)],[i i+1])
    end
end
subplot(2,1,2)
Value = [];
Bins = [min(allSpks):0.01:max(allSpks)];
for i = 1:length(Bins)
    Value(i) = sum(normpdf(allSpks,Bins(i),0.01));
end
hold off
plot(Bins,Value)
hold on
plot([0 0], [0 max(Value)], 'color','red')
plot([StimDur StimDur], [0 max(Value)], 'color','green')
% hObject   
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    



function animal_Callback(hObject, eventdata, handles)
% hObject    handle to animal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of animal as text
%        str2double(get(hObject,'String')) returns contents of animal as a double


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



function funit_Callback(hObject, eventdata, handles)
% hObject    handle to funit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of funit as text
%        str2double(get(hObject,'String')) returns contents of funit as a double


% --- Executes during object creation, after setting all properties.
function funit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to funit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function exp_Callback(hObject, eventdata, handles)
% hObject    handle to exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exp as text
%        str2double(get(hObject,'String')) returns contents of exp as a double


% --- Executes during object creation, after setting all properties.
function exp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadSpks.
function LoadSpks_Callback(hObject, eventdata, handles)
clear UnitType Data CondInfo Variables Values Analyzer Mapping
global UnitType Data CondInfo Variables Values Analyzer Mapping Spikes
eval(strcat('load(''','out\AnalyzedEphys\',get(handles.animal,'string'),'_',get(handles.funit,'string'),'-',get(handles.exp,'string'),' Data.mat'',''-mat'')'))
eval(strcat('load(''','out\SpikesEphys\',get(handles.animal,'string'),'_',get(handles.funit,'string'),'-',get(handles.exp,'string'),' Spikes.mat'',''-mat'')'))
eval(strcat('load(''Z:\Ephys\AnalyzerFiles\',get(handles.animal,'string'), '\', get(handles.animal,'string'),'_',get(handles.funit,'string'),'_',get(handles.exp,'string'),'.analyzer'',''-mat'')'))
SiteMenuStr{1} = 'site';
for i = 1:length(Data)
    SiteMenuStr{i+1} = int2str(i);
end
set(handles.SiteMenu,'String',SiteMenuStr);
set(handles.SiteMenu,'Value',1);
set(handles.unit,'Value',1);

% hObject    handle to LoadSpks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in siteMenu.
function SiteMenu_Callback(hObject, eventdata, handles)
global UnitType Spikes Mapping
input = get(handles.SiteMenu,'Value')-1;
if input<1
    return
end
UnitMenuStr{1} = 'unit';
for i = 1:length(UnitType{input})
    UnitMenuStr{i+1} = int2str(i);
end
Chs = Mapping(1,(find(Mapping(2,:) == input)));
Str = [];
Str{1} = 'Channel';
for i =1:length(Chs)
Str{i+1} = int2str(Chs(i));
end
set(handles.Channel,'String',Str);
set(handles.unit,'String',UnitMenuStr);
set(handles.unit,'Value',1);
% hObject    handle to siteMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns siteMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from siteMenu


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
