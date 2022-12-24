function varargout = PIDSimulation(varargin)
% PIDSIMULATION MATLAB code for PIDSimulation.fig
%      PIDSIMULATION, by itself, creates a new PIDSIMULATION or raises the existing
%      singleton*.
%
%      H = PIDSIMULATION returns the handle to a new PIDSIMULATION or the handle to
%      the existing singleton*.
%
%      PIDSIMULATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PIDSIMULATION.M with the given input arguments.
%
%      PIDSIMULATION('Property','Value',...) creates a new PIDSIMULATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PIDSimulation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PIDSimulation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PIDSimulation

% Last Modified by GUIDE v2.5 03-Apr-2018 19:09:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PIDSimulation_OpeningFcn, ...
                   'gui_OutputFcn',  @PIDSimulation_OutputFcn, ...
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


% --- Executes just before PIDSimulation is made visible.
function PIDSimulation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PIDSimulation (see VARARGIN)

% Choose default command line output for PIDSimulation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PIDSimulation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PIDSimulation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;


% --- Executes on button press in simplepid.
function simplepid_Callback(hObject, eventdata, handles)
% hObject    handle to simplepid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in fuzzypid.
function fuzzypid_Callback(hObject, eventdata, handles)
% hObject    handle to fuzzypid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in simplepidbtn.
function simplepidbtn_Callback(hObject, eventdata, handles)
% hObject    handle to simplepidbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[E,EE,iii,A,B,C]= simplepid();
axes(handles.auv)
plot3(EE(iii,1),EE(iii,2),EE(iii,3),'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Simple PID AUV Track');
plot3(A,B,C,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
A1=EE(:,1);

A2=A;
B1=EE(:,2);
B2=B;
C1=EE(:,3);
C2=C;
axes(handles.axes4)
plot(A1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Direct Graph');
plot(A2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
axes(handles.axes5)
plot(B1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Trim Graph');
plot(B2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
axes(handles.axes6)
plot(C1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Depth Graph');
plot(C2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;



% --- Executes on button press in fuzzypidbtn.
function fuzzypidbtn_Callback(hObject, eventdata, handles)
% hObject    handle to fuzzypidbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[E,EE,iii,A,B,C]= fuzzypid();
plot3(EE(iii,1),EE(iii,2),EE(iii,3),'-');
axes(handles.auv)
plot3(EE(iii,1),EE(iii,2),EE(iii,3),'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('Fuzzy logics PID AUV Track');
plot3(A,B,C,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
A1=EE(:,1);

A2=A;
B1=EE(:,2);
B2=B;
C1=EE(:,3);
C2=C;
axes(handles.axes4)
plot(A1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Direct Graph');
plot(A2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
axes(handles.axes5)
plot(B1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Trim Graph');
plot(B2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
axes(handles.axes6)
plot(C1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Depth Graph');
plot(C2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;




% --- Executes when figure1 is resized.


% --- Executes on button press in fuzzynnpidbtn.
function fuzzynnpidbtn_Callback(hObject, eventdata, handles)
% hObject    handle to fuzzynnpidbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[E,EE,iii,A,B,C]= nnfuzzypid();
axes(handles.auv)
plot3(EE(iii,1),EE(iii,2),EE(iii,3),'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Fuzzy logics & Neural Network PID AUV Track');
plot3(A,B,C,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
A1=EE(:,1);

A2=A;
B1=EE(:,2);
B2=B;
C1=EE(:,3);
C2=C;
axes(handles.axes4)
plot(A1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Direct Graph');
plot(A2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
axes(handles.axes5)
plot(B1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Trim Graph');
plot(B2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
axes(handles.axes6)
plot(C1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Depth Graph');
plot(C2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;




% --- Executes on button press in nnpidbtn.
function nnpidbtn_Callback(hObject, eventdata, handles)
% hObject    handle to nnpidbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[E,EE,iii,A,B,C]= nnpid();
axes(handles.auv)
plot3(EE(iii,1),EE(iii,2),EE(iii,3),'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Neural Network PID AUV Track');
plot3(A,B,C,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
A1=EE(:,1);

A2=A;
B1=EE(:,2);
B2=B;
C1=EE(:,3);
C2=C;
axes(handles.axes4)
plot(A1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Direct Graph');
plot(A2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
axes(handles.axes5)
plot(B1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Trim Graph');
plot(B2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;
axes(handles.axes6)
plot(C1,'-');
grid on;
set(gca,'zdir','reverse');
hold on;
title('	Depth Graph');
plot(C2,'r');
legend('actual','expected');%Actual ',' Expected
hold off;


% --- Executes on button press in Directfisbtn.
function Directfisbtn_Callback(hObject, eventdata, handles)
% hObject    handle to Directfisbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fuzzyLogicDesigner('fuzzy_piddirect')




% --- Executes on button press in nnbtn.
function nnbtn_Callback(hObject, eventdata, handles)
% hObject    handle to nnbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Status = nn_training();


% --- Executes on button press in trimfisbtn.
function trimfisbtn_Callback(hObject, eventdata, handles)
% hObject    handle to trimfisbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fuzzyLogicDesigner('fuzzy_pid3')


% --- Executes on button press in deepfisbtn.
function deepfisbtn_Callback(hObject, eventdata, handles)
% hObject    handle to deepfisbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fuzzyLogicDesigner('fuzzy_piddepth')


% --- Executes on button press in nnkpbtn.
function nnkpbtn_Callback(hObject, eventdata, handles)
% hObject    handle to nnkpbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nnkp_training();



% --- Executes on button press in nnkibtn.
function nnkibtn_Callback(hObject, eventdata, handles)
% hObject    handle to nnkibtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nnki_training();



% --- Executes on button press in nnkdbtn.
function nnkdbtn_Callback(hObject, eventdata, handles)
% hObject    handle to nnkdbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nnkd_training();
