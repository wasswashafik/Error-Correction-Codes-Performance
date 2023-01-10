function varargout = projectgui(varargin)
% PROJECTGUI MATLAB code for projectgui.fig
%      PROJECTGUI, by itself, creates a new PROJECTGUI or raises the existing
%      singleton*.
%
%      H = PROJECTGUI returns the handle to a new PROJECTGUI or the handle to
%      the existing singleton*.
%
%      PROJECTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECTGUI.M with the given input arguments.
%
%      PROJECTGUI('Property','Value',...) creates a new PROJECTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before projectgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to projectgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help projectgui

% Last Modified by GUIDE v2.5 30-Apr-2015 20:02:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @projectgui_OpeningFcn, ...
                   'gui_OutputFcn',  @projectgui_OutputFcn, ...
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


% --- Executes just before projectgui is made visible.
function projectgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to projectgui (see VARARGIN)

% Choose default command line output for projectgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projectgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = projectgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function nsym_Callback(hObject, eventdata, handles)
% hObject    handle to nsym (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nsym as text
%        str2double(get(hObject,'String')) returns contents of nsym as a double


% --- Executes during object creation, after setting all properties.
function nsym_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nsym (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in modulation.
function modulation_Callback(hObject, eventdata, handles)
% hObject    handle to modulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modulation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modulation


% --- Executes during object creation, after setting all properties.
function modulation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in awgn.
function awgn_Callback(hObject, eventdata, handles)
% hObject    handle to awgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of awgn
set(handles.rician,'Value',0.0);
set(handles.rayleigh,'Value',0.0);


% --- Executes on button press in rayleigh.
function rayleigh_Callback(hObject, eventdata, handles)
% hObject    handle to rayleigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of rayleigh
set(handles.rician,'Value',0.0);
set(handles.awgn,'Value',0.0);



% --- Executes on button press in rician.
function rician_Callback(hObject, eventdata, handles)
% hObject    handle to rician (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of rician
set(handles.rayleigh,'Value',0.0);
set(handles.awgn,'Value',0.0);


% --- Executes on button press in convolutional.
function convolutional_Callback(hObject, eventdata, handles)
% hObject    handle to convolutional (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of convolutional
set(handles.ldpc,'Value',0.0);




% --- Executes on button press in ldpc.
function ldpc_Callback(hObject, eventdata, handles)
% hObject    handle to ldpc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of ldpc
set(handles.convolutional,'Value',0.0);

% --- Executes on button press in btn_simulation.
function btn_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to btn_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.convolutional,'Value')==get(handles.convolutional,'Max') && get(handles.bpsk,'Value')==get(handles.bpsk,'Max') && get(handles.awgn,'Value')==get(handles.awgn,'Max'))
nbitpersym  = 52;   % number of bits per OFDM symbol
nsym        = str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 64;   % fft size
sub_car     = 52;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(sub_car/len_fft)+ 10*log10(64/80); % symbol to noise ratio
snr=EsNo - 10*log10(52/64); % snr as to be used by awgn fn.
% Generation of DATA
t_data=randi([0 1],nbitpersym*nsym,1);
%Convolution Coder 
M = 4; codeRate = 1/2; constlen = 7; k = log2(M); codegen = [171 133];
tblen = 32;     % traceback length
trellis = poly2trellis(constlen, codegen);
msg_enc = convenc(t_data, trellis);
nsym=nsym/codeRate;
% BPSK Modulation 
hMod = modem.pskmod('M',2,'InputType','Bit');
mod_data = modulate(hMod, msg_enc);
% OFDM Transmitter
% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';
% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;
%Inverse FFT
IFFT_data = (64/sqrt(52))*ifft(fftshift(pilot_ins_data.')).';
% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
% parallel to serial conversion
ser_data = reshape(cylic_add_data,80*nsym,1);
%Sending data through AWGN channel
nChnlErrs = zeros(1,length(EbNo));
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;
for ii=1:length(snr)
chan_awgn = sqrt(80/52)*awgn(ser_data,snr(ii),'measured'); % awgn addition
%Reception OFDM
% serial to parallel conversion
ser_to_para = reshape(chan_awgn,80,nsym).'; 
%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:80]);   
% freq domain transform
FFT_recdata =(sqrt(52)/64)*fftshift(fft(cyclic_pre_rem.')).';   
%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 
% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  
% BPSK Demodulation
hDemod = modem.pskdemod(hMod);
demod_data = demodulate(hDemod,ser_data_1);
%Convolution Decoder 
msg_dec = vitdec(demod_data, trellis, tblen, 'cont', 'hard');
% BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, demod_data);
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data(1:end-tblen), msg_dec(1+tblen:end)); 
end
%Plotting Graphical Results
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With CONVOLUTIONAL Coding');
xlabel('Eb/No');
ylabel('BER');
title('  Modulation:BPSK   Channel:AWGN ');

elseif (get(handles.convolutional,'Value')==get(handles.convolutional,'Max') && get(handles.bpsk,'Value')==get(handles.bpsk,'Max') && get(handles.rayleigh,'Value')==get(handles.rayleigh,'Max'))
nbitpersym  = 52;   % number of bits per OFDM symbol
nsym        = str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 64;   % fft size
sub_car     = 52;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(sub_car/len_fft)+ 10*log10(64/80); % symbol to noise ratio
snr=EsNo - 10*log10(52/64); % snr as to be used by awgn fn.
% Generation of DATA
t_data=randi([0 1],nbitpersym*nsym,1);
%Convolution Coder 
M = 4; codeRate = 1/2; constlen = 7; k = log2(M); codegen = [171 133];
tblen = 32;     % traceback length
trellis = poly2trellis(constlen, codegen);
msg_enc = convenc(t_data, trellis);
nsym=nsym/codeRate;
% BPSK Modulation 
hMod = modem.pskmod('M',2,'InputType','Bit');
mod_data = modulate(hMod, msg_enc);
% OFDM Transmitter
% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';
% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;
%Inverse FFT
IFFT_data = (64/sqrt(52))*ifft(fftshift(pilot_ins_data.')).';
% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
% parallel to serial conversion
ser_data = reshape(cylic_add_data,80*nsym,1);

%Sending data through AWGN channel
nChnlErrs = zeros(1,length(EbNo));
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

rh=rayleighchan(1/10000,1/1000);
changain1=filter(rh,ones(nsym*80,1));
a=max(max(abs(changain1)));
changain1=changain1./a;
chan_data = changain1.*ser_data;

for ii=1:length(snr)
chan_awgn = sqrt(80/52)*awgn(chan_data,snr(ii),'measured'); % awgn addition
chan_awgn =a* chan_awgn./changain1;


%Reception OFDM
% serial to parallel conversion
ser_to_para = reshape(chan_awgn,80,nsym).'; 
%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:80]);   
% freq domain transform
FFT_recdata =(sqrt(52)/64)*fftshift(fft(cyclic_pre_rem.')).';   
%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 
% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  

% BPSK Demodulation
hDemod = modem.pskdemod(hMod);
demod_data = demodulate(hDemod,ser_data_1);  

%Convolution Decoder 
msg_dec = vitdec(demod_data, trellis, tblen, 'cont', 'hard');

% BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, demod_data);
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data(1:end-tblen), msg_dec(1+tblen:end)); 
end
%Plotting Graphical Results
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With CONVOLUTIONAL Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:BPSK   Channel:Rayleigh');

elseif (get(handles.convolutional,'Value')==get(handles.convolutional,'Max') && get(handles.bpsk,'Value')==get(handles.bpsk,'Max') && get(handles.rician,'Value')==get(handles.rician,'Max'))
nbitpersym  = 52;   % number of bits per OFDM symbol
nsym        = str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 64;   % fft size
sub_car     = 52;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(sub_car/len_fft)+ 10*log10(64/80); % symbol to noise ratio
snr=EsNo - 10*log10(52/64); % snr as to be used by awgn fn.

% Generation of DATA

t_data=randi([0 1],nbitpersym*nsym,1);

%Convolution Coder 

M = 4; codeRate = 1/2; constlen = 7; k = log2(M); codegen = [171 133];
tblen = 32;     % traceback length
trellis = poly2trellis(constlen, codegen);


msg_enc = convenc(t_data, trellis);
nsym=nsym/codeRate;

% BPSK Modulation 
hMod = modem.pskmod('M',2,'InputType','Bit');
mod_data = modulate(hMod, msg_enc);

% OFDM Transmitter
% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';

% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;
%Inverse FFT
IFFT_data = (64/sqrt(52))*ifft(fftshift(pilot_ins_data.')).';
% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
% parallel to serial conversion
ser_data = reshape(cylic_add_data,80*nsym,1);

%Sending data through AWGN channel
nChnlErrs = zeros(1,length(EbNo));
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

rh=ricianchan(1/10000,1/800,0);
changain1=filter(rh,ones(nsym*80,1));
a=max(max(abs(changain1)));
changain1=changain1./a;
chan_data = changain1.*ser_data;

for ii=1:length(snr)
chan_awgn = sqrt(80/52)*awgn(chan_data,snr(ii),'measured'); % awgn addition
chan_awgn =a* chan_awgn./changain1;


%Reception OFDM
% serial to parallel conversion
ser_to_para = reshape(chan_awgn,80,nsym).'; 
%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:80]);   
% freq domain transform
FFT_recdata =(sqrt(52)/64)*fftshift(fft(cyclic_pre_rem.')).';   
%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 
% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  

% BPSK Demodulation
hDemod = modem.pskdemod(hMod);
demod_data = demodulate(hDemod,ser_data_1);  

%Convolution Decoder 
msg_dec = vitdec(demod_data, trellis, tblen, 'cont', 'hard');

% BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, demod_data);
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data(1:end-tblen), msg_dec(1+tblen:end)); 
end
%Plotting Graphical Results
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With CONVOLUTIONAL Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:BPSK   Channel:Rician');
elseif (get(handles.convolutional,'Value')==get(handles.convolutional,'Max') && get(handles.qam,'Value')==get(handles.qam,'Max') && get(handles.awgn,'Value')==get(handles.awgn,'Max'))
nbitpersym  = 52;   % number of bits per OFDM symbol
nsym        = str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 64;   % fft size
sub_car     = 52;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(sub_car/len_fft)+ 10*log10(64/80); % symbol to noise ratio
snr=EsNo - 10*log10(52/64); % snr as to be used by awgn fn.

% Generation of DATA

t_data=randi([0 1],nbitpersym*nsym,1);

%Convolution Coder 

M = 4; codeRate = 1/2; constlen = 7; k = log2(M); codegen = [171 133];
tblen = 32;     % traceback length
trellis = poly2trellis(constlen, codegen);
msg_enc = convenc(t_data, trellis);
nsym=nsym/codeRate;

% QAM Modulation 
hMod = modem.qammod('M',2,'InputType','Bit');
mod_data = modulate(hMod, msg_enc);

% OFDM Transmitter
% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';

% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;
%Inverse FFT
IFFT_data = (64/sqrt(52))*ifft(fftshift(pilot_ins_data.')).';
% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
% parallel to serial conversion
ser_data = reshape(cylic_add_data,80*nsym,1);

%Sending data through AWGN channel
nChnlErrs = zeros(1,length(EbNo));
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;


for ii=1:length(snr)
chan_awgn = sqrt(80/52)*awgn(ser_data,snr(ii),'measured'); % awgn addition
%Reception OFDM
% serial to parallel conversion
ser_to_para = reshape(chan_awgn,80,nsym).'; 
%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:80]);   
% freq domain transform
FFT_recdata =(sqrt(52)/64)*fftshift(fft(cyclic_pre_rem.')).';   
%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 
% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  

% BPSK Demodulation
hDemod = modem.qamdemod(hMod);
demod_data = demodulate(hDemod,ser_data_1);


%Convolution Decoder 
msg_dec = vitdec(demod_data, trellis, tblen, 'cont', 'hard');

% BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, demod_data);
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data(1:end-tblen), msg_dec(1+tblen:end)); 
end
%Plotting Graphical Results
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With CONVOLUTIONAL Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:QAM   Channel:AWGN');


elseif (get(handles.convolutional,'Value')==get(handles.convolutional,'Max') && get(handles.qam,'Value')==get(handles.qam,'Max') && get(handles.rayleigh,'Value')==get(handles.rayleigh,'Max'))
nbitpersym  = 52;   % number of bits per OFDM symbol
nsym        = str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 64;   % fft size
sub_car     = 52;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(sub_car/len_fft)+ 10*log10(64/80); % symbol to noise ratio
snr=EsNo - 10*log10(52/64); % snr as to be used by awgn fn.

% Generation of DATA

t_data=randi([0 1],nbitpersym*nsym,1);

%Convolution Coder 

M = 4; codeRate = 1/2; constlen = 7; k = log2(M); codegen = [171 133];
tblen = 32;     % traceback length
trellis = poly2trellis(constlen, codegen);


msg_enc = convenc(t_data, trellis);
nsym=nsym/codeRate;

% QAM Modulation 
hMod = modem.qammod('M',2,'InputType','Bit');
mod_data = modulate(hMod, msg_enc);

% OFDM Transmitter
% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';

% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;
%Inverse FFT
IFFT_data = (64/sqrt(52))*ifft(fftshift(pilot_ins_data.')).';
% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
% parallel to serial conversion
ser_data = reshape(cylic_add_data,80*nsym,1);

%Sending data through channel
nChnlErrs = zeros(1,length(EbNo));
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

rh=rayleighchan(1/10000,1/1000);
changain1=filter(rh,ones(nsym*80,1));
a=max(max(abs(changain1)));
changain1=changain1./a;
chan_data = changain1.*ser_data;

for ii=1:length(snr)
chan_awgn = sqrt(80/52)*awgn(chan_data,snr(ii),'measured'); % awgn addition
chan_awgn =a* chan_awgn./changain1;


%Reception OFDM
% serial to parallel conversion
ser_to_para = reshape(chan_awgn,80,nsym).'; 
%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:80]);   
% freq domain transform
FFT_recdata =(sqrt(52)/64)*fftshift(fft(cyclic_pre_rem.')).';   
%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 
% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  

% QAM Demodulation
hDemod = modem.qamdemod(hMod);
demod_data = demodulate(hDemod,ser_data_1);  

%Convolution Decoder 
msg_dec = vitdec(demod_data, trellis, tblen, 'cont', 'hard');

% BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, demod_data);
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data(1:end-tblen), msg_dec(1+tblen:end)); 
end
%Plotting Graphical Results
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With CONVOLUTIONAL Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:QAM  Channel:Rayleigh');


elseif (get(handles.convolutional,'Value')==get(handles.convolutional,'Max') && get(handles.qam,'Value')==get(handles.qam,'Max') && get(handles.rician,'Value')==get(handles.rician,'Max'))
nbitpersym  = 52;   % number of bits per OFDM symbol
nsym        = str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 64;   % fft size
sub_car     = 52;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(sub_car/len_fft)+ 10*log10(64/80); % symbol to noise ratio
snr=EsNo - 10*log10(52/64); % snr as to be used by awgn fn.

% Generation of DATA

t_data=randi([0 1],nbitpersym*nsym,1);

%Convolution Coder 

M = 4; codeRate = 1/2; constlen = 7; k = log2(M); codegen = [171 133];
tblen = 32;     % traceback length
trellis = poly2trellis(constlen, codegen);


msg_enc = convenc(t_data, trellis);
nsym=nsym/codeRate;

% BPSK Modulation 
hMod = modem.pskmod('M',2,'InputType','Bit');
mod_data = modulate(hMod, msg_enc);

% OFDM Transmitter
% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';

% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;
%Inverse FFT
IFFT_data = (64/sqrt(52))*ifft(fftshift(pilot_ins_data.')).';
% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
% parallel to serial conversion
ser_data = reshape(cylic_add_data,80*nsym,1);

%Sending data through AWGN channel
nChnlErrs = zeros(1,length(EbNo));
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

rh=ricianchan(1/10000,1/800,0);
changain1=filter(rh,ones(nsym*80,1));
a=max(max(abs(changain1)));
changain1=changain1./a;
chan_data = changain1.*ser_data;

for ii=1:length(snr)
chan_awgn = sqrt(80/52)*awgn(chan_data,snr(ii),'measured'); % awgn addition
chan_awgn =a* chan_awgn./changain1;


%Reception OFDM
% serial to parallel conversion
ser_to_para = reshape(chan_awgn,80,nsym).'; 
%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:80]);   
% freq domain transform
FFT_recdata =(sqrt(52)/64)*fftshift(fft(cyclic_pre_rem.')).';   
%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 
% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  

% BPSK Demodulation
hDemod = modem.pskdemod(hMod);
demod_data = demodulate(hDemod,ser_data_1);  

%Convolution Decoder 
msg_dec = vitdec(demod_data, trellis, tblen, 'cont', 'hard');

% BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, demod_data);
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data(1:end-tblen), msg_dec(1+tblen:end)); 
end
%Plotting Graphical Results
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With CONVOLUTIONAL Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:QAM   Channel:Rician');


elseif (get(handles.ldpc,'Value')==get(handles.ldpc,'Max') && get(handles.bpsk,'Value')==get(handles.bpsk,'Max') && get(handles.awgn,'Value')==get(handles.awgn,'Max'))
nbitpersym  = 324;   % number of bits per OFDM symbol
nsym        = str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 336;   % fft size
sub_car     = 648;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(nbitpersym/len_fft)+ 10*log10(len_fft/352); % symbol to noise ratio
SNR=EsNo - 10*log10(nbitpersym/len_fft); % snr as to be used by awgn fn.

Nsamp = 1;

% Generation of DATA
t_data=randi([0 1],1,nbitpersym*nsym); %enc.NumInfoBits

%LDPC Coding
codeRate=1/2;
enc = fec.ldpcenc;
dec = fec.ldpcdec;
dec.DecisionType = 'Hard decision';
dec.OutputFormat = 'Information part';
dec.NumIterations = 20;

% Stop if all parity-checks are satisfied
dec.DoParityChecks = 'Yes';
msg_enc = encode(enc,t_data); 
paritychecks = mod(enc.ParityCheckMatrix * msg_enc', 2);
nsym=nsym/codeRate;


% %  BPSK Modulation
hMod = modem.pskmod('M', 2, 'PhaseOffset', 0, 'SymbolOrder',...
    'binary', 'InputType', 'Bit');
% Modulate the signal (map bit 0 to 1 + 0i, bit 1 to -1 + 0i)
mod_data = modulate(hMod, msg_enc);


%OFDM Transmitter

% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';

% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;

% Inverse FFT
IFFT_data = (len_fft/sqrt(nbitpersym))*ifft(fftshift(pilot_ins_data.')).';

% Adding cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';

% parallel to serial conversion
ser_data = reshape(cylic_add_data,352*nsym,1);

%Sending data through AWGN channel
lg=length(SNR);
nChnlErrs = zeros(1,lg);
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;
for ii=1:length(SNR)
chan_awgn = sqrt(352/nbitpersym)*awgn(ser_data,SNR(ii),'measured'); % awgn addition


%Reception OFDM

% serial to parallel coversion
ser_to_para = reshape(chan_awgn,352,nsym).'; 

%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:352]);   

% freq domain transform
FFT_recdata =(sqrt(nbitpersym)/len_fft)*fftshift(fft(cyclic_pre_rem.')).';   

%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 

% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  



%BPSK Demodulation
hDemod1 = modem.pskdemod(hMod,'DecisionType','llr', ...
     'NoiseVariance',3.91);
hDemod2 = modem.pskdemod(hMod);
sortieHard = demodulate(hDemod2,ser_data_1);
llr = demodulate(hDemod1, ser_data_1);
demod_data=llr';

%LDPC Decoder 
msg_dec = decode(dec, demod_data);

%BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, sortieHard');
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data, msg_dec);
end



%Results on Graph
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With LDPC Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:BPSK    Channel:AWGN');

elseif (get(handles.ldpc,'Value')==get(handles.ldpc,'Max') && get(handles.bpsk,'Value')==get(handles.bpsk,'Max') && get(handles.rayleigh,'Value')==get(handles.rayleigh,'Max'))

nbitpersym  = 324;   % number of bits per OFDM symbol
nsym        = str2double(get(handles.nsym,'String'));% number of symbols
len_fft     = 336;   % fft size
sub_car     = 648;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(324/336)+ 10*log10(336/352); % symbol to noise ratio
SNR=EsNo - 10*log10(324/336); % snr as to be used by awgn fn.

Nsamp = 1;

% Generation of DATA
t_data=randi([0 1],1,nbitpersym*nsym); %enc.NumInfoBits

%LDPC Coding
codeRate=1/2;
enc = fec.ldpcenc;
dec = fec.ldpcdec;
dec.DecisionType = 'Hard decision';
dec.OutputFormat = 'Information part';
dec.NumIterations = 20;

% Stop if all parity-checks are satisfied
dec.DoParityChecks = 'Yes';

%n = 64800; k = 32400;
msg_enc = encode(enc,t_data); %On doit avoir : size(t_data)=enc.NumInfoBits
paritychecks = mod(enc.ParityCheckMatrix * msg_enc', 2);
nsym=nsym/codeRate;


% %  BPSK Modulation
hMod = modem.pskmod('M', 2, 'PhaseOffset', 0, 'SymbolOrder',...
    'binary', 'InputType', 'Bit');
% Modulate the signal (map bit 0 to 1 + 0i, bit 1 to -1 + 0i)
mod_data = modulate(hMod, msg_enc);


%OFDM Transmitter

% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';

% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;

% Inverse FFT
IFFT_data = (336/sqrt(324))*ifft(fftshift(pilot_ins_data.')).';

% Adding cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';

% parallel to serial conversion
ser_data = reshape(cylic_add_data,352*nsym,1);

%Sending data through AWGN channel
lg=length(SNR);
nChnlErrs = zeros(1,lg);
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

rh=rayleighchan(1/10000,1/1000);
changain1=filter(rh,ones(nsym*352,1));
a=max(max(abs(changain1)));
changain1=changain1./a;
chan_data = changain1.*ser_data;

for ii=1:length(SNR)
chan_awgn = sqrt(352/324)*awgn(chan_data,SNR(ii),'measured'); % awgn addition
chan_awgn =a* chan_awgn./changain1;


%Reception OFDM

% serial to parallel coversion
ser_to_para = reshape(chan_awgn,352,nsym).'; 

%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:352]);   

% freq domain transform
FFT_recdata =(sqrt(324)/336)*fftshift(fft(cyclic_pre_rem.')).';   

%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 

% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  



%BPSK Demodulation
hDemod1 = modem.pskdemod(hMod,'DecisionType','llr', ...
     'NoiseVariance',3.91);
hDemod2 = modem.pskdemod(hMod);
sortieHard = demodulate(hDemod2,ser_data_1);
llr = demodulate(hDemod1, ser_data_1);
demod_data=llr';

%LDPC Decoder 
msg_dec = decode(dec, demod_data);

%BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, sortieHard');
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data, msg_dec);
end
%Results on Graph
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With LDPC Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:BPSK    Channel:Rayleigh');

elseif (get(handles.ldpc,'Value')==get(handles.ldpc,'Max') && get(handles.bpsk,'Value')==get(handles.bpsk,'Max') && get(handles.rician,'Value')==get(handles.rician,'Max'))
nbitpersym  = 324;   % number of bits per OFDM symbol
nsym        =  str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 336;   % fft size
sub_car     = 648;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(324/336)+ 10*log10(336/352); % symbol to noise ratio
SNR=EsNo - 10*log10(324/336); % snr as to be used by awgn fn.

Nsamp = 1;

% Generation of DATA
t_data=randi([0 1],1,nbitpersym*nsym); %enc.NumInfoBits

%LDPC Coding
codeRate=1/2;
enc = fec.ldpcenc;
dec = fec.ldpcdec;
dec.DecisionType = 'Hard decision';
dec.OutputFormat = 'Information part';
dec.NumIterations = 200;

% Stop if all parity-checks are satisfied
dec.DoParityChecks = 'Yes';

%n = 64800; k = 32400;
msg_enc = encode(enc,t_data); %On doit avoir : size(t_data)=enc.NumInfoBits
paritychecks = mod(enc.ParityCheckMatrix * msg_enc', 2);
nsym=nsym/codeRate;


% %  BPSK Modulation
hMod = modem.qammod('M', 2, 'PhaseOffset', 0, 'SymbolOrder',...
    'binary', 'InputType', 'Bit');
% Modulate the signal (map bit 0 to 1 + 0i, bit 1 to -1 + 0i)
mod_data = modulate(hMod, msg_enc);


%OFDM Transmitter

% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';

% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;

% Inverse FFT
IFFT_data = (336/sqrt(324))*ifft(fftshift(pilot_ins_data.')).';

% Adding cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';

% parallel to serial conversion
ser_data = reshape(cylic_add_data,352*nsym,1);

%Sending data through AWGN channel
lg=length(SNR);
nChnlErrs = zeros(1,lg);
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

rh=ricianchan(1/10000,1/800,0);
changain1=filter(rh,ones(nsym*352,1));
a=max(max(abs(changain1)));
changain1=changain1./a;
chan_data = changain1.*ser_data;

for ii=1:length(SNR)
chan_awgn = sqrt(352/324)*awgn(chan_data,SNR(ii),'measured'); % awgn addition
chan_awgn =a* chan_awgn./changain1;


%Reception OFDM

% serial to parallel coversion
ser_to_para = reshape(chan_awgn,352,nsym).'; 

%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:352]);   

% freq domain transform
FFT_recdata =(sqrt(324)/336)*fftshift(fft(cyclic_pre_rem.')).';   

%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 

% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  



%BPSK Demodulation
hDemod1 = modem.qamdemod(hMod,'DecisionType','llr', ...
     'NoiseVariance',3.91);
hDemod2 = modem.qamdemod(hMod);
sortieHard = demodulate(hDemod2,ser_data_1);
llr = demodulate(hDemod1, ser_data_1);
demod_data=llr';

%LDPC Decoder 
msg_dec = decode(dec, demod_data);

%BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, sortieHard');
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data, msg_dec);
end



%Results on Graph
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With LDPC Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:BPSK    Channel:Rician');



elseif (get(handles.ldpc,'Value')==get(handles.ldpc,'Max') && get(handles.qam,'Value')==get(handles.qam,'Max') && get(handles.awgn,'Value')==get(handles.awgn,'Max'))
nbitpersym  = 324;   % number of bits per OFDM symbol
nsym        = str2double(get(handles.nsym,'String'));% number of symbols
len_fft     = 336;   % fft size
sub_car     = 648;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(324/336)+ 10*log10(336/352); % symbol to noise ratio
SNR=EsNo - 10*log10(324/336); % snr as to be used by awgn fn.

%SNR = -19:5:16;
%SNR = 5:2:20;
Nsamp = 1;

% Generation of DATA
t_data=randi([0 1],1,nbitpersym*nsym); %enc.NumInfoBits

%Coder LDPC
codeRate=1/2;
enc = fec.ldpcenc;
dec = fec.ldpcdec;
dec.DecisionType = 'Hard decision';
dec.OutputFormat = 'Information part';
dec.NumIterations = 20;
% Stop if all parity-checks are satisfied
dec.DoParityChecks = 'Yes';

%n = 64800; k = 32400;
msg_enc = encode(enc,t_data); %On doit avoir : size(t_data)=enc.NumInfoBits
paritychecks = mod(enc.ParityCheckMatrix * msg_enc', 2);
nsym=nsym/codeRate;

% %  BPSK Modulation
hMod = modem.qammod('M', 2, 'PhaseOffset', 0, 'SymbolOrder','binary', 'InputType', 'Bit');
% Modulate the signal (map bit 0 to 1 + 0i, bit 1 to -1 + 0i)
mod_data = modulate(hMod, msg_enc);

%Emission OFDM
% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';
% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;
% fourier transform time doamain data and normalizing the data
IFFT_data = (336/sqrt(324))*ifft(fftshift(pilot_ins_data.')).';
% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
% parallel to serial coversion

ser_data = reshape(cylic_add_data,352*nsym,1);

%Sending data through AWGN channel
lg=length(SNR);
nChnlErrs = zeros(1,lg);
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

for ii=1:length(SNR)
chan_awgn = sqrt(352/324)*awgn(ser_data,SNR(ii),'measured'); % awgn addition
%Reception OFDM
% serial to parallel coversion
ser_to_para = reshape(chan_awgn,352,nsym).'; 
%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:352]);   
% freq domain transform
FFT_recdata =(sqrt(324)/336)*fftshift(fft(cyclic_pre_rem.')).';   
%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 
% serial coversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  


hDemod1 = modem.qamdemod(hMod,'DecisionType','llr', ...
     'NoiseVariance',3.91);
hDemod2 = modem.qamdemod(hMod);
sortieHard = demodulate(hDemod2,ser_data_1);
llr = demodulate(hDemod1, ser_data_1);
demod_data=llr';

%Decoder LDPC
msg_dec = decode(dec, demod_data);

% BER
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, sortieHard');
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data, msg_dec);
end
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With LDPC Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:QAM   Channel:AWGN');
elseif (get(handles.ldpc,'Value')==get(handles.ldpc,'Max') && get(handles.qam,'Value')==get(handles.qam,'Max') && get(handles.rayleigh,'Value')==get(handles.rayleigh,'Max'))
nbitpersym  = 324;   % number of bits per OFDM symbol
nsym= str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 336;   % fft size
sub_car     = 648;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(324/336)+ 10*log10(336/352); % symbol to noise ratio
SNR=EsNo - 10*log10(324/336); % snr as to be used by awgn fn.
Nsamp = 1;
% Generation of DATA
t_data=randi([0 1],1,nbitpersym*nsym); %enc.NumInfoBits

%LDPC Coding
codeRate=1/2;
enc = fec.ldpcenc;
dec = fec.ldpcdec;
dec.DecisionType = 'Hard decision';
dec.OutputFormat = 'Information part';
dec.NumIterations = 20;

% Stop if all parity-checks are satisfied
dec.DoParityChecks = 'Yes';

%n = 64800; k = 32400;
msg_enc = encode(enc,t_data); %On doit avoir : size(t_data)=enc.NumInfoBits
paritychecks = mod(enc.ParityCheckMatrix * msg_enc', 2);
nsym=nsym/codeRate;


% %  QAM Modulation
hMod = modem.qammod('M', 2, 'PhaseOffset', 0, 'SymbolOrder',...
    'binary', 'InputType', 'Bit');
% Modulate the signal (map bit 0 to 1 + 0i, bit 1 to -1 + 0i)
mod_data = modulate(hMod, msg_enc);


%OFDM Transmitter

% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';

% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;

% Inverse FFT
IFFT_data = (336/sqrt(324))*ifft(fftshift(pilot_ins_data.')).';

% Adding cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';

% parallel to serial conversion
ser_data = reshape(cylic_add_data,352*nsym,1);

%Sending data through AWGN channel
lg=length(SNR);
nChnlErrs = zeros(1,lg);
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

rh=rayleighchan(1/10000,1/1000);
changain1=filter(rh,ones(nsym*352,1));
a=max(max(abs(changain1)));
changain1=changain1./a;
chan_data = changain1.*ser_data;

for ii=1:length(SNR)
chan_awgn = sqrt(352/324)*awgn(chan_data,SNR(ii),'measured'); % awgn addition
chan_awgn =a* chan_awgn./changain1;


%Reception OFDM

% serial to parallel coversion
ser_to_para = reshape(chan_awgn,352,nsym).'; 

%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:352]);   

% freq domain transform
FFT_recdata =(sqrt(324)/336)*fftshift(fft(cyclic_pre_rem.')).';   

%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 

% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  



%BPSK Demodulation
hDemod1 = modem.qamdemod(hMod,'DecisionType','llr', ...
     'NoiseVariance',3.91);
hDemod2 = modem.qamdemod(hMod);
sortieHard = demodulate(hDemod2,ser_data_1);
llr = demodulate(hDemod1, ser_data_1);
demod_data=llr';

%LDPC Decoder 
msg_dec = decode(dec, demod_data);

%BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, sortieHard');
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data, msg_dec);
end



%Results on Graph
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With LDPC Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:QAM    Channel:Rayleigh');

elseif (get(handles.ldpc,'Value')==get(handles.ldpc,'Max') && get(handles.qam,'Value')==get(handles.qam,'Max') && get(handles.rician,'Value')==get(handles.rician,'Max'))
nbitpersym  = 324;   % number of bits per OFDM symbol
nsym= str2double(get(handles.nsym,'String')); % number of symbols
len_fft     = 336;   % fft size
sub_car     = 648;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(324/336)+ 10*log10(336/352); % symbol to noise ratio
SNR=EsNo - 10*log10(324/336); % snr as to be used by awgn fn.
Nsamp = 1;

% Generation of DATA
t_data=randi([0 1],1,nbitpersym*nsym); %enc.NumInfoBits

%LDPC Coding
codeRate=1/2;
enc = fec.ldpcenc;
dec = fec.ldpcdec;
dec.DecisionType = 'Hard decision';
dec.OutputFormat = 'Information part';
dec.NumIterations = 20;

% Stop if all parity-checks are satisfied
dec.DoParityChecks = 'Yes';

%n = 64800; k = 32400;
msg_enc = encode(enc,t_data); %On doit avoir : size(t_data)=enc.NumInfoBits
paritychecks = mod(enc.ParityCheckMatrix * msg_enc', 2);
nsym=nsym/codeRate;


% %  BPSK Modulation
hMod = modem.qammod('M', 2, 'PhaseOffset', 0, 'SymbolOrder',...
    'binary', 'InputType', 'Bit');
% Modulate the signal (map bit 0 to 1 + 0i, bit 1 to -1 + 0i)
mod_data = modulate(hMod, msg_enc);
%OFDM Transmitter
% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';
% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;
% Inverse FFT
IFFT_data = (336/sqrt(324))*ifft(fftshift(pilot_ins_data.')).';
% Adding cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
% parallel to serial conversion
ser_data = reshape(cylic_add_data,352*nsym,1);
%Sending data through AWGN channel
lg=length(SNR);
nChnlErrs = zeros(1,lg);
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

rh=ricianchan(1/10000,1/800,0);
changain1=filter(rh,ones(nsym*352,1));
a=max(max(abs(changain1)));
changain1=changain1./a;
chan_data = changain1.*ser_data;

for ii=1:length(SNR)
chan_awgn = sqrt(352/324)*awgn(chan_data,SNR(ii),'measured'); % awgn addition
chan_awgn =a* chan_awgn./changain1;


%Reception OFDM

% serial to parallel coversion
ser_to_para = reshape(chan_awgn,352,nsym).'; 

%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:352]);   

% freq domain transform
FFT_recdata =(sqrt(324)/336)*fftshift(fft(cyclic_pre_rem.')).';   

%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]); 

% serial conversion
ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);  



%BPSK Demodulation
hDemod1 = modem.qamdemod(hMod,'DecisionType','llr', ...
     'NoiseVariance',3.91);
hDemod2 = modem.qamdemod(hMod);
sortieHard = demodulate(hDemod2,ser_data_1);
llr = demodulate(hDemod1, ser_data_1);
demod_data=llr';

%LDPC Decoder 
msg_dec = decode(dec, demod_data);

%BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, sortieHard');
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data, msg_dec);
end
%Results on Graph
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('Without Coding','With LDPC Coding');
xlabel('Eb/No');
ylabel('BER');
title('Modulation:QAM   Channel:Rician');
end
% --- Executes on button press in bpsk.
function bpsk_Callback(hObject, eventdata, handles)
% hObject    handle to bpsk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of bpsk
set(handles.qam,'Value',0.0);


% --- Executes on button press in qam.
function qam_Callback(hObject, eventdata, handles)
% hObject    handle to qam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of qam
set(handles.bpsk,'Value',0.0);



% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
