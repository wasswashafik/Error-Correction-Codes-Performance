close all
clear all
clc

nbitpersym  = 52;   % number of bits per OFDM symbol
nsym        = 10^4; % number of symbols
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
%demod_data = demodulate(modem.pskdemod(2),ser_data_1);  

%Convolution Decoder 
msg_dec = vitdec(demod_data, trellis, tblen, 'cont', 'hard');

% BER Calculation
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enc, demod_data);
[nCodErrs(ii) BERCoded(ii)] = biterr(t_data(1:end-tblen), msg_dec(1+tblen:end)); 
end
%Plotting Graphical Results
BERwithCode  = mean(BERCoded)
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
grid on
axis([-2 10 10^-5 .1]);
legend('OFDM without Coding','OFDM with CONVOLUTIONAL Coding');
xlabel('Eb/No');
ylabel('BER');
title('Bit Error Rate of OFDM  signal with BPSK modulation');
