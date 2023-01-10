close all
clear all
clc

nbitpersym  = 324;   % number of bits per OFDM symbol
nsym        = 10^2; % number of symbols
len_fft     = 336;   % fft size
sub_car     = 567;   % number of data subcarriers
EbNo        = -2:1:10;
EsNo= EbNo + 10*log10(324/336)+ 10*log10(336/352); % symbol to noise ratio
SNR=EsNo - 10*log10(324/336); % snr as to be used by awgn fn.

%SNR = -19:5:16;
%SNR = 0:2:15;
Nsamp = 1;

%Codeur BCH
n = 7; k =4;
%enc = fec.bchenc(n,k);
[gp,t] = bchgenpoly(n,k); % t is error-correction capability.
nmc1=nbitpersym*nsym/k;
%dec = fec.bchdec(enc);

% Generation of DATA
%t_data=randi([0 1],1,nbitpersym*nsym);
%t_data = gf(randi([0 1],nsym,k)); % Random k-symbol messages
%t_data = randi([0 1],nmc1,k);
t_data = gf(randi([0 1],nmc1,k));

%msg_enc = encode(enc,t_data); %
msg_enc = bchenc(t_data,n,k); % Encode the data.
codeRate = k/n;
nsym=nsym/codeRate;
%nmc2=nbitpersym*nsym/k;

% % Modulation BPSK
hMod = modem.pskmod('M',2,'InputType','Bit');
% Modulate the signal (map bit 0 to 1 + 0i, bit 1 to -1 + 0i)
msg_enco = double(msg_enc.x);
mod_data = modulate(hMod, msg_enco);

%Emission OFDM
% serial to parallel conversion
par_data = reshape(mod_data,nbitpersym,nsym).';
% pilot insertion
pilot_ins_data=[zeros(nsym,6) par_data(:,[1:nbitpersym/2]) zeros(nsym,1) par_data(:,[nbitpersym/2+1:nbitpersym]) zeros(nsym,5)] ;
% fourier transform time doamain data and normalizing the data
%IFFT_data = (64/sqrt(324))*ifft(fftshift(pilot_ins_data.')).';
IFFT_data = (336/sqrt(324))*ifft(pilot_ins_data.').';



% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[49:64]) IFFT_data].';
% parallel to serial conversion
%ser_data = reshape(cylic_add_data,80*nsym,1);
%ser_data = reshape(cylic_add_data,352*nsym,1);
ser_data = reshape(cylic_add_data,352*nsym,1);

%Transmission sur le canal AWGN
nChnlErrs = zeros(1,length(EbNo));
BERChnl = nChnlErrs;
nCodErrs = nChnlErrs;
BERCoded = nChnlErrs;

for ii=1:length(SNR)
%sigma = sqrt(10^(SNR(ii)/10));
%chan_awgn = sqrt(80/52)*awgn(ser_data,sigma,'measured'); % awgn addition
%chan_awgn = awgn(ser_data,sigma,'measured'); % awgn addition

chan_awgn = sqrt(352/324)*awgn(ser_data,SNR(ii),'measured'); % awgn addition

%Reception OFDM
% serial to parallel coversion
ser_to_para = reshape(chan_awgn,352,nsym).'; 

%cyclic prefix removal
cyclic_pre_rem = ser_to_para(:,[17:352]); 
%cyclic_pre_rem = IFFT_data;

% freq domain transform
%FFT_recdata =(sqrt(324)/64)*fftshift(fft(cyclic_pre_rem.')).'; 
FFT_recdata =(sqrt(324)/336)*fft(cyclic_pre_rem.').';
%pilot removal
rem_pilot = FFT_recdata (:,[6+[1:nbitpersym/2] 7+[nbitpersym/2+1:nbitpersym] ]);
% serial conversion
%ser_data_1 = reshape(rem_pilot.',nbitpersym*nsym,1);
ser_data_1 = reshape(rem_pilot.',nmc1,n);

% Visualize received signal
%pp=reshape(ser_data_1,1,nmc1*n);
%scatterplot(pp);
 % scatterplot(ser_data_1)

%msg_rx_int = intdump(ser_data_1, Nsamp);

%demodulation BPSK
hDemod = modem.pskdemod(hMod);
demod_data = demodulate(hDemod,ser_data_1);


%Decodeur BCH

msg_dec = bchdec(gf(demod_data), n, k);


% Calcule du BER
%t_data_ok = double(t_data.x);
%msg_dec_ok = double(msg_dec.x);
[nChnlErrs(ii) BERChnl(ii)] = biterr(msg_enco, demod_data);
[nCodErrs(ii) BERCoded(ii)] = biterr(double(t_data.x), double(msg_dec.x));
end

% Affichage des résultats sur graphique
figure
semilogy(EbNo,BERChnl,'--*b','linewidth',2);
hold on;
semilogy(EbNo,BERCoded,'ro--','linewidth',2);
%semilogy(EbNo,ratio,'--or','linewidth',2);
grid on
%axis([-19 16 10^-5 .1]);
axis([-2 10 10^-5 .1]);
legend('OFDM without Coding','OFDM with BCH Coding');
xlabel('Eb/No');
ylabel('BER');
title('Bit Error Rate of OFDM  signal with BPSK modulation');
