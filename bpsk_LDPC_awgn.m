close all
clear all
clc

nbitpersym  = 324;   % number of bits per OFDM symbol
nsym        = 10^2; % number of symbols
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

%n = sub_car00; k = nbitpersym00;
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
legend('OFDM without Coding','OFDM with LDPC Coding');
xlabel('Eb/No');
ylabel('BER');
title('Bit Error Rate of OFDM  signal with BPSK modulation');
