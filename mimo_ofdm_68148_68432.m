function trab2_mimo_ofdm_68148_68432
%CSF 18/19
% Daniel Cunha 68148 & Tiago Marques 68432
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference ofdm chain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
% allocating memory & Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_OFDM_SYM =1e3;
Nc=768; % exemple, number of subcarriers
EbN0=0:2:20 ;
load ('pdp.mat','pdp');
fftSize=1024;
fa=15.36e6;
cp_time=5.21e-6;
m=4;
PR=1;
N_bits = Nc*2;
ber_freq21 = zeros(1,length(N_OFDM_SYM ));
ber_freq22 = ber_freq21;
ber_each_freq21 = zeros(1,length(EbN0));
ber_each_freq22= ber_each_freq21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:length(EbN0)    % EbN0_dB range
    
    for k=1:N_OFDM_SYM   % number of OFDM symbols to be simulated
        
        %################# Transmitter #################
        % data generation
        data = data_gen(N_bits);
        %data modulation
        data_mod = mod_data(data,m);
        % Space-Frequency Coding
        [ant1,ant2]= sf_coding(Nc,data_mod);
     
        %################# Channel Model #################
        %awgn noise
        var = PR*10.^(-EbN0(p)/10)/log2(m);
        Noise = sqrt(var/2).*(randn(1,Nc) + 1j*randn(1,Nc));
        %multipath channel
        [~, H11f] = channel_gen(pdp, fa, fftSize);
        [~, H21f] = channel_gen(pdp, fa, fftSize);
        [~, H12f] = channel_gen(pdp, fa, fftSize);
        [~, H22f] = channel_gen(pdp, fa, fftSize);
        %################# Receiver ######################
        sf_noise= ant1.*H11f(129:Nc+128) + ant2.*H21f(129:Nc+128)+Noise;
        sf_noise2= ant1.*H12f(129:Nc+128) + ant2.*H22f(129:Nc+128)+Noise;
        
        % Space-Frequency Decoding
        s_21 = sf_decoding2x1(sf_noise,H11f(129:Nc+128),H21f(129:Nc+128));
        s_22 = sf_decoding2x2(sf_noise,sf_noise2,H11f(129:Nc+128),H12f(129:Nc+128),H21f(129:Nc+128),H22f(129:Nc+128));
        % Data demodulation
        SF_21 = demod_data(s_21,m,N_bits);
        SF_22 = demod_data(s_22,m,N_bits);
        % BER computation
        ber_freq21(k) = compute_ber(SF_21,data);
        ber_freq22(k) = compute_ber(SF_22,data);
    end
    ber_each_freq21(p)= mean(ber_freq21);
    ber_each_freq22(p)= mean(ber_freq22);
end
[ber_each_eb,~,ber_each_ebf]=trab1(EbN0);
figure(1)
semilogy(EbN0,ber_each_freq21,'r',EbN0,ber_each_eb,'k');
axis([0 20 10^-5 0.5]),grid on,xlabel('EbN0 (dB)'),ylabel('BER'),
legend('2x1 MISO','awgn');
figure(2)
semilogy(EbN0,ber_each_freq21,'r',EbN0,ber_each_freq22,'g',EbN0,ber_each_eb,'k',EbN0,ber_each_ebf,'y');
axis([0 20 10^-5 0.5]),grid on,xlabel('EbN0 (dB)'),ylabel('BER'),
legend('2x1 MISO','2x2 MIMO','awgn','freq ofdm');

end
%%%%%%%%%%%%%%% Auxiliary Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************** Data Generation *********************
function [randbits] = data_gen(N_bits)
randbits = randi([0 1],1,N_bits);
end

function data_symbol=mod_data(data, m)
% Data Modulation
%BPSK, m=2
%QPSK, m=4,
%16-QAM, m=16

switch m
    case 2
        data_symbol = -data*2+1; %mapping of 1->-1 and 0->1
        
    case 4
        % Coding of data bits in QPSK symbols - using Grey coding
        % (00->1+i; 01->1-i; 10->-1+i; 11->-1-i)
        % bit MS defines real polarity
        % bit LS defines imag polarity
        data_temp = reshape(data,2,length(data)/2);
        data_real = data_temp(1,:);
        data_imag = data_temp(2,:);
        data_symbol = sqrt(2)/2*((-1).^(data_real)+1i*(-1).^(data_imag));
        
    case 16
        data_temp = reshape(data,4,length(data)/4);
        data_r1 = data_temp(1,:);
        data_i1 = data_temp(2,:);
        data_r2 = data_temp(3,:);
        data_i2 = data_temp(4,:);
        data_symbol = 2/sqrt(10).*(0.5*(-1).^(data_r2).*(-1).^(data_r1)+(-1).^( data_r1)+1i.*(0.5*(-1).^(data_i2).*(-1).^(data_i1)+ (-1).^(data_i1)));
    otherwise
        helpdlg('Constellation size (m) not available');
end
end

function decoded_data=demod_data(data_symbol, m, N_Data)

vect_IMAG = imag(data_symbol);
vect_REAL = real(data_symbol);
coder_type_value=0;
switch m
    case 2
        if (coder_type_value==0)
            %hard decision
            decoded_data= ceil(-vect_REAL./(1.0000000000001.*max(abs(vect_REAL))));
            
        else
            %soft decision
            decoded_data = vect_REAL;
        end
        
    case 4
        % Decoding of data bits in QPSK symbols - using Grey coding
        % (1+i->00; 1-i->01; -1+i->10; -1-i->11)
        % real polarity defines bit MS
        % imag polarity defines bit LS
        if (coder_type_value==0)
            %hard decision
            vect_REAL_1 = ceil(-vect_REAL./(1.0000000000001.*max(abs(vect_REAL))));
            vect_IMAG_1 = ceil(-vect_IMAG./(1.0000000000001.*max(abs(vect_IMAG))));
            decoded_data = reshape([vect_REAL_1; vect_IMAG_1],1,N_Data);
        else
            %soft decision
            decoded_data = reshape([vect_REAL; vect_IMAG],1,N_Data);
        end
    case 16
        P_1= vect_REAL;
        P_2= vect_IMAG;
        P_3= abs(vect_REAL)-2/sqrt(10);
        P_4= abs(vect_IMAG)-2/sqrt(10);
        if (coder_type_value==0)
            %hard decision
            vect_IMAG_1 = ceil(-P_2./(1.0000000000001.*max(abs(P_2))));
            vect_IMAG_2 = ceil(-P_4./(1.0000000000001.*max(abs(P_4))));
            vect_REAL_1 = ceil(-P_1./(1.0000000000001.*max(abs(P_1))));
            vect_REAL_2 = ceil(-P_3./(1.0000000000001.*max(abs(P_3))));
            decoded_data =  reshape([vect_REAL_1; vect_IMAG_1; vect_REAL_2; vect_IMAG_2],1,N_Data);
        else
            %soft decision
            decoded_data =  reshape([P_1; P_2; P_3; P_4],1,N_Data);
        end
        
    otherwise
        helpdlg('Constellation size (m) not available');
end
end

function y=conv_s_h(s,h,pdp,Nc,samp_freq,tg)
%tg=5.21e-6;                       % guerad time -> 80 samples

delays=pdp(:,1);
delta_t=1/samp_freq;
Npaths = length(delays);                      % No. of paths considered for the channel
deltans=delta_t/1e-9;             % Sampling interval in ns


delays = round(delays/(deltans))+1;
Ng= round(tg/(delta_t))+1;
f_zeros=find(h==0);
h(f_zeros)=[];

aux = zeros(Npaths,Nc+Ng);

for n=1:Npaths
    
    
    conv_sh=h(n)*s;
    
    aux(n,delays(n):Nc+Ng-1+delays(n)-1)=conv_sh;
    
end


y=sum(aux);
y=y(1:Nc+Ng-1);

end

function [ber]=compute_ber(x,y)
ber = sum(x~=y)/length(x);
end

function [ht, Hf]=channel_gen(pdp,samp_freq, Nc)

delta_t=1/samp_freq;              %sample duratiion
Npaths = length(pdp(:,1));        % No. of paths considered for the channel
deltans=delta_t/1e-9;             % Sampling interval in ns

path_pot_lin=10.^(pdp(:,2)/10);
path_pot_lin=path_pot_lin./sum(path_pot_lin);

delays = pdp(:,1);
delays = round(delays./(deltans))+1;

multipath = zeros(1,Npaths);

for n=1:Npaths
    pot_componente=0.5;
    multipath(n)=sqrt(pot_componente)*randn(1,1)+j*sqrt(pot_componente)*randn(1,1);
    multipath(n)=multipath(n).*sqrt(path_pot_lin(n));
    
end

RI=zeros(1,Nc);
RI(delays) = RI(delays) + multipath;

ht=RI;
Hf=fft(ht);
end

%****************SF_CODING ******************
function [ant1, ant2]=sf_coding(Nc, symbols)
ant1 = zeros(1,Nc);
ant1(1:2:end) = 1/sqrt(2) * symbols(1:2:end);
ant1(2:2:end) = 1/sqrt(2) * (-conj(symbols(2:2:end)));

ant2 = zeros(1,Nc);
ant2(1:2:end) = 1/sqrt(2) * symbols(2:2:end);
ant2(2:2:end) = 1/sqrt(2) * conj(symbols(1:2:end));
end

%****************SF_DECODING ******************
function s = sf_decoding2x1(rs1,channel11,channel21)
s = zeros(1,768);
s(1:2:end)= 1/sqrt(2).*(conj(channel11(2:2:end)).*rs1(1:2:end)+channel21(1:2:end).*conj(rs1(2:2:end)));
s(2:2:end)= 1/sqrt(2).*(conj(channel21(2:2:end)).*rs1(1:2:end)-channel11(1:2:end).*conj(rs1(2:2:end)));
end

function s=sf_decoding2x2(rx1, rx2,channel11,channel12,channel21,channel22)
s1 = zeros(1,768);
s2 = zeros(1,768);
s1(1:2:end)= 1/sqrt(2).*(conj(channel11(2:2:end)).*rx1(1:2:end)+channel21(1:2:end).*conj(rx1(2:2:end)));
s1(2:2:end)= 1/sqrt(2).*(conj(channel21(2:2:end)).*rx1(1:2:end)-channel11(1:2:end).*conj(rx1(2:2:end)));
s2(1:2:end)= 1/sqrt(2).*(conj(channel12(2:2:end)).*rx2(1:2:end)+channel22(1:2:end).*conj(rx2(2:2:end)));
s2(2:2:end)= 1/sqrt(2).*(conj(channel22(2:2:end)).*rx2(1:2:end)-channel12(1:2:end).*conj(rx2(2:2:end)));
s = s1 + s2;
end

function[ber_each_eb,ber_each_ebeq,ber_each_ebf] =trab1(EbN0)
%%
%CSF 18/19
% Daniel Cunha 68148 & Tiago Marques 68432
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference ofdm chain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocating memory & Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_OFDM_SYM =1e3;
Nc=768; % exemple, number of subcarriers
load ('pdp.mat','pdp');
fftSize=1024;
fa=15.36e6;
cp_time=5.21e-6;
m=4;
PR=1;
N_bits = Nc*2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:length(EbN0)    % EbN0_dB range
    
    
    
    for k=1:N_OFDM_SYM   % number of OFDM symbols to be simulated
        
        %################# Transmitter #################
        % data generation
        data = data_gen(N_bits);
        %data modulation
        data_mod = mod_data(data,m);
        Frame_freq=[data_mod];
        %OFDM FRAMING +CP
        dif=fftSize-length(data_mod);
        IFFT_Frame = ifft([zeros(1,dif/2) data_mod zeros(1,dif/2)]).*sqrt(fftSize);
        ta=1/(fa);
        cp_value= round(cp_time/ta);
        CP= [IFFT_Frame(fftSize-cp_value+1:fftSize)];
        Frame=[CP IFFT_Frame];
        
        
        %################# Channel Model #################
        %awgn noise
        var = PR*10^(-EbN0(p)/10)/log2(m);
        Noise = sqrt(var/2).*(randn(1,length(Frame)) + 1j*randn(1,length(Frame))) ;
        Noise_Frame=Frame+Noise;
        %multipath channel
        [ht, Hf] = channel_gen(pdp, fa, fftSize);
        Conv_sinal = conv_s_h(Frame, ht, pdp, fftSize, fa, cp_time)+Noise;
        %Equalization
        G = conj(Hf)./((abs(Hf)).^2);
        
        %Freq
        NoiseF = sqrt(var/2)*(randn(1,length(Frame_freq)) + 1j*randn(1,length(Frame_freq)));
        Freq_noise=Frame_freq.*Hf(129:896)+NoiseF;
        
        %################# Receiver #################
        %time
        %CP removal + FFT
        removCPFrame = [Conv_sinal(cp_value+1:end)];
        FFT_Frame= fft(removCPFrame)./sqrt(fftSize);
        % Equalization
        Equalized= FFT_Frame.*G;
        FrameRecived = [Equalized(129:896)] ;
        %data demodulation
        DesmodData = demod_data(FrameRecived,m,N_bits);
        
        %BER computation
        ber_ofdmeq(k) = compute_ber(DesmodData,data);
        
        %awgn
        removeCPframe2 = fft(Noise_Frame(81:end))./sqrt(fftSize);
        FrameRecived2 = removeCPframe2(129:896);
        DesmodData2 = demod_data(FrameRecived2,m,N_bits);
        ber_ofdm(k) = compute_ber(DesmodData2,data);
        
        
        %freq
        Equalized_f= Freq_noise.*G(129:896);
        DesmodDataF=demod_data(Equalized_f, m, N_bits);
        ber_ofdmf(k)=compute_ber(DesmodDataF,data);
        
    end
    
    ber_each_eb(p)= mean(ber_ofdm);
    ber_each_ebeq(p)= mean(ber_ofdmeq);
    ber_each_ebf(p)= mean(ber_ofdmf);
    
end


%  Ploting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% semilogy(EbN0,ber_each_eb,'r',EbN0,ber_each_ebeq,'k');
% grid on
% xlabel(['EbN0[dB]']), ylabel(['BER']);
% legend(['AWGN'], ['SISO OFDM']);
% figure(2)
% semilogy(EbN0,ber_each_ebeq,'b',EbN0,ber_each_ebf,'r');
% grid on
% xlabel(['EbN0[dB]']), ylabel(['BER']);
% legend(['SISO OFDM'], ['FREQ']);
%semilogy(EbN0,ber_each_eb)
%semilogy(EbN0,ber_each_ebf)
end