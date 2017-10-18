clear
close all;

delete kraftFilterDataSigma_improved.dat;

%%  initialize parameters
%   In the plugin version, these parameters should be controlled by the host
fs = 44100;

window_size = 2048;
overlap = window_size*0.75;
fft_size = 2*window_size;

decorrelation_strength = 1;

%%  generate the decorrelation filters
R = 2*rand(window_size,1)-1;

d = designfilt('bandpassiir', 'FilterOrder', 20, 'HalfPowerFrequency1', 300, 'HalfPowerFrequency2', 10000, 'SampleRate', fs);
sigmaLR = abs(freqz(d, window_size));

decorrelation_strength_filter_LR = sigmaLR*decorrelation_strength;

RL =(1/pi)*atan(decorrelation_strength_filter_LR.*R) + (1/2);
RR = 1.- RL;

H_A_L = RL.*complex(ones(size(RL)), zeros(size(RL)));
H_A_R = RR.*complex(zeros(size(RR)), ones(size(RR)));

R = 2*rand(window_size,1)-1;
d = designfilt('bandpassiir', 'FilterOrder', 20, 'HalfPowerFrequency1', 300, 'HalfPowerFrequency2', 10000, 'SampleRate', fs);
sigma = abs(freqz(d, window_size));
decorrelation_strength_filter = sigma*decorrelation_strength;

RLr =(1/pi)*atan(decorrelation_strength_filter.*R) + (1/2);
RRr = 1.- RLr;

S_Back = ones(window_size,1);

H_A_Rr = RRr.*S_Back;
H_A_Fr = 1-H_A_Rr;

d_Lo = designfilt('lowpassiir', 'FilterOrder', 20, 'HalfPowerFrequency', 1000, 'SampleRate', fs);
d_Hi = designfilt('highpassiir', 'FilterOrder', 20, 'HalfPowerFrequency', 1000, 'SampleRate', fs);
S_Lo = abs(freqz(d_Lo, window_size));
S_Hi = abs(freqz(d_Hi, window_size));
S_Hi = S_Hi.^2;
S_Lo = S_Lo.^2;
H_A_Lo = S_Lo;
H_A_Hi = S_Hi;

sigmas = [sigmaLR, S_Back, S_Lo];
save kraftFilterDataSigma_improved.dat sigmas -ascii;
