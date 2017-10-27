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

RLr =(1/pi)*atan(decorrelation_strength_filter_LR.*R) + (1/2);
RRr = 1.- RLr;

%   dS_Fr = designfilt('lowpassiir', 'FilterOrder', 20, 'HalfPowerFrequency', 10000, 'SampleRate', fs);
dS_Fr = shelving(-12, 10000, fs, sqrt(0.5), 'Treble_Shelf');
S_fr = abs(freqz(dS_Fr, window_size));
S_fr = S_fr.^4;
H_A_Rr = RRr.*S_fr;
H_A_Fr = 1-H_A_Rr;

%   dS_LH = designfilt('lowpassiir', 'FilterOrder', 20, 'HalfPowerFrequency', 7000, 'SampleRate', fs);
dS_LH
S_LH = abs(freqz(dS_LH, window_size));
S_LH = S_LH.^4;
H_A_Lo = S_LH;
H_A_Hi = 1-H_A_Lo;

sigmas = [sigmaLR, S_fr, S_LH];
save kraftFilterDataSigma_original.dat sigmas -ascii;
