% Script to run all the simulations:

close all; clear all; clc;

%% PREPARE SIMULATION:

% PULSE:
load pulse.mat
t = 1/fs:1/fs:(1/fs*length(pulse)); % time axis of the pulse
pulse_freq = fft(hilbert(pulse));
pulse_freq_abs = abs(pulse_freq);
freq_axis = linspace(0,fs,length(pulse_freq_abs));
% Centre frequency of the transducer:
[~, ind] = max(pulse_freq_abs);
f0 = freq_axis(ind);

% SCATTERING IN THE VESSEL:
d = 10*10^-3; % vessel diameter
ang = pi/4; % angle of the transducer
c = 1500; % m/s
lambda = c/fs; % wave length
depth = sqrt(2*d^2); % depth in the vessel
nsamples = round(depth/lambda); % 943 wavelengths -> penetration in the vessel

% SIMULATION OF SCATTERING:
scatter = randn(nsamples,1); % scatter (column vector 943x1)
rf_signal = conv(pulse,scatter); % RF signal
time_rf_signal = 0:1/fs:1/fs*(length(rf_signal)-1); % time axis of RF signal
rf_signal_spectrum = fft(hilbert(rf_signal)); % Spectrum of RF signal
freq_rf_signal = linspace(0,fs,length(rf_signal_spectrum)); % frequency axis of the spectrum

% SIMULATION OF PLUG FLOW WITH 100 SIGNALS RECEIVED (contant velocity):
num_signals = 100; % number of signals received
vz = 0.15; % velocity of the flow in the vessel
fprf = 5*10^3; % pulse repetition frequnecy
Tprf = 1/fprf; % period of pulse repition frequency
timeShift = 2*vz*Tprf/c; % time shift generated in the pulse by the velocity of the blood

% PLOTS:
%Plot pulse:
figure();
subplot(1,2,1);
plot(t,pulse)
xlabel("Time [s]")
ylabel("Voltage [V]")
title("Pulse")
% Plot spectrum of the pulse:
subplot(1,2,2);
plot(freq_axis,pulse_freq_abs)
xlabel("Frequency [Hz]")
title("Pulse spectrum")
ylabel("Power [W]")
xlim([0 3e7])

signal = delayseq(rf_signal,(1:num_signals)*timeShift,fs); % signals received
figure();
imagesc([0 Tprf*num_signals],[0 depth],signal);
colormap (hot)
colorbar
title("Emitted signal")
xlabel("Time [s]")
ylabel("Depth [m]")

figure();
plot(time_rf_signal,rf_signal)
xlabel("Time [s]")
title("Model of RF signal behaviour in the blood")
ylabel("Voltage [V]")

%% ECHO CANCELLING TEST:
% Prepare signals with white noise and stationary signals:
stationary_echo = 20*randn(length(rf_signal),1); % between 10 and 100 times 
noise_signal = signal + stationary_echo;
for i = 1:num_signals
    white_noise = randn(length(rf_signal),1);
    noise_signal(:,i) = noise_signal(:,i) + white_noise;
end
% Cancel echo:
signalRecovered = echoCancelling (noise_signal);

% Calculate results:
signal_index = 20;
[snr_before_echo_filter,snr_dB_before_echo_filter] = calculateSNR(signal,noise_signal,0);
[snr_after_echo_filter,snr_dB_after_echo_filter] = calculateSNR(signal,signalRecovered,0);
snr_echo_improve = snr_dB_after_echo_filter - snr_dB_before_echo_filter;

g = [repmat({'SNR dB before filter'},100,1); repmat({'SNR dB after filter'},100,1); repmat({'SNR imporve'},100,1)];
figure;
boxplot([snr_dB_before_echo_filter'; snr_dB_after_echo_filter'; snr_echo_improve'],g,'color','brg')
ylabel("SNR [dB]")
title("SNR before and afer echo filter")

% Plot results:
figure;
subplot(1,2,1);
imagesc([0 Tprf*num_signals],[0 depth],noise_signal);
colormap (hot)
colorbar
title("Signal received with noise")
xlabel("Time [s]")
ylabel("Depth [m]")

subplot(1,2,2);
imagesc([0 Tprf*num_signals],[0 depth],signalRecovered);
colormap (hot)
colorbar
title("Signal after echo cancelling")
xlabel("Time [s]")
ylabel("Depth [m]")

% --------- Represent one line (pulse)----------
% Time domain:
num = 20;
figure;
subplot(1,2,1);
plot(time_rf_signal, signal(:,num), 'b')
hold on
plot(time_rf_signal, noise_signal(:,num), 'g')
plot(time_rf_signal, signalRecovered(:,num), 'r')
title("Time domain")
xlabel("Time [s]")
ylabel("Voltage [V]")
xlim([0.2e-5 0.45e-5])
legend('Original signal','Signal received','Signal filtered')
% Frequnecy domain:
freq = linspace(0,fs,length(signal(:,num)));
subplot(1,2,2);
plot(freq, abs(fft(hilbert(signal(:,num)))), 'b')
hold on
plot(freq, abs(fft(hilbert(noise_signal(:,num)))), 'g')
plot(freq, abs(fft(hilbert(signalRecovered(:,num)))), 'r')
title("Frequency domain")
xlabel("Frequency [Hz]")
ylabel("Power [W]")
xlim([0 1e7])
legend('Original signal','Signal received','Signal filtered')


%% MATCHED FILTER TEST:
% Prepare signals with white noise and stationary signals:
whiteNoise_signal = signal;
whiteNoise_matrix = zeros(size(signal,1),size(signal,2));
for i = 1:num_signals
    white_noise = randn(length(rf_signal),1);
    whiteNoise_signal(:,i) = whiteNoise_signal(:,i) + white_noise;
    whiteNoise_matrix (:,i) = white_noise;
end

% Generate tranfer funtion of the filter:
h_filter = flip(pulse);
% Filter the data:
data_matched = filtfilt(h_filter,1,whiteNoise_signal);
noise_matched = filtfilt(h_filter,1,whiteNoise_matrix); % filt the noise
signal_matched = filtfilt(h_filter,1,signal); % filt the original signal

% Calculate results:
[snr_before_matched_filter,snr_dB_before_matched_filter] = calculateSNR(signal,whiteNoise_matrix,1);
[snr_after_matched_filter,snr_dB_after_matched_filter] = calculateSNR(signal_matched,noise_matched,1);
snr_matched_improve = snr_dB_after_matched_filter - snr_dB_before_matched_filter;

gMatch = [repmat({'SNR dB before filter'},100,1); repmat({'SNR dB after filter'},100,1); repmat({'SNR imporve'},100,1)];
figure;
boxplot([snr_dB_before_matched_filter'; snr_dB_after_matched_filter'; snr_matched_improve'],gMatch,'color','brg')
ylabel("SNR [dB]")
title("SNR before and afer matched filter")

% Plot results:
figure;
subplot(1,2,1);
imagesc([0 Tprf*num_signals],[0 depth],whiteNoise_signal);
colormap (hot)
colorbar
title("Signal received with noise")
xlabel("Time [s]")
ylabel("Depth [m]")

subplot(1,2,2);
imagesc([0 Tprf*num_signals],[0 depth],data_matched);
colormap (hot)
colorbar
title("Signal after matched filter")
xlabel("Time [s]")
ylabel("Depth [m]")

% --------- Represent one line (pulse)----------
% Time domain:
num = 20;
figure;
subplot(1,2,1);
plot(time_rf_signal, signal(:,num)/max(signal(:,num)), 'b')
hold on
plot(time_rf_signal, whiteNoise_signal(:,num)/max(whiteNoise_signal(:,num)), 'g')
plot(time_rf_signal, data_matched(:,num)/max(data_matched(:,num)), 'r')
title("Time domain")
xlabel("Time [s]")
ylabel("Voltage [V]")
xlim([0.2e-5 0.45e-5])
legend('Original signal','Signal received','Signal filtered')
% Frequnecy domain:
freq = linspace(-fs/2,fs/2,length(signal(:,num)));
subplot(1,2,2);
plot(freq, abs(fftshift(fft(signal(:,num))))/max(abs(fftshift(fft(signal(:,num))))), 'b')
hold on
plot(freq, abs(fftshift(fft(whiteNoise_signal(:,num))))/max(abs(fftshift(fft(whiteNoise_signal(:,num))))), 'g')
plot(freq, abs(fftshift(fft(data_matched(:,num))))/max(abs(fftshift(fft(data_matched(:,num))))), 'r')
title("Frequency domain")
xlabel("Frequency [Hz]")
ylabel("Power [W]")
xlim([0 2e7])
legend('Original signal','Signal received','Signal filtered')


%% VELOCITY ESTIMATOR TEST:
segmentSize = 10;
numPointsCorr = 100;
velocityRange = 1;
columnsOverlaped = 1;

velocity_matrix = crossCorrelationVariableSpan(signal,segmentSize,numPointsCorr,fs,c,Tprf,velocityRange,columnsOverlaped);

figure;
imagesc([0 Tprf*num_signals],[0 depth],velocity_matrix);
colormap (jet)
colorbar
title("")
xlabel("Time [s]")
ylabel("Depth [m]")


figure;
plot(linspace(0,depth,length(velocity_matrix(:,10))),velocity_matrix(:,10))
title("Velocity profile (Segment: 10; Span: 100; Line Overlap: 1)")
xlabel("Depth [m]")
ylabel("Velocity [m/s]")

% probability of detection:
figure;
hist_velocity = histogram(velocity_matrix,15);
ylabel("Number of points")
xlabel("Velocity [m/s]")

detectionProbability = max(hist_velocity.Values)/sum(hist_velocity.Values)

%% VELOCITY DETECTION SYSTEM TEST:
segmentSize = 5;
numPointsCorr = 126;
velocityRange = 1;
columnsOverlaped = 3;

velocity_detected = mainFunction(noise_signal,fs,f0,pulse,c,Tprf,segmentSize,numPointsCorr,velocityRange,columnsOverlaped);

figure;
imagesc([0 Tprf*num_signals],[0 depth],velocity_detected);
colormap (jet)
colorbar
title("Velocity map (Probability of detection: 64.03%)")
xlabel("Time [s]")
ylabel("Depth [m]")


figure;
plot(linspace(0,depth,length(velocity_detected(:,14))),velocity_detected(:,14))
title("Velocity profile (Segment: 5; Span: 126; Line Overlap: 3)")
xlabel("Depth [m]")
ylabel("Velocity [m/s]")

% probability of detection:
figure;
hist_velocity_system = histogram(velocity_matrix,15);
ylabel("Number of points")
xlabel("Velocity [m/s]")

detectionProbabilitySystem = max(hist_velocity_system.Values)/sum(hist_velocity_system.Values)
