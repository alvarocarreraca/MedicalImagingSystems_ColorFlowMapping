function [snr,snr_dB] = calculateSNR(originalSignal,noiseSignal,isFilter)
%// Method that calculates the signal to noise ratio.
% ------------ Input parameters:
% - originalSignal -> signal without noise
% - noiseSignal -> signal with noise or only noise
% - isFilter -> boolean that defines if the 'noiseSignal' is signal with noise or only noise
% ------------ Output paramenters:
% - snr -> signal-to-noise ratio
% - snr_dB -> signal-to-noise ratio in dB

    if (isFilter)
        expected_signal = mean(originalSignal.^2);
        expected_noise = mean(noiseSignal.^2);
    else
        noise = noiseSignal - originalSignal;
        expected_signal = mean(originalSignal.^2);
        expected_noise = mean(noise.^2);
    end
    snr = expected_signal./expected_noise;
    snr_dB = 10*log10(snr);
end

