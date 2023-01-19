function [data_matched] = matchedFilter(data,fs,f0,cycles)
%// Method that filters the white noise from the signal.
% ------------ Input parameters:
% - data -> signal with noise
% - fs -> sampling frequency
% - f0 -> center frequency
% - cycles -> number of cycles in a pusle
% ------------ Output paramenters:
% - data_matched -> signal after filtration

    if(length(cycles) == 1)
        % Prepare the original pulse:
        t_pulse = linspace(0,cycles/f0,cycles*fs/f0);
        signal = sin(2*pi*f0*t_pulse);
        hanning_window = hann(length(signal));
        signal_sent = signal.*hanning_window'; 
    else
        % In the case of the simulation, cycles variable will be the pulse
        % vector:
        signal_sent = cycles; 
    end

    % Generate tranfer funtion of the filter:
    h_filter = flip(signal_sent);
    
    % Filter the data:
    data_matched = filtfilt(h_filter,1,data);
end

