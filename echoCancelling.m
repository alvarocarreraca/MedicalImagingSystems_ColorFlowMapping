function [data_filtered] = echoCancelling(data)
%// Method that filters the stationary echo from the signal.
% ------------ Input parameters:
% - data -> signal with statoinary echo
% ------------ Output paramenters:
% - data_filtered -> signal after removing the echo

    stationary_signal = mean(data,2);
    data_filtered = data - stationary_signal;
end

