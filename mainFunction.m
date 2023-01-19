function [velocity_matrix] = mainFunction(data,fs,f0,cycles,c,T_prf,segmentSize,numPointsCorr,velocityRange,columnsOverlaped)
%// Method that runs all the blocks in the system and returns to the
%   original script the velocity map obtained after processing.

    data_matched = matchedFilter(data,fs,f0,cycles); %time domain
    data_matched_echo = echoCancelling(data_matched); %time domain
    velocity_matrix = crossCorrelationVariableSpan(data_matched_echo,segmentSize,numPointsCorr,fs,c,T_prf,velocityRange,columnsOverlaped);
end

