function [velocity_matriz] = crossCorrelationVariableSpan(data_matched_echo,segmentSize,numPointsCorr,fs,c,T_prf,velocityRange,columnsOverlaped) 
%// Method that defines a velocity map depending the resolution of the system, calculates
%   the cross-correlation for each 'pixel' of the velocity map and call the method 
%   'velocityEstimator.m' to obtain the velocity from one cross-correlation.
% ------------ Input parameters:
% - data_matched_echo -> signal after filtering out the noise
% - segmentSize -> resolution of the system in depth
% - numPointsCorr -> number of points considered to calculate the
%                    cross-correlation for one segment
% - fs -> sampling frequency
% - c -> velocity of propagation in the mean
% - T_prf -> period of pulse repetition frequency
% - velocityRange -> maximum velocity in absolut value that the system can
%                    detect
% - columnsOverlaped -> number of lines used in order to calculate one
%                       average cross-correlation, resolution in time
% ------------ Output paramenters:
% - velocity_matriz -> velocity map obtained

    % Calculate all the correlations:
    if (numPointsCorr <= segmentSize)
        diff = segmentSize - numPointsCorr;
        number_of_segments = floor(size(data_matched_echo,1)/segmentSize);
        corr_3D_matrix = zeros(2*numPointsCorr-1,size(data_matched_echo,2)-1,number_of_segments);
        for j = 1:size(data_matched_echo,2)-1
            signal2 = data_matched_echo(:,j+1);
            signal1 = data_matched_echo(:,j);
            for g = 0:number_of_segments-1
                corr = xcorr(signal2((1+g*segmentSize+ceil(diff/2)):((1+g)*segmentSize-floor(diff/2))),... 
                    signal1((1+g*segmentSize+ceil(diff/2)):((1+g)*segmentSize-floor(diff/2))));
                corr_3D_matrix(:,j,g+1) = corr;
            end
        end
    else
        diff = numPointsCorr - segmentSize;
        number_of_segments = floor((size(data_matched_echo,1)-diff)/segmentSize);
        corr_3D_matrix = zeros(2*numPointsCorr-1,size(data_matched_echo,2)-1,number_of_segments);
        for j = 1:size(data_matched_echo,2)-1
            signal2 = data_matched_echo(:,j+1);
            signal1 = data_matched_echo(:,j);
            for g = 0:number_of_segments-1
                corr = xcorr(signal2((1+g*segmentSize):(numPointsCorr+g*segmentSize)),... 
                    signal1((1+g*segmentSize):(numPointsCorr+g*segmentSize)));
                corr_3D_matrix(:,j,g+1) = corr;
            end
        end
    end
    % Calculate all the velocities:
    timeResolution = (size(data_matched_echo,2)-1)/columnsOverlaped;
    velocity_matriz = zeros(number_of_segments,timeResolution);
    for j = 1 : number_of_segments
        for i = 1 : timeResolution
            mean_corr = sum(corr_3D_matrix(:,1+(i-1)*columnsOverlaped:i*columnsOverlaped,j),2)/columnsOverlaped;
            velocity_matriz(j,i) = velocityEstimator(mean_corr,fs,c,T_prf,velocityRange);
        end
    end
end