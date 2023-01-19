function [velocity_matriz] = crossCorrelation(data_matched_echo,segmentSize,numPointsCorr,fs,c,T_prf)
    number_of_segments = floor(size(data_matched_echo,1)/segmentSize);
    velocity_matriz = zeros(number_of_segments,size(data_matched_echo,2)-1);
    for j = 1:size(data_matched_echo,2)
        signal2 = data_matched_echo(:,j+1);
        signal1 = data_matched_echo(:,j);
        for g = 0:number_of_segments-1
            corr = xcorr(signal2((1+g*segmentSize):(segmentSize+g*segmentSize)),... 
                signal1((1+g*segmentSize):(segmentSize+g*segmentSize)));
            velocity_matriz(g+1,j) = velocityEstimator(corr,fs,c,T_prf);
        end
    end
end

