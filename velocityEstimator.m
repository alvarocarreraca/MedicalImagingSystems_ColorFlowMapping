function [velocity] = velocityEstimator(crossCorrelation_vector,fs,c,T_prf,velocityRange)
%// Method that calculates the velocity from one cross-correlation vector
%   by looking for the maximum value in the vector.
% ------------ Input parameters:
% - crossCorrelation_vector -> vector of cross-correlation
% - fs -> sampling frequency
% - c -> velocity of propagation in the mean
% - T_prf -> period of pulse repetition frequency
% - velocityRange -> maximum velocity in absolut value that the system can
%                    detect
% ------------ Output paramenters:
% - velocity -> estimation of the velocity from the cross-correlation

    if(max(crossCorrelation_vector) ~= 0 || min(crossCorrelation_vector) ~= 0)
        % Generate a mask for the velocity range that we are looking for:
        mask = zeros(1,length(crossCorrelation_vector));
        centre_Value_Mask = round(length(mask)/2);
        range_points_shifted = round((velocityRange*2*fs*T_prf)/c);
        if (centre_Value_Mask - 1 < range_points_shifted)
            % If the range of velocities is bigger than the mask generated
            mask = ones(1,length(crossCorrelation_vector));
        else
            mask_in_range = ones(1,range_points_shifted);
            mask(centre_Value_Mask) = 1;
            mask(centre_Value_Mask + 1:centre_Value_Mask + length(mask_in_range)) = mask_in_range;
            mask(centre_Value_Mask - length(mask_in_range):centre_Value_Mask - 1) = mask_in_range;
        end
        % Apply mask:
        crossCorrelation_vector_cut = crossCorrelation_vector.*mask';
        % Vector that centers the cross-correlation in the origin:
        x_n = -(length(crossCorrelation_vector)-1)/2:(length(crossCorrelation_vector)-1)/2;
        [~,ind] = max(crossCorrelation_vector_cut); % index of the maximum value in the cross-correlation
        points_shifted = x_n(ind);
        % Conversion from points to time:
        t_shift = points_shifted/fs;
        if (ind == 1)
            % In case after applying the maks there is a vertor of 0
            velocity = 0;
        else
            velocity = (t_shift*c)/(2*T_prf);
        end
    else
        % In case the cross-correlation is a vector of 0
        velocity = 0;
    end
end

