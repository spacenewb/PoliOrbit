function filtered_data = filterData(data, cutoff_freq, f_s)

% data          Must be 1D array
% Cutoff_freq   Lower limit of the filter in Cycles/sec [Hz]
% f_s           Sampling frequency in Datapoints/sec [Hz]

% All freq. above cutoff will be attenuated

F_co = cutoff_freq / f_s; % Normalised Cut-off frequency
N = round( sqrt( 0.442947 + F_co^2 ) / F_co ); % Moving Window Length of Box-Car filter

if length(data) <= N
    filtered_data = data;
    disp('Cut-off Frequency too low...moving window length is larger than the data size')
else
    filtered_data = movmean(data, N);
end

end