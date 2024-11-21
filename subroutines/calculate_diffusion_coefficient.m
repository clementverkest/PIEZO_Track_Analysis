function diffusion_coefficient = calculate_diffusion_coefficient(msd, time_lags)
    num_time_lags = min(10, length(msd));
    if num_time_lags < 2
        diffusion_coefficient = NaN;
        return;
    end
    
    % Perform linear regression on the first num_time_lags time lags
    p = polyfit(time_lags(1:num_time_lags), msd(1:num_time_lags), 1);
    diffusion_coefficient = p(1) / 4;
end