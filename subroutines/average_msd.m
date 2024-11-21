% Function to calculate average MSD
function [avg_msd, std_msd] = average_msd(msd_cell_array)
    % Find the maximum length of MSD
    max_length = max(cellfun(@length, msd_cell_array));
    
    % Initialize sum and count arrays
    sum_msd = zeros(max_length, 1);
    count = zeros(max_length, 1);
    sum_squared_msd = zeros(max_length, 1);


    % Sum up all MSDs, ignoring empty arrays
    for i = 1:numel(msd_cell_array)
        msd = msd_cell_array{i};
        if ~isempty(msd)
            sum_msd(1:length(msd)) = sum_msd(1:length(msd)) + msd;
            sum_squared_msd(1:length(msd)) = sum_squared_msd(1:length(msd)) + msd.^2;
            count(1:length(msd)) = count(1:length(msd)) + 1;
        end
    end
    
    % Calculate average and standard deviation, avoiding division by zero
    avg_msd = sum_msd ./ max(count, 1);  % Use max to avoid division by zero
    std_msd = sqrt(sum_squared_msd ./ max(count, 1) - avg_msd.^2);
    
    % Set to NaN where count is zero
    avg_msd(count == 0) = NaN;
    std_msd(count == 0) = NaN;
end
