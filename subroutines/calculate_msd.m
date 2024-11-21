function msd = calculate_msd(track, min_frames)
    
 % If min_frames is not provided, set a default value
    if nargin < 2
        min_frames = 40;  % Default minimum number of frames
    end
    
    % Check if the track meets the minimum frame requirement
    if height(track) < min_frames
        msd =[];
        return;
    end



%  Get columns for X and Y from track table
    x = track.POSITION_X;
    y = track.POSITION_Y;
    
    n = height(track);
    msd = zeros(n-1, 1);
    
    for tau = 1:(n-1)
        dx = x(tau+1:end) - x(1:end-tau);
        dy = y(tau+1:end) - y(1:end-tau);
        msd(tau) = mean(dx.^2 + dy.^2);
    end
end