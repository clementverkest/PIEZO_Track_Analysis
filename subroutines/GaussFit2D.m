%% 2D Gauss fit function for PIEZO cluster - Data from Trackmate - Adjuster for TrackSubtypes
function results = GaussFit2D(subtype_tracks, tiffstack, imageInfo);

results = cell(size(subtype_tracks));

% Loop through each track
 for i = 1:numel(subtype_tracks);
     % Check if the current cell is empty
    if isempty(subtype_tracks{i})
        continue;  % Skip to the next iteration if the cell is empty
    end
     track = subtype_tracks{i};
    
    % Get the first position and frame
    x_um = track.POSITION_X(1);
    y_um = track.POSITION_Y(1);
    frame = track.FRAME(1) + 1;  % Add 1 to convert from 0-299 to 1-300
    
    % Convert position from micrometers to pixels
    x_px = round(x_um / 0.110);
    y_px = round(y_um / 0.110);



%  half-width of the region (10 pixels for a 20x20 square or 5 for a 10x10)
half_width = 10;

% Define the start and end indices for the x and y dimensions
x_start = max(1, x_px - half_width + 1);
x_end = min(size(tiffstack, 2), x_px + half_width);
y_start = max(1, y_px - half_width + 1);
y_end = min(size(tiffstack, 1), y_px + half_width);

% Extract the 10X10 region from the first frame
extracted_region = tiffstack(y_start:y_end, x_start:x_end, frame);

% % Display the extracted region
% figure;
% imagesc(extracted_region);
% title('Extracted 20x20 Region');
% axis image; % Maintain aspect ratio
% colorbar;



% Assuming the extracted 20x20 region is stored in 'extracted_region'

% Get the size of the extracted region
[rows, cols] = size(extracted_region);

% Create a meshgrid for the x and y coordinates
[X_fit, Y_fit] = meshgrid(1:cols, 1:rows);

% Flatten the grid and the data for fitting
xy = [X_fit(:), Y_fit(:)];
z = double(extracted_region(:));

% Define the 2D Gaussian function
gaussian2D = @(params, xy) params(1) * exp(-((xy(:,1) - params(2)).^2 / (2 * params(3)^2) + ...
                                             (xy(:,2) - params(4)).^2 / (2 * params(5)^2))) + params(6);

% Initial guess for parameters [A, x0, sigmaX, y0, sigmaY, bg]
p0 = [max(z)-min(z), cols/2, 0.6, rows/2, 0.6, min(z)];

% Set lower and upper bounds for the parameters
lb = [0, 1, 0.1, 1, 0.1, 0]; % Lower bounds

ub = [Inf, cols, 7, rows, 7, max(z)]; % Upper bounds

% Perform the fit
options = optimoptions('lsqcurvefit', 'Display', 'off');
[params, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(gaussian2D, p0, xy, z, lb, ub, options);

% Reshape the fitted data back into a 2D array
fitted_gaussian = reshape(gaussian2D(params, xy), size(extracted_region));

  % Estimate cluster diameter (FWHM)
    fwhm_x = 2 * sqrt(2 * log(2)) * params(3) * imageInfo.pixelSize;
    fwhm_y = 2 * sqrt(2 * log(2)) * params(5) * imageInfo.pixelSize;
    cluster_diameter = mean([fwhm_x, fwhm_y]);
    

% Calculate R-squared to assess goodness of fit
    ssresid = sum(residual.^2);
    sstotal = sum((z - mean(z)).^2);
    rsq = 1 - (ssresid / sstotal);


    cluster_diameters(i) = cluster_diameter;
    results{i} = [cluster_diameter];



%% Option to display fit on cluster ! Will slow down processing !
% Display the original data and the fitted Gaussian
% figure;
% subplot(1,2,1);
% imagesc(extracted_region);
% title('Original 20x20 Region');
% colorbar;
% axis image;

% subplot(1,2,2);
% imagesc(fitted_gaussian);
% title('Fitted Gaussian');
% colorbar;
% axis image;


% subplot(1,2,2);
% imagesc(extracted_region);
% title('Fitted Gaussian');
% colormap gray
% colorbar;
% axis image;


% Calculate the center and sigma values in pixels for the overlay
center_x = params(2); % x-center of Gaussian fit
center_y = params(4); % y-center of Gaussian fit
sigma_x = params(3); % sigma in the x-direction
sigma_y = params(5); % sigma in the y-direction

% Convert to pixel coordinates (optional, if necessary)
% Assuming your `xy` is in pixel units, no conversion needed here.

% Overlay the sigma circle
% hold on;
% viscircles([center_x, center_y], 2.25*sigma_x, 'Color', 'r', 'LineWidth', 0.6,'LineStyle','--');
% viscircles([center_x, center_y], 2*sigma_x, 'Color', 'r', 'LineWidth', 0.6,'LineStyle','--');
% viscircles([center_x, center_y], 1.75*sigma_x, 'Color', 'y', 'LineWidth', 0.6,'LineStyle','--');
% viscircles([center_x, center_y], 1.5*sigma_x, 'Color', 'y', 'LineWidth', 0.6,'LineStyle','--');
% viscircles([center_x, center_y], sigma_x, 'Color', 'g', 'LineWidth', 0.6,'LineStyle','--');
% viscircles([center_x, center_y], 1.25*sigma_x, 'Color', 'g', 'LineWidth', 0.6,'LineStyle','--');
% viscircles([center_x, center_y], 0.5*sigma_x, 'Color', 'b', 'LineWidth', 0.6,'LineStyle','--');
% viscircles([center_x, center_y], 0.75*sigma_x, 'Color', 'b', 'LineWidth', 0.6,'LineStyle','--');
% hold off;



 % % Print fit results
 %    fprintf('Track %d: Diameter = %.2f μm, R^2 = %.3f\n', i, cluster_diameter, rsq);

 end




end
