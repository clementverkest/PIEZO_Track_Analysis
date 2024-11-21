%% CV -- Code to analyze Trackmate and TraJclassifier data of TIRF difusion of PIEZO cluster

% Trackmate data required and needed to be exported from the ImageJ plugin
% : allspots csv table, spots and tracks (optional, not used here) csv tables, and Track in
% xml format to be read/used by TraJclassifier

% TrajClassifier tables required: all four track classes and Parents tables
% All tables stored in one folder per cell/condition

%  Perform 2D Gauss fit on Cluster/particle identified in
% Trackmate: required the tiff stack where the tracks are coming from and
% to indicate the pixel size (set here at 0.110 µm)
% Note : this can be also done directly in trackmate

% Required subroutines/functions : ExtractAndAssignTracks.m, PlotAllTracks.m, calculate_msd.m,
% average_msd.m, calculate_diffusion_coefficient.m, GaussFit2D



%%
clearvars -except allResults
close all
clc


all_tracks = struct();

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data -- Loop through selected folder and Load trajectories/spots data information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
pixelSize = 0.110 % indicate the pixel size in your images. Must match with the one used when generating trackmate tables !

% Specify the folder where the files are
myfolderpath = uigetdir();

% Get a list of all files in the folder with the desired file name pattern
% CSV files
filePattern = fullfile(myfolderpath, '**/*.csv'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

% XML files
filePattern_xml = fullfile(myfolderpath, '**/*.xml'); % Change to whatever pattern you need.
theFiles_xml = dir(filePattern_xml);
tracks_data_xml = theFiles_xml(contains({theFiles_xml.name}, 'Tracks'));

% Use fileparts to get information about the tested construct/cell (if
% csv trajectories are sorted by construct/cell)
testedconstruct = strsplit(myfolderpath, '\');
testedconstruct = testedconstruct{end};


% Isolate trajectories and spots information .csv files
tracks_data = theFiles(contains({theFiles.name}, 'trajectories')); %option to have all the trajectories in one folder
confined_data = theFiles(contains({theFiles.name}, 'CONFINED'));
directed_data = theFiles(contains({theFiles.name}, 'DIRECTED'));
norm_data = theFiles(contains({theFiles.name}, 'NORM'));
subdif_data = theFiles(contains({theFiles.name}, 'SUBDIFFUSION'));
allspots_data = theFiles(contains({theFiles.name}, 'allspots'));
Parent_data = theFiles(contains({theFiles.name}, 'Parents'));
XY_data = theFiles(contains({theFiles.name}, '_spots')); %trackmate spot file containing selected tracks XY info


% Load trajectories files
% individual loading
 confined_table_path = fullfile(confined_data.folder, confined_data.name);
 confined_table=readtable(confined_table_path,'VariableNamingRule','preserve');
 confined_table.Properties.VariableNames{'PARENT-ID'} = 'PARENT_ID';
 
directed_table_path = fullfile(directed_data.folder, directed_data.name);
    % Check if the file is empty
fileInfo = dir(directed_table_path);
    if fileInfo.bytes == 0
    % If the file is empty, create an empty table or handle the situation
    directed_table = table(); % Create an empty table
    disp('Directed file is empty. Created an empty table.');
    else
 directed_table=readtable(directed_table_path,'VariableNamingRule','preserve');
 directed_table.Properties.VariableNames{'PARENT-ID'} = 'PARENT_ID';
    end

 norm_table_path = fullfile(norm_data.folder, norm_data.name);
 norm_table=readtable(norm_table_path,'VariableNamingRule','preserve');
norm_table.Properties.VariableNames{'PARENT-ID'} = 'PARENT_ID';
 
 subdif_table_path = fullfile(subdif_data.folder, subdif_data.name);
 subdif_table=readtable(subdif_table_path,'VariableNamingRule','preserve');
subdif_table.Properties.VariableNames{'PARENT-ID'} = 'PARENT_ID';
 
% Load spots (XY of tracks), allspots and Parent files
 allspots_table_path = fullfile(allspots_data.folder, allspots_data.name);
 allspots_table=readtable(allspots_table_path,'VariableNamingRule','preserve');

Parents_table_path = fullfile(Parent_data.folder, Parent_data.name);
Parents_table=readtable(Parents_table_path,'VariableNamingRule','preserve');

XY_table_path = fullfile(XY_data.folder, XY_data.name);
XY_table=readtable(XY_table_path,'VariableNamingRule','preserve');
XY_table = sortrows(XY_table, {'TRACK_ID', 'FRAME'});



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prepare for assigning and extracting tracks %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xy_array = table2array(XY_table(:, {'POSITION_X', 'POSITION_Y'}));
original_track_ids = XY_table.TRACK_ID;

% Initialize variables
    newTrackID = 1;
    start = 1;
    n = height(XY_table);
    new_track_ids = zeros(n, 1);

% Get the base name for new tracks
    [~, base_name, ~] = fileparts(XY_table_path);

 % Main loop to renumber tracks
    for i = 2:n
        if original_track_ids(i) > original_track_ids(i-1)
            % Create new track
            trackName = sprintf('%s_track%d.csv', base_name, newTrackID);
            newTrack = array2table(xy_array(start:i-1, :), 'VariableNames', {'X', 'Y'});
           

            % Update variables
            newTrackID = newTrackID + 1;
            start = i;
        end
        new_track_ids(i) = newTrackID;
    end

% Handle the last track
    trackName = sprintf('%s_track%d.csv', base_name, newTrackID);
   newTrack = array2table(xy_array(start:end, :), 'VariableNames', {'X', 'Y'});
  ;

% Update the original data with new track IDs
    XY_table.TrackID = new_track_ids ;
     XY_table.TrackID(1) = 1;



 % Extract and assign subtracks
 folderIdx = 1
 [confined_tracks, directed_tracks, diffusion_tracks, subdiffusion_tracks, all_tracks] = ExtractAndAssignTracks(folderIdx,confined_table, directed_table, norm_table, subdif_table, XY_table, all_tracks, testedconstruct);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Option to plot tracks %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PlotAllTracks(confined_tracks, directed_tracks, diffusion_tracks, subdiffusion_tracks);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract and calculate particle Intensity - average cell %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract intensity for subdiffusion tracks
    subdiffusion_intensity = zeros(numel(subdiffusion_tracks), 1);
     subdiffusion_tracklength = zeros(numel(subdiffusion_tracks), 1);
    for i = 1:numel(subdiffusion_tracks)
    current_table = subdiffusion_tracks{i};
    subdiffusion_intensity(i) = mean(current_table.MEAN_INTENSITY_CH1);
    subdiffusion_tracklength(i) = numel(current_table(:,1));
    end


% Extract intensity for diffusion tracks
    diffusion_intensity = zeros(numel(diffusion_tracks), 1);
    diffusion_tracklength = zeros(numel(diffusion_tracks), 1);
for i = 1:numel(diffusion_tracks)
    current_table = diffusion_tracks{i};
    diffusion_intensity(i) = mean(current_table.MEAN_INTENSITY_CH1);
    diffusion_tracklength(i) = numel(current_table(:,1));
end


% Extract intensity for confined tracks
    confined_intensity = zeros(numel(confined_tracks), 1);
confined_tracklength = zeros(numel(confined_tracks), 1);
for i = 1:numel(confined_tracks)
    current_table = confined_tracks{i};
    confined_intensity(i) = mean(current_table.MEAN_INTENSITY_CH1);
 confined_tracklength(i) = numel(current_table(:,1));
end


% Extract intensity for directed tracks
    directed_intensity = zeros(numel(directed_tracks), 1);
directed_tracklength = zeros(numel(directed_tracks), 1);
for i = 1:numel(directed_tracks)
    current_table = directed_tracks{i};
    directed_intensity(i) = mean(current_table.MEAN_INTENSITY_CH1);
directed_tracklength(i) = numel(current_table(:,1));
end


% Mean
meansubdiff_intensity = mean(subdiffusion_intensity);
meansubdiff_tracklength = mean(subdiffusion_tracklength);
meandiff_intensity = mean(diffusion_intensity);
meandiff_tracklength = mean(diffusion_tracklength);
meanconf_intensity = mean(confined_intensity);
meanconf_tracklength = mean(confined_tracklength);
meandir_intensity = mean(directed_intensity);
meandir_tracklength = mean(directed_tracklength);
%Median
mediansubdiff_intensity = median(subdiffusion_intensity);
mediansubdiff_tracklength = median(subdiffusion_tracklength);
mediandiff_intensity = median(diffusion_intensity);
mediandiff_tracklength = median(diffusion_tracklength);
medianconf_intensity = median(confined_intensity);
medianconf_tracklength = median(confined_tracklength);
mediandir_intensity = median(directed_intensity);
mediandir_tracklength = median(directed_tracklength);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate MSD for all the tracks and plot average and individual MSD %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize cell arrays to store MSD results
subdiffusion_msd = cell(size(subdiffusion_tracks));
diffusion_msd = cell(size(diffusion_tracks));
directed_msd = cell(size(directed_tracks));
confined_msd = cell(size(confined_tracks));

% Calculate MSD for subdiffusion tracks
for i = 1:numel(subdiffusion_tracks)
    subdiffusion_msd{i} = calculate_msd(subdiffusion_tracks{i}, 30);
end

% Calculate MSD for normal diffusion tracks
for i = 1:numel(diffusion_tracks)
    diffusion_msd{i} = calculate_msd(diffusion_tracks{i}, 30);
end

% Calculate MSD for directed tracks
for i = 1:numel(directed_tracks)
    directed_msd{i} = calculate_msd(directed_tracks{i}, 30);
end

% Calculate MSD for confined tracks
for i = 1:numel(confined_tracks)
    confined_msd{i} = calculate_msd(confined_tracks{i}, 30);
end


% Calculate average MSD for each type
[avg_subdiffusion_msd, std_subdiffusion_msd] = average_msd(subdiffusion_msd);
[avg_diffusion_msd, std_diffusion_msd] = average_msd(diffusion_msd);
[avg_directed_msd, std_directed_msd] = average_msd(directed_msd);
[avg_confined_msd, std_confined_msd] = average_msd(confined_msd);


% plot average MSD - Use "if" statement for the case of missing track type
num_MSDvalues_to_plot = 20; % adjust depending on recording time/track length
figure;

% Plot Subdiffusion
if ~isempty(avg_subdiffusion_msd)
    plot(avg_subdiffusion_msd(1:num_MSDvalues_to_plot), 'b-', 'DisplayName', 'Subdiffusion', 'LineWidth', 2);
end
hold on;

% Plot Normal Diffusion
if ~isempty(avg_diffusion_msd)
    plot(avg_diffusion_msd(1:num_MSDvalues_to_plot), 'k-', 'DisplayName', 'Normal Diffusion', 'LineWidth', 2);
end

% Plot Directed
if ~isempty(avg_directed_msd)
    plot(avg_directed_msd(1:num_MSDvalues_to_plot), 'r-', 'DisplayName', 'Directed', 'LineWidth', 2);
end

% Plot Confined
if ~isempty(avg_confined_msd)
    plot(avg_confined_msd(1:num_MSDvalues_to_plot), 'g-', 'DisplayName', 'Confined', 'LineWidth', 2);
end


% Plot the individual points with standard deviation
if ~isempty(avg_subdiffusion_msd)
errorbar(1:num_MSDvalues_to_plot, avg_subdiffusion_msd(1:num_MSDvalues_to_plot), ...
        std_subdiffusion_msd(1:num_MSDvalues_to_plot), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Subdiffusion');
end
if ~isempty(avg_diffusion_msd)
errorbar(1:num_MSDvalues_to_plot, avg_diffusion_msd(1:num_MSDvalues_to_plot), ...
        std_diffusion_msd(1:num_MSDvalues_to_plot), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Normal Diffusion');
end
if ~isempty(avg_confined_msd)
errorbar(1:num_MSDvalues_to_plot, avg_confined_msd(1:num_MSDvalues_to_plot), ...
        std_confined_msd(1:num_MSDvalues_to_plot), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Confined');
end
if ~isempty(avg_directed_msd)
errorbar(1:num_MSDvalues_to_plot, avg_directed_msd(1:num_MSDvalues_to_plot), ...
        std_directed_msd(1:num_MSDvalues_to_plot), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Directed');
end

xlabel('Time Lag');
ylabel('MSD');
title(sprintf('Average Mean Square Displacement', num_MSDvalues_to_plot));
legend('show');
hold off;


% Subdiffusion tracks
figure;
for i = 1:numel(subdiffusion_tracks)
    msd = calculate_msd(subdiffusion_tracks{i}, 40);
    if ~isempty(msd)
        plot(msd(1:num_MSDvalues_to_plot), 'b-', 'LineWidth', 1);
       hold on
    end
end
xlabel('Time Lag');
ylabel('MSD');
title('Individual Mean Square Displacement - Subdiffusion');
hold off;

% Normal Diffusion tracks
figure;
for i = 1:numel(diffusion_tracks)
    msd = calculate_msd(diffusion_tracks{i}, 40);
    if ~isempty(msd)
        plot(msd(1:num_MSDvalues_to_plot), 'k-', 'LineWidth', 1);
       hold on
    end
end
xlabel('Time Lag');
ylabel('MSD');
title('Individual Mean Square Displacement - Normal Diffusion');
hold off;

% Confined tracks
figure;
for i = 1:numel(confined_tracks)
    msd = calculate_msd(confined_tracks{i}, 40);
    if ~isempty(msd)
        plot(msd(1:num_MSDvalues_to_plot), 'g-', 'LineWidth', 0.5);
   hold on
    end
end
xlabel('Time Lag');
ylabel('MSD');
title('Individual Mean Square Displacement - Confined');
hold off;

% Directed tracks
figure;
for i = 1:numel(directed_tracks)
    msd = calculate_msd(directed_tracks{i}, 40);
    if ~isempty(msd)
        plot(msd(1:num_MSDvalues_to_plot), 'r-', 'LineWidth', 0.5);
       hold on
    end
end
xlabel('Time Lag');
ylabel('MSD');
title('Individual Mean Square Displacement - Directed');
hold off;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate diffusion coefficient %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate diffusion coefficients (i.e slope of MSD curve) for each track type
subdiffusion_coeffs = cellfun(@(track) calculate_diffusion_coefficient(track, 1:num_MSDvalues_to_plot), subdiffusion_msd, 'UniformOutput', false);
diffusion_coeffs = cellfun(@(track) calculate_diffusion_coefficient(track, 1:num_MSDvalues_to_plot), diffusion_msd, 'UniformOutput', false);
directed_coeffs = cellfun(@(track) calculate_diffusion_coefficient(track, 1:num_MSDvalues_to_plot), directed_msd, 'UniformOutput', false);
confined_coeffs = cellfun(@(track) calculate_diffusion_coefficient(track, 1:num_MSDvalues_to_plot), confined_msd, 'UniformOutput', false);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and overlay tiff image with tracks %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit identified trackmate particle/cluster with a Gauss function to
% calculate their radius/diameter
% Note : this part can also be done directly in latest trackmate version
% (v7 and after)


% Load Tiff stack
% Automatically find and load the TIFF file in the previously selected
% folder

cd(myfolderpath);
    tiffFiles = dir('*.tif');
    if isempty(tiffFiles)
        warning('No TIFF file found in folder: %s', currentFolder);
        
    end
    
    % Load the first TIFF file found
    fullFilePath = fullfile(myfolderpath, tiffFiles(1).name);
    
    % existing code for loading and processing the TIFF file 
    info = imfinfo(fullFilePath);
    num_images = numel(info);
    width = info(1).Width;
    height = info(1).Height;
    tiffstack = zeros(height, width, num_images, 'uint16');
    
    for k = 1:num_images
        tiffstack(:,:,k) = imread(fullFilePath, k);
    end




% Calculate dimensions in micrometers
width_um = width * pixelSize;
height_um = height * pixelSize;


% Create a structure to hold both pixel and micrometer information
imageInfo = struct();
imageInfo.pixelSize = pixelSize;
imageInfo.dimensions.pixels = [height, width, num_images];
imageInfo.dimensions.micrometers = [height_um, width_um, num_images];


% Overlay tracks with tiff image
% % Display the first image of the stack
figure;
imshow(tiffstack(:,:,1), [], 'InitialMagnification', 'fit'); % 1st frame
hold on;
% Set up the axis to match the image dimensions in micrometers
% axis([0 imageInfo.dimensions.micrometers(2) 0 imageInfo.dimensions.micrometers(1)]);
% Function to convert micrometers to pixels
um2px = @(x) x / imageInfo.pixelSize;

% Loop through each track in directed_tracks
for i = 1:numel(directed_tracks)
    % Get the current track
    track = directed_tracks{i};

    % Extract X and Y coordinates
    x = track.POSITION_X;
    y = track.POSITION_Y;

    % Convert coordinates from micrometers to pixels
    x_px = um2px(x);
    y_px = um2px(y);

    % Plot the track
    plot(x_px, y_px, 'LineWidth', 2);

    % Optionally, mark the start and end points
    plot(x_px(1), y_px(1), 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g');  % Start point
    plot(x_px(end), y_px(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');  % End point
end


% Loop through each track in confined_tracks
for i = 1:numel(confined_tracks)
    % Get the current track
    track = confined_tracks{i};

    % Extract X and Y coordinates
    x = track.POSITION_X;
    y = track.POSITION_Y;

    % Convert coordinates from micrometers to pixels
    x_px = um2px(x);
    y_px = um2px(y);

    % Plot the track
    plot(x_px, y_px, 'LineWidth', 2);

    % Optionally, mark the start and end points
    plot(x_px(1), y_px(1), 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g');  % Start point
    plot(x_px(end), y_px(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');  % End point
end

% Loop through each track in diffusion_tracks
for i = 1:numel(diffusion_tracks)
    % Get the current track
    track = diffusion_tracks{i};

    % Extract X and Y coordinates
    x = track.POSITION_X;
    y = track.POSITION_Y;

    % Convert coordinates from micrometers to pixels
    x_px = um2px(x);
    y_px = um2px(y);

    % Plot the track
    plot(x_px, y_px, 'LineWidth', 2);

    % Optionally, mark the start and end points
    plot(x_px(1), y_px(1), 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g');  % Start point
    plot(x_px(end), y_px(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');  % End point
end

% Loop through each track in subdiffusion_tracks
for i = 1:numel(subdiffusion_tracks)
    % Get the current track
    track = subdiffusion_tracks{i};

    % Extract X and Y coordinates
    x = track.POSITION_X;
    y = track.POSITION_Y;

    % Convert coordinates from micrometers to pixels
    x_px = um2px(x);
    y_px = um2px(y);

    % Plot the track
    plot(x_px, y_px, 'LineWidth', 2);

    % Optionally, mark the start and end points
    plot(x_px(1), y_px(1), 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g');  % Start point
    plot(x_px(end), y_px(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');  % End point
end

% Add title and labels
title('First Image of Stack with Tracks');
xlabel('X position (μm)');
ylabel('Y position (μm)');
axis equal;
% Add colorbar
colorbar;

% Adjust colormap if needed
% colormap('gray');  % Uncomment this line if you want a grayscale image

hold off;




%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional -- Create and export a movie for tracks display %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Parameters
% frameInterval = 0.1; % 100ms per frame
% pixelSize = 0.110; % µm per pixel
% scaleBarLength = 4; % desired length in µm
% scaleBarPixels = round(scaleBarLength / 0.073); % convert to pixels
% 
% % Get image dimensions
% [height, width, ~] = size(tiffstack);
% centerY = height/2;
% centerX = width/2;
% 
% % Define zoom parameters
% zoomFactor = 0.7; % How much to zoom in (smaller number = more zoom)
% finalHeight = round(height * zoomFactor);
% finalWidth = round(width * zoomFactor);
% 
% % offset - if cell is not in the center of the image
% verticalOffset = -round(height*0.1);
% 
% 
% % Calculate zoom region
% y1 = round(centerY - finalHeight/2 +verticalOffset);
% y2 = round(centerY + finalHeight/2 +verticalOffset);
% x1 = round(centerX - finalWidth/2);
% x2 = round(centerX + finalWidth/2);
% 
% % Ensure bounds are within image
% y1 = max(1, y1); y2 = min(height, y2);
% x1 = max(1, x1); x2 = min(width, x2);
% 
% % Create figure and video writer
% fig = figure('Position', [100 100 800 600]);
% videoFileName = 'tracks_movie.avi';
% writerObj = VideoWriter(videoFileName);
% writerObj.FrameRate = 20;
% open(writerObj);
% 
% % Function to convert micrometers to pixels
% um2px = @(x) x / pixelSize;
% 
% % Create axes that fill the figure
% ax = axes('Position', [0 0 1 1]);
% 
% % Loop through each frame in the stack
% for frameIdx = 1:num_images
%     % Clear the current axes
%     cla(ax);
% 
%     % Display the current frame without any margins
%     % Use only the zoomed region
%     croppedFrame = tiffstack(y1:y2, x1:x2, frameIdx);
%     imshow(croppedFrame, [], 'Parent', ax);
%     hold on;
% 
%     % Ensure the image is shown at full scale
%     axis equal;
%     axis tight;
% 
%     % Loop through each track
%     for trackIdx = 1:numel(directed_tracks)
%         track = directed_tracks{trackIdx};
%         frameNums = track.FRAME;
%         validIdx = frameNums <= frameIdx & frameNums > max(1, 0);%frameIdx - 20);
% 
%         if any(validIdx)
%             % Extract and convert coordinates
%             x = track.POSITION_X(validIdx);
%             y = track.POSITION_Y(validIdx);
%             framesValid = frameNums(validIdx);
% 
% 
% 
% 
%             x_px = um2px(x);
%             y_px = um2px(y);
% 
%             % Adjust coordinates for zoom
%             x_px = x_px - x1;
%             y_px = y_px - y1;
% 
%             % Plot only if points are within view
%             inView = x_px >= 0 & x_px <= (x2-x1) & y_px >= 0 & y_px <= (y2-y1);
% 
%             if any(inView)
%                 x_visible = x_px(inView);
%                 y_visible = y_px(inView);
%                 frames_visible = framesValid(inView);
% 
% 
%                 %fading option
%                 % if length(x_visible) > 1
%                 %     % Create color gradient for fade effect
%                 %     numPoints = length(x_visible);
%                 %     alphas = linspace(0.2, 1, numPoints);
%                 % 
%                 %     % Plot segments with increasing alpha
%                 %     for i = 1:numPoints-1
%                 %         plot(x_visible(i:i+1), y_visible(i:i+1), 'LineWidth', 4, ...
%                 %             'Color', [1 0 1 alphas(i)]);
%                 %     end
%                 % end
% 
%                 %Without fading
%                 if length(x_visible) > 1
%                     % Plot the entire visible track segment with full opacity
%                     plot(x_visible, y_visible, 'Color', [0 1 0], 'LineWidth', 1);
%                 end
% 
%                 % Mark current position if the last visible point is from current frame
%                 % if frames_visible(end) == frameIdx
%                 %     plot(x_visible(end), y_visible(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
%                 % end
%             end
%         end
%     end
% 
% % Diffusion
%     for trackIdx = 1:numel(diffusion_tracks)
%         track = diffusion_tracks{trackIdx};
%         frameNums = track.FRAME;
%         validIdx = frameNums <= frameIdx & frameNums > max(1, 0);%frameIdx - 20);
% 
%         if any(validIdx)
%             % Extract and convert coordinates
%             x = track.POSITION_X(validIdx);
%             y = track.POSITION_Y(validIdx);
%             framesValid = frameNums(validIdx);
% 
% 
% 
% 
%             x_px = um2px(x);
%             y_px = um2px(y);
% 
%             % Adjust coordinates for zoom
%             x_px = x_px - x1;
%             y_px = y_px - y1;
% 
%             % Plot only if points are within view
%             inView = x_px >= 0 & x_px <= (x2-x1) & y_px >= 0 & y_px <= (y2-y1);
% 
%             if any(inView)
%                 x_visible = x_px(inView);
%                 y_visible = y_px(inView);
%                 frames_visible = framesValid(inView);
% 
% 
%                 %fading option
%                 % if length(x_visible) > 1
%                 %     % Create color gradient for fade effect
%                 %     numPoints = length(x_visible);
%                 %     alphas = linspace(0.2, 1, numPoints);
%                 % 
%                 %     % Plot segments with increasing alpha
%                 %     for i = 1:numPoints-1
%                 %         plot(x_visible(i:i+1), y_visible(i:i+1), 'LineWidth', 4, ...
%                 %             'Color', [1 0 1 alphas(i)]);
%                 %     end
%                 % end
% 
% 
%                   %Without fading
%                 if length(x_visible) > 1
%                     % Plot the entire visible track segment with full opacity
%                     plot(x_visible, y_visible, 'Color', [0 1 0], 'LineWidth', 1);
%                 end
% 
%                 % Mark current position if the last visible point is from current frame
%                 % if frames_visible(end) == frameIdx
%                 %     plot(x_visible(end), y_visible(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
%                 % end
%             end
%         end
%     end
% 
% 
%     %Confined
% for trackIdx = 1:numel(confined_tracks)
%         track = confined_tracks{trackIdx};
%         frameNums = track.FRAME;
%         validIdx = frameNums <= frameIdx & frameNums > max(1, 0);%frameIdx - 20);
% 
%         if any(validIdx)
%             % Extract and convert coordinates
%             x = track.POSITION_X(validIdx);
%             y = track.POSITION_Y(validIdx);
%             framesValid = frameNums(validIdx);
% 
% 
% 
% 
%             x_px = um2px(x);
%             y_px = um2px(y);
% 
%             % Adjust coordinates for zoom
%             x_px = x_px - x1;
%             y_px = y_px - y1;
% 
%             % Plot only if points are within view
%             inView = x_px >= 0 & x_px <= (x2-x1) & y_px >= 0 & y_px <= (y2-y1);
% 
%             if any(inView)
%                 x_visible = x_px(inView);
%                 y_visible = y_px(inView);
%                 frames_visible = framesValid(inView);
% 
% 
%                 % %fading option
%                 % if length(x_visible) > 1
%                 %     % Create color gradient for fade effect
%                 %     numPoints = length(x_visible);
%                 %     alphas = linspace(0.2, 1, numPoints);
%                 % 
%                 %     % Plot segments with increasing alpha
%                 %     for i = 1:numPoints-1
%                 %         plot(x_visible(i:i+1), y_visible(i:i+1), 'LineWidth', 4, ...
%                 %             'Color', [0 1 0 alphas(i)]);
%                 %     end
%                 % end
% 
%                   %Without fading
%                 if length(x_visible) > 1
%                     % Plot the entire visible track segment with full opacity
%                     plot(x_visible, y_visible, 'Color', [1 0 1], 'LineWidth', 1);
%                 end
% 
%                 % Mark current position if the last visible point is from current frame
%                 % if frames_visible(end) == frameIdx
%                 %     plot(x_visible(end), y_visible(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
%                 % end
%             end
%         end
%     end
% 
%        %Subdiffusion
% for trackIdx = 1:numel(subdiffusion_tracks)
%         track = subdiffusion_tracks{trackIdx};
%         frameNums = track.FRAME;
%         validIdx = frameNums <= frameIdx & frameNums > max(1, 0);%frameIdx - 20);
% 
%         if any(validIdx)
%             % Extract and convert coordinates
%             x = track.POSITION_X(validIdx);
%             y = track.POSITION_Y(validIdx);
%             framesValid = frameNums(validIdx);
% 
% 
% 
% 
%             x_px = um2px(x);
%             y_px = um2px(y);
% 
%             % Adjust coordinates for zoom
%             x_px = x_px - x1;
%             y_px = y_px - y1;
% 
%             % Plot only if points are within view
%             inView = x_px >= 0 & x_px <= (x2-x1) & y_px >= 0 & y_px <= (y2-y1);
% 
%             if any(inView)
%                 x_visible = x_px(inView);
%                 y_visible = y_px(inView);
%                 frames_visible = framesValid(inView);
% 
% 
% 
%                 % %fading option
%                 % if length(x_visible) > 1
%                 %     % Create color gradient for fade effect
%                 %     numPoints = length(x_visible);
%                 %     alphas = linspace(0.2, 1, numPoints);
%                 % 
%                 %     % Plot segments with increasing alpha
%                 %     for i = 1:numPoints-1
%                 %         plot(x_visible(i:i+1), y_visible(i:i+1), 'LineWidth', 4, ...
%                 %             'Color', [0 1 0 alphas(i)]);
%                 %     end
%                 % end
% 
%                   %Without fading
%                 if length(x_visible) > 1
%                     % Plot the entire visible track segment with full opacity
%                     plot(x_visible, y_visible, 'Color', [1 0 1], 'LineWidth', 1);
%                 end
% 
%                 % Mark current position if the last visible point is from current frame
%                 % if frames_visible(end) == frameIdx
%                 %     plot(x_visible(end), y_visible(end), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
%                 % end
%             end
%         end
%     end
% 
% 
% 
%     % Add timestamp (white text with black edge)
%     time = frameIdx * frameInterval;
%     text(10, 20, sprintf('%.1f s', time), 'Color', 'white', ...
%         'FontSize', 12, 'FontWeight', 'bold', ...
%         'EdgeColor', 'black');
% 
%     % Add scale bar (white with black edge)
%     % Position at bottom right with some padding
%     barY = finalHeight - 30;
%     barX = finalWidth - scaleBarPixels - 20;
% 
%     % Draw scale bar
%     line([barX, barX + scaleBarPixels], [barY, barY], ...
%         'Color', 'white', 'LineWidth', 3);
%     % Add black edge
%     % line([barX, barX + scaleBarPixels], [barY, barY], ...
%     %     'Color', 'black', 'LineWidth', 4, 'Alpha', 0.5);
% 
%     % Add scale bar label
%     text(barX, barY - 20, sprintf('%d μm', scaleBarLength), ...
%         'Color', 'white', 'HorizontalAlignment', 'left', ...
%         'FontSize', 12, 'FontWeight', 'bold');
% 
%     % Remove unnecessary elements
%     set(ax, 'XTick', [], 'YTick', []);
%     axis off;
% 
%     % Ensure proper sizing
%     set(fig, 'PaperPositionMode', 'auto');
% 
%     % Capture just the image area
%     frame = getframe(ax);
%     writeVideo(writerObj, frame);
% 
%     % Optional: display progress
%     if mod(frameIdx, 10) == 0
%         fprintf('Processing frame %d/%d\n', frameIdx, num_images);
%     end
% end
% 
% % Close the video writer and figure
% close(writerObj);
% close(fig);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perfom 2D gauss fit on identified particle/cluster by trackmate %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Possibility to visualize fit: unmarked directly in GaussFit2D
% Fit done on first frame of the track
% Contained estimated diameters and R² of the fit

results1 = GaussFit2D(confined_tracks, tiffstack, imageInfo);
results2 = GaussFit2D(directed_tracks, tiffstack, imageInfo);
results3 = GaussFit2D(diffusion_tracks, tiffstack, imageInfo);
results4 = GaussFit2D(subdiffusion_tracks, tiffstack, imageInfo);



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export results as .csv tables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average MSD curves for each track subtype, with STD

    % Create the folder if it doesn't exist
    foldername = ['Results_', testedconstruct]
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

    % Find the maximum length of the arrays
max_length = max([length(avg_confined_msd), length(std_confined_msd), ...
                  length(avg_diffusion_msd), length(std_diffusion_msd), ...
                  length(avg_subdiffusion_msd), length(std_subdiffusion_msd), ...
                  length(avg_directed_msd), length(std_directed_msd)]);

    % Pad the shorter arrays with NaN values
avg_confined_msd = [avg_confined_msd; nan(max_length - length(avg_confined_msd), 1)];
std_confined_msd = [std_confined_msd; nan(max_length - length(std_confined_msd), 1)];
avg_diffusion_msd = [avg_diffusion_msd; nan(max_length - length(avg_diffusion_msd), 1)];
std_diffusion_msd = [std_diffusion_msd; nan(max_length - length(std_diffusion_msd), 1)];
avg_subdiffusion_msd = [avg_subdiffusion_msd; nan(max_length - length(avg_subdiffusion_msd), 1)];
std_subdiffusion_msd = [std_subdiffusion_msd; nan(max_length - length(std_subdiffusion_msd), 1)];
avg_directed_msd = [avg_directed_msd; nan(max_length - length(avg_directed_msd), 1)];
std_directed_msd = [std_directed_msd; nan(max_length - length(std_directed_msd), 1)];

    %Combine data for export
confined_msd_curve = [avg_confined_msd, std_confined_msd];
directed_msd_curve = [avg_directed_msd, std_directed_msd];
diffusion_msd_curve = [avg_diffusion_msd, std_diffusion_msd];
subdiffusion_msd_curve = [avg_subdiffusion_msd, std_subdiffusion_msd];
All_msd_curve = [diffusion_msd_curve,directed_msd_curve,subdiffusion_msd_curve,confined_msd_curve];
All_msd_curve_table = array2table(All_msd_curve);
headers = {'Diffusion_MSD','Diffusion_MSD_STD','Directed_MSD','Directed_MSD_STD', 'Subdiffusion_MSD', 'Subdiffusion_MSD_STD', 'Confined_MSD','Confined_MSD_STD'};
   All_msd_curve_table.Properties.VariableNames=headers;
    writetable(All_msd_curve_table, fullfile(foldername, 'AverageMSDandStdDev.csv'))

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track and cluster informations table

        %Number of tracks
         num_diffusiontrack = size(norm_table,1);
        num_subdiffusiontrack = size(subdif_table,1);
        num_directedtrack = size(directed_table,1);
        num_confinedtrack = size(confined_table,1);

        %number of clusters on 1st frame
        num_cluster = sum(allspots_table.FRAME == 0);

        %Cluster size
        average_cluster_diameter = mean(allspots_table.RADIUS(allspots_table.FRAME == 0));
        
        %Diffusion coefficient average
        diffcoef_norm = mean(cell2mat(diffusion_coeffs), 'omitnan')*10; % x10 to get a um²/s value
        diffcoef_subdif = mean(cell2mat(subdiffusion_coeffs), 'omitnan')*10;
        diffcoef_confined = mean(cell2mat(confined_coeffs), 'omitnan')*10;
        diffcoef_directed = mean(cell2mat(directed_coeffs), 'omitnan')*10;

        %Diameter Gauss fit - average per cell
        cluster_diameter = zeros(length(results1), 1);
        for i= 1:length(results1);
            cluster_diameter(i)=results1{i}(1);
        end
        diameter_conf = mean(cluster_diameter);

          cluster_diameter = zeros(length(results2), 1);
        for i= 1:length(results2);
            cluster_diameter(i)=results2{i}(1);
        end
        diameter_dir = mean(cluster_diameter);


         cluster_diameter = zeros(length(results3), 1);
        for i= 1:length(results3);
            cluster_diameter(i)=results3{i}(1);
        end
        diameter_dif = mean(cluster_diameter);




          cluster_diameter = zeros(length(results4), 1);
        for i= 1:length(results4);
            cluster_diameter(i)=results4{i}(1);
        end
        diameter_subdif = mean(cluster_diameter);
        
% Option - Get all the cluster diameter       
% Create structure to store all the determined cluster diameter
       % Check if the structure array already exists in the workspace
if ~exist('allResults', 'var')
    allResults = struct('results1', {}, 'results2', {}, 'results3', {}, 'results4', {});
end

% Assuming result1, result2, result3, and result4 are generated in this run
% Assign the results into the structure
newEntry.results1 = results1;
newEntry.results2 = results2;
newEntry.results3 = results3;
newEntry.results4 = results4;

% Append the new results to the structure array
allResults(end + 1) = newEntry;

% Initialize an empty array to store the extracted values
extractedValues = [];

% Loop through each structure in the array
for i = 1:numel(allResults)
    % Loop through each field in the structure
    for fieldName = ["results1", "results2", "results3", "results4"]
        % Get the current cell array
        cellArray = allResults(i).(fieldName);

        % Loop through each cell in the cell array
        for j = 1:numel(cellArray)
            % Get the first number from the 1x2 array (first row, first column)
            extractedNumber = cellArray{j}(1);

            % Append the number to the extractedValues array
            extractedValues(end + 1) = extractedNumber;
        end
    end
end



        track_cluter_info = [num_diffusiontrack, num_subdiffusiontrack, num_directedtrack, num_confinedtrack, num_cluster, average_cluster_diameter, diffcoef_norm, diffcoef_subdif, diffcoef_confined, diffcoef_directed, meansubdiff_intensity, meansubdiff_tracklength, meandiff_intensity, meandiff_tracklength, meanconf_intensity, meanconf_tracklength, meandir_intensity, meandir_tracklength, diameter_conf, diameter_dir, diameter_dif, diameter_subdif];
        track_cluter_info_table = array2table(track_cluter_info);
headers = {'Nb norm_diff track', 'Nb subdiff track','Nb directed track','Nb confined track', 'Nb cluster 1st frame', 'Avg cluster diameter', 'diffusion coeff norm_diff','diffusion coeff subdif', 'diffusion coeff confined', 'diffusion coeff directed', 'mean subdif intensity','mean subdif tracklength', 'mean dif intensity','mean dif tracklength', 'mean conf intensity','mean conf tracklength', 'mean dir intensity','mean dir tracklength', 'Avg diameter conf','Avg diameter dir','Avg diameter dif','Avg diameter subdif'};
   track_cluter_info_table.Properties.VariableNames=headers;
    writetable(track_cluter_info_table, fullfile(foldername, 'TrackAndClusterInformation.csv'))


   
