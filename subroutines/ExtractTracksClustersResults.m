

%% Loop and extract tracks/cluster data

clear all
close all
%%
% Select the main folder
mainFolder = uigetdir('Select the main folder');

% Get a list of all subfolders
subFolders = dir(mainFolder);
subFolders = subFolders([subFolders.isdir]);  % Keep only directories
subFolders = subFolders(~ismember({subFolders.name},{'.','..'}));  % Remove . and ..

% Initialize an empty table to store all data
allData = table();

% Loop through each subfolder
for i = 1:length(subFolders)
    % Construct the full path to the current subfolder
    currentFolder = fullfile(mainFolder, subFolders(i).name);
    
    % Look for the CSV file in the current folder
    csvFile = fullfile(currentFolder, 'TrackAndClusterInformation.csv');
    
    % Check if the file exists
    if exist(csvFile, 'file')
        % Read the CSV file
        currentData = readtable(csvFile);
        
        % Add a column to identify the source folder (optional)
        currentData.SourceFolder = repmat({subFolders(i).name}, height(currentData), 1);
        
        % Append to the main table
        allData = [allData; currentData];
    else
        warning('CSV file not found in folder: %s', subFolders(i).name);
    end
end

% Ask the user where to save the Excel file
[fileName, filePath] = uiputfile('*.xlsx', 'Save Excel File');

if fileName ~= 0  % Check if user didn't cancel the file dialog
    fullFilePath = fullfile(filePath, fileName);
    
    % Write the table to an Excel file
    writetable(allData, fullFilePath, 'Sheet', 'CombinedData');
    
    fprintf('Data successfully exported to: %s\n', fullFilePath);
else
    fprintf('Excel file export cancelled by user.\n');
end