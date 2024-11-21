%% Extract MSD info and do average for all cells

% Select the main folder
mainFolder = uigetdir('Select the main folder');

% Get a list of all subfolders
subFolders = dir(mainFolder);
subFolders = subFolders([subFolders.isdir]);  % Keep only directories
subFolders = subFolders(~ismember({subFolders.name},{'.','..'}));  % Remove . and ..

% Initialize matrix to store data - adjust number of rows depending on
% recording and track length
allData = zeros(20, 4, length(subFolders));  % 20 rows, 4 columns, number of files

% Loop through each subfolder
for i = 1:length(subFolders)
    % Construct the full path to the current subfolder
    currentFolder = fullfile(mainFolder, subFolders(i).name);
    
    % Look for the CSV file in the current folder
    csvFile = fullfile(currentFolder, 'AverageMSDandStdDev.csv');
    
    % Check if the file exists
    if exist(csvFile, 'file')
        % Read the CSV file
        currentData = readtable(csvFile);
        
        % Extract rows 1 to 20 of columns 1, 3, 5, and 7 - adjust number of rows depending on
% recording and track length
        extractedData = currentData{1:min(20,height(currentData)), [1 3 5 7]};
        
        % Store the extracted data
        allData(1:size(extractedData,1), :, i) = extractedData;
    else
        warning('CSV file not found in folder: %s', subFolders(i).name);
    end
end
%%
% Calculate average for each row across all files
avgData = mean(allData, 3,"omitmissing");
avgData_std = std(allData, 0, 3, "omitmissing");
avgAllData = [avgData avgData_std];
%%
% Create the final table
finalTable = array2table(avgAllData, 'VariableNames', {...
    'Avg_MSD_NormDiff', 'Avg_MSD_DIrected', 'Avg_MSD_Subdiff', 'Avg_MSD_Confined', 'STD_MSD_NormDiff', 'STD_MSD_DIrected', 'STD_MSD_Subdiff', 'STD_MSD_Confined'});

% Display the resulting table
disp(finalTable);

% Ask the user where to save the Excel file
[fileName, filePath] = uiputfile('*.xlsx', 'Save Excel File');

if fileName ~= 0  % Check if user didn't cancel the file dialog
    fullFilePath = fullfile(filePath, fileName);
    
    % Write the table to an Excel file
    writetable(finalTable, fullFilePath, 'Sheet', 'AverageMSDandStdDev');
    
    fprintf('Data successfully exported to: %s\n', fullFilePath);
else
    fprintf('Excel file export cancelled by user.\n');
end