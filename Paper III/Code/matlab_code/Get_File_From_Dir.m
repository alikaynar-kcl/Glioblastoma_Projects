% Define the root directory path
rootDir = '/Users/alikaynar/Desktop/UpgrdPostDoc/GEM/GEM_R_Second_Project/';

% Initialize an empty structure for storing results
All_RES = struct();

% Get a list of all folders starting with 'model_'
folders = dir(fullfile(rootDir, 'model_*'));
folders = folders([folders.isdir]); % Filter only directories

% Loop through each folder and process the specific file
for i = 1:length(folders)
    % Construct the folder name
    folderName = folders(i).name;
    folderPath = fullfile(rootDir, folderName);
    
    % Specify the file you want to access inside this folder
    fileName = 'Res_SGD_C.mat';
    filePath = fullfile(folderPath, fileName);
    
    % Check if the file exists before trying to read or process it
    if exist(filePath, 'file')
        fprintf('Processing file: %s\n', filePath);
        % Load the .mat file
        load(filePath); % Assuming this loads Res_SGD_C variable
        
        % Use dynamic field names to store the data in the structure
        % Replace invalid characters in folderName to make it a valid field name
        validFieldName = matlab.lang.makeValidName(folderName);
        All_RES.(validFieldName) = Res_SGD_C;
    else
        fprintf('File not found: %s\n', filePath);
    end
end

% Now, All_RES contains fields named after each folder, storing their respective Res_SGD_C data
% Save the All_RES structure to a .mat file in the rootDir
saveFilePath = fullfile(rootDir, 'All_RES.mat');
save(saveFilePath, 'All_RES');

fprintf('All_RES structure saved to %s\n', saveFilePath);

