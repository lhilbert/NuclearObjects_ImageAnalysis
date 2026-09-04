%% Collect experiment image metadata
%
% This script collects the metadata from the experiment image data files and
% writes it to a CSV file. Metadata is only read from the file headers, the
% image data itself is not read.
%
% The metadata gives an overview of the experiment data and can help identify
% corrupt image files.
%
% The script requires the BioFormats toolbox and the OMEImageReaderLazy class to
% be accessible on the MATLAB path. They can be found here:
% - BioFormats toolbox: www.openmicroscopy.org/bio-formats/
% - OMEImageReaderLazy: <to be done>
%
% If the MATLAB path can't be adjusted in the MATLAB client software (e.g. on a
% compute cluster), add the respective paths during runtime.

% addpath(fullfile(".", "bfmatlab"));
% addpath(fullfile(".", "OMEImageReaderLazy"));

%% Process parameter section

% image data root directory
dataRootDirectory = fullfile("path", "to", "data");

% data file path pattern (for absolute file paths)
% can contain placeholders "{ParameterName}" for filepath-encoded experiment
% parameters; those must then be specified as shown below
fileSelector = fullfile(dataRootDirectory, "Data_{Day}", "{Condition}{CellLine}*.nd2");

% experiment data using the DatasetFileManager class
fileManager = DatasetFileManager();
fileManager = fileManager.addParameter("Day");
fileManager = fileManager.addLabel("Day", "2026-Jan-01", "010126", "Day1");
fileManager = fileManager.addLabel("Day", "2026-Feb-28", "280226", "Day2");
fileManager = fileManager.addParameter("Condition");
fileManager = fileManager.addLabel("Condition", "Control"    , "Ctrl"  , "Cond0");
fileManager = fileManager.addLabel("Condition", "Condition 1", "Shake" , "Cond1");
fileManager = fileManager.addLabel("Condition", "Condition 2", "Rattle", "Cond2");
fileManager = fileManager.addLabel("Condition", "Condition 3", "Roll"  , "Cond3");
fileManager = fileManager.addLabel("Condition", "Condition 4", "Heat"  , "Cond4");
fileManager = fileManager.addLabel("Condition", "Condition 5", "Freeze", "Cond5");
fileManager = fileManager.addParameter("CellLine");
fileManager = fileManager.addLabel("CellLine", "Cell line 1", "A", "Cell1");
fileManager = fileManager.addLabel("CellLine", "Cell line 2", "B", "Cell2");
fileManager.OriginalFilepathPattern = fileSelector;

% metadata output file
metadataFile = fullfile(dataRootDirectory, "metadata.csv");

%% Helper functions

% creates the result metadata table structure (including experiment parameters)
function tab = init_metadata_table(tab_par_combos)
    par_names = string(tab_par_combos.Properties.VariableNames');
    par_types = tab_par_combos.Properties.VariableTypes';
    tabVars = [ ...
        "Filepath"   , "string" ; ...
        par_names    , par_types; ...
        "SeriesTotal", "int32"  ; ...
        "SeriesInd"  , "int32"  ; ...
        "SizeC"      , "int32"  ; ...
        "SizeX"      , "int32"  ; ...
        "SizeY"      , "int32"  ; ...
        "SizeZ"      , "int32"  ; ...
        "SizeT"      , "int32"  ; ...
        "VoxelSizeX" , "double" ; ...
        "VoxelSizeY" , "double" ; ...
        "VoxelSizeZ" , "double" ; ...
        ]';
    tabSize = [0, size(tabVars, 2)];
    tab = table(Size=tabSize, VariableNames=tabVars(1,:), VariableTypes=tabVars(2,:));
end

function tab = append_metadata_entry(tab, file_list_row, image_reader, series_ind)
    par_vals = arrayfun(@(a) a, file_list_row{1,:}, UniformOutput=false);
    num_series = image_reader.getNumSeries();
    num_channels = image_reader.getNumChannels(series_ind);
    image_size = image_reader.getStackSizeXYZ(series_ind);
    num_time = image_reader.getSizeT(series_ind);
    voxel_size = [image_reader.getPixelSizeXY(series_ind), image_reader.getZStepSize(series_ind)];
    tab(end + 1, :) = [ ...
        par_vals, ...       % Filepath and experiment parameters
        { ...
        num_series, ...     % SeriesTotal
        series_ind, ...     % SeriesInd
        num_channels, ...   % SizeC
        image_size(1), ...  % SizeX
        image_size(2), ...  % SizeY
        image_size(3), ...  % SizeZ
        num_time, ...       % SizeT
        voxel_size(1), ...  % VoxelSizeX
        voxel_size(2), ...  % VoxelSizeY
        voxel_size(3), ...  % VoxelSizeZ
        }];
end

%% Main script section

% metadata result table
metadataTable = init_metadata_table(fileManager.compileParameterCombinations());

% file list
fileList = fileManager.compileOriginalFileTable();

numFiles = height(fileList);

for ff = 1:numFiles

    filepath = fileList.Filepath(ff);
    assert(isfile(filepath))
    reader = OMEImageReaderLazy(filepath);
    numSeries = reader.getNumSeries(); % Fast, uses metadata only

    fprintf("File %d of %d (%s, contains %d series)\n", ff, numFiles, filepath, numSeries)

    for ss = 1:numSeries

        metadataTable = append_metadata_entry(metadataTable, fileList(ff,:), reader, ss);

    end

    reader.close(); % Always close reader after use!

end

folder = fileparts(metadataFile);
if ~isfolder(folder)
    mkdir(folder);
end

writetable(metadataTable, metadataFile);
