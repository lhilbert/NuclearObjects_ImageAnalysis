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
% can contain placeholders "{P1}", "{P2}", "{P3}", ... for filepath-encoded
% experiment parameters; those must then be specified as shown below
fileSelector = fullfile(dataRootDirectory, "Data_{P1}", "{P2}{P3}*.nd2");

% experiment parameter specifications as structs with the following fields:
%   - name: (string scalar, without spaces)
%       descriptive name of the parameter (appears in the metadata table head)
%   - items: string array with three columns
%       column 1: descriptive name for each parameter expression
%       column 2: filepath-encoded representation in the source data files
%       column 3: filepath-encoded representation in the extracted stack files
experimentParameters(1).name = "Day";
experimentParameters(1).items = [ ...
	"2026-Jan-01", "010126", "Day1"; ...
	"2026-Feb-28", "280226", "Day2"; ...
	];
experimentParameters(2).name = "Condition";
experimentParameters(2).items = [ ...
	"Control"    , "Ctrl"  , "Cond0"; ...
	"Condition 1", "Shake" , "Cond1"; ...
	"Condition 2", "Rattle", "Cond2"; ...
	"Condition 3", "Roll"  , "Cond3"; ...
	"Condition 4", "Heat"  , "Cond4"; ...
	"Condition 5", "Freeze", "Cond5"; ...
	];
experimentParameters(3).name = "CellLine";
experimentParameters(3).items = [ ...
	"Cell line 1", "A", "Cell1"; ...
	"Cell line 2", "B", "Cell2"; ...
	];

% metadata output file
metadataFile = fullfile(dataRootDirectory, "metadata.csv");

%% Helper functions

% creates the result metadata table structure (including experiment parameters)
function tab = init_metadata_table(par_struct)
    parVars = [ ...
        [par_struct.name]; ...
        compose("%sFileDesc", [par_struct.name]); ...
        ];
    tabVars = [ ...
        "Filepath",    "string"; ...
        parVars(:),    repmat("string",numel(parVars),1); ...
        "SeriesTotal", "int32"; ...
        "SeriesInd",   "int32"; ...
        "SizeC",       "int32"; ...
        "SizeX",       "int32"; ...
        "SizeY",       "int32"; ...
        "SizeZ",       "int32"; ...
        "SizeT",       "int32"; ...
        "VoxelSizeX",  "double"; ...
        "VoxelSizeY",  "double"; ...
        "VoxelSizeZ",  "double"; ...
        ]';
    tabSize = [0, size(tabVars, 2)];
    tab = table(Size=tabSize, VariableNames=tabVars(1,:), VariableTypes=tabVars(2,:));
end

function vals = retrieve_parameter_values(par_struct, row_selectors, col)
    assert(all(size(row_selectors) == size(par_struct)))
    vals = strings(size(row_selectors));
    for pp = 1:length(par_struct)
        vals(pp) = par_struct(pp).items(row_selectors(pp), col);
    end
end

function tab = compile_file_table(par_struct, filepath_pattern)
    n_pars = length(par_struct);
    par_ind_func = @(m) 1:size(m, 1);
    par_inds = cellfun(par_ind_func, {par_struct.items}, UniformOutput=false);
    par_ind_combos = combinations(par_inds{:});
    par_ind_combos.Properties.VariableNames = compose("Par%dInd", 1:n_pars);
    filepath_placeholders = compose("{P%d}", 1:n_pars);
    tab = table( ...
        Size=[0,n_pars+1], ...
        VariableTypes=["string",repmat("double",1,n_pars)], ...
        VariableNames=["Filepath",par_ind_combos.Properties.VariableNames]);
    for pp = 1:height(par_ind_combos)
        filepath_values = retrieve_parameter_values(par_struct, par_ind_combos{pp,:}, 2);
        file_pattern = replace(filepath_pattern, filepath_placeholders, filepath_values);
        files = dir(file_pattern);
        files = fullfile(string({files.folder}), string({files.name}))';
        new_rows = [table(files, VariableNames="Filepath"), repmat(par_ind_combos(pp,:), length(files), 1)];
        tab = [tab; new_rows]; %#ok<AGROW>
    end
end

function tab = append_metadata_entry(tab, file_list_row, par_struct, image_reader, series_ind)
    filepath = file_list_row.Filepath;
    file_parameters = [ ...
        retrieve_parameter_values(par_struct, file_list_row{1, 2:end}, 1); ...
        retrieve_parameter_values(par_struct, file_list_row{1, 2:end}, 3);];
    par_vals = mat2cell(file_parameters(:), ones(numel(file_parameters), 1));
    num_series = image_reader.getNumSeries();
    num_channels = image_reader.getNumChannels(series_ind);
    image_size = image_reader.getStackSizeXYZ(series_ind);
    num_time = image_reader.getSizeT(series_ind);
    voxel_size = [image_reader.getPixelSizeXY(series_ind), image_reader.getZStepSize(series_ind)];
    tab(end + 1, :) = { ...
        filepath, ...       % FileName
        par_vals{:}, ...    % Experiment parameters
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
        };
end

%% Main script section

% metadata result table
metadataTable = init_metadata_table(experimentParameters);

% file list
fileList = compile_file_table(experimentParameters, fileSelector);

numFiles = height(fileList);

for ff = 1:numFiles

    filepath = fileList.Filepath(ff);
    assert(isfile(filepath))
    reader = OMEImageReaderLazy(filepath);
    numSeries = reader.getNumSeries(); % Fast, uses metadata only

    fprintf("File %d of %d (%s, contains %d series)\n", ff, numFiles, filepath, numSeries)

    for ss = 1:numSeries

        metadataTable = append_metadata_entry(metadataTable1, fileList(ff,:), experimentParameters, reader, ss);

    end

    reader.close(); % Always close reader after use!

end

folder = fileparts(metadataFile);
if ~isfolder(folder)
    mkdir(folder);
end

writetable(metadataTable, metadataFile);
