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

% image data source directories
sourceDirectories = [ ...
	fullfile(".", "ImageData", "Folder_1"), ...
	fullfile(".", "ImageData", "Folder_2"), ...
	];

% condition labels (number must match source directories)
condLabels = [ ...
	"Example 1", ...
	"Example 2", ...
	];

% condition indexing (number must match source directories)
condInds = (1:numel(condLabels))';

% file name pattern to look for
fileSelector = "*.nd2";

% metadata output file
metadataFile = fullfile(".", "ExtractedStacks", "metadata_temp.csv");

%% Main script section

% metadata result table
tabVars = [ ...
    "FileName",    "string"; ...
    "CondName",    "string"; ...
    "CondInd",     "int32"; ...
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
metadataTable = table(Size=tabSize, VariableNames=tabVars(1,:), VariableTypes=tabVars(2,:));

numDirs = numel(sourceDirectories);

for cc = 1:numDirs

	thisDir = sourceDirectories(cc);

	fprintf("Directory %d of %d (%s)\n", cc, numDirs, thisDir)

	listing = dir(fullfile(thisDir, "**", fileSelector));

	numFiles = numel(listing);

    condName = condLabels(cc);
	condInd = condInds(cc);

	for ff = 1:numFiles

		filepath = fullfile(listing(ff).folder, listing(ff).name);

		fprintf("File %d of %d (%s)\n", ff, numFiles, filepath)

        reader = OMEImageReaderLazy(filepath);
        numSeries = reader.getNumSeries(); % Fast, uses metadata only

        for ss = 1:numSeries
            
            fprintf("Series %d of %d\n", ss, numSeries)

            numChannels = reader.getNumChannels(ss);
            imgSize = reader.getStackSizeXYZ(ss);
            numTime = reader.getSizeT(ss);
		 	voxelSize = [reader.getPixelSizeXY(ss), reader.getZStepSize(ss)];

            metadataTable(end + 1, :) = { ...
                string(filepath), ...   % FileName
                string(condName), ...   % CondName
                condInd, ...            % CondInd
                numSeries, ...          % SeriesTotal
                ss, ...                 % SeriesInd
                numChannels, ...        % SizeC
                imgSize(1), ...         % SizeX
                imgSize(2), ...         % SizeY
                imgSize(3), ...         % SizeZ
                numTime, ...            % SizeT
                voxelSize(1), ...       % VoxelSizeX
                voxelSize(2), ...       % VoxelSizeY
                voxelSize(3), ...       % VoxelSizeZ
                };

        end
        reader.close(); % Always close reader after use!

	end

end

folder = fileparts(metadataFile);
if ~isfolder(folder)
    mkdir(folder);
end

writetable(metadataTable, metadataFile);
