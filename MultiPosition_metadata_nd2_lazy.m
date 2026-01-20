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

% metadata output file
metadataFile = fullfile(".", "ExtractedStacks", "metadata_temp.csv");

%% Main script section

% metadata result table
metadataTable = table(Size=[0,13], ...
    VariableNames=["FileName","CondName","CondInd","SeriesTotal","SeriesInd","SizeC","SizeX","SizeY","SizeZ","SizeT","VoxelSizeX","VoxelSizeY","VoxelSizeZ"], ...
    VariableTypes=["string","string","int32","int32","int32","int32","int32","int32","int32","int32","double","double","double"]);

numDirs = numel(sourceDirectories);

for cc = 1:numDirs

	thisDir = sourceDirectories(cc);

	fprintf("Directory %d of %d (%s)\n", cc, numDirs, thisDir)

	listing = rdir(char(fullfile(thisDir, "*.nd2")), '~contains(name,''._'')');

	numFiles = numel(listing);

    condName = condLabels(cc);
	condInd = condInds(cc);

	for ff = 1:numFiles

		combined_filepath = fullfile(listing(ff).name);

		fprintf("File %d of %d (%s)\n", ff, numFiles, combined_filepath)

        reader = OMEImageReaderLazy(combined_filepath);
        numSeries = reader.getNumSeries(); % Fast, uses metadata only

        for ss = 1:numSeries
            
            fprintf("Series %d of %d\n", ss, numSeries)

            numChannels = reader.getNumChannels(ss);
            imgSize = reader.getStackSizeXYZ(ss);
            numTime = reader.getSizeT(ss);
		 	voxelSize = [reader.getPixelSizeXY(ss), reader.getZStepSize(ss)];

            metadataTable(end + 1, :) = {string(combined_filepath), string(condName), condInd, numSeries, ss, numChannels, imgSize(1), imgSize(2), imgSize(3), numTime, voxelSize(1), voxelSize(2), voxelSize(3)};

        end
        reader.close(); % Always close reader after use!

	end

end

folder = fileparts(metadataFile);
if ~isfolder(folder)
    mkdir(folder);
end

writetable(metadataTable, metadataFile);
