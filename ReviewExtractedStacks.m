%% Review extracted experiment image stacks
%
% This script loads image data stacks one by one and shows the central xy
% section (from each extracted channel) for visual review. It works on
% single-stack MATLAB files that were previously extracted from the original
% source files.
%
% The section images can either be reviewed during script execution, or they can
% be created and saved to a separate folder for later review. The review images
% can be useful for sanity-checking the extracted data.

clear all

%% Process parameter section

% switch figure plot export on/off (script pauses if off)
batchExportFlag = true;

% source directory containing the extraced single-stack MATLAB files
sourceDirectory = fullfile('.', 'ExtractedStacks', '**');

% file name pattern for the single-stack MATLAB files
imageFileSelector = '*Image*.mat';

% review plot output root folder (sub-folders are created as in image folder)
exportDirectory = fullfile('.', 'ReviewPlots');

%% Main script section

% --- file read and display procedure begins here

listing = rdir(fullfile(sourceDirectory, '*Image*.mat'));
numFiles = numel(listing);

% Condition index retrieval
condInds = zeros(1, numFiles);
condNames = cell(1, numFiles);
for ff = 1:numFiles
	thisFilePath = listing(ff).name;
	thisCondInd = load(thisFilePath,'condInd');
	condInds(ff) = thisCondInd.condInd;
	thisCondName = load(thisFilePath,'condName');
	condNames{ff} = thisCondName.condName;
end

% --- analyze image stacks one by one

% Variables to store properties of nuclei
numNuclei_vec = zeros(1,numFiles);
nuc_intCell = cell(1,numFiles);
cyto_intCell = cell(1,numFiles);
nuc_medianVolCell = cell(1,numFiles);
perNuc_countCell = cell(1,numFiles);

% Variables to store properties of objects inside nuclei
S5P_volCell = cell(1,numFiles);
S5P_solCell = cell(1,numFiles);
S5P_eloCell = cell(1,numFiles);
S5P_intCell = cell(1,numFiles);
S5P_centCell = cell(1,numFiles);
S5P_imgCell = cell(1,numFiles);

S2P_volCell = cell(1,numFiles);
S2P_solCell = cell(1,numFiles);
S2P_eloCell = cell(1,numFiles);
S2P_intCell = cell(1,numFiles);
S2P_centCell = cell(1,numFiles);
S2P_imgCell = cell(1,numFiles);

for ff = 1:numFiles
	
	fprintf('Processing file %d of %d\n',ff,numFiles)
	
	thisCondInd = condInds(ff);	
    thisCondName = condNames{ff};
	thisFilePath = listing(ff).name;
	
    fprintf('File name: %s\n',thisFilePath)

    clear imgStack imgSize pixelSize zStepSize;
	load(thisFilePath,'imgStack','imgSize','pixelSize','zStepSize');
    centralSlice = round(imgSize(3)./2);
    numChannels = numel(imgStack);
	
    for nn = 1:numChannels
        
        subplot(1,numChannels,nn)
        imagesc([0,imgSize(2)].*pixelSize,...
            [0,imgSize(1)].*pixelSize,...
            squeeze(imgStack{nn}(:,:,centralSlice)))
        if nn == 2
            title(sprintf('File name: %s',thisFilePath),...
                'interpreter','none')
        elseif nn == 1
            title(sprintf('Condition name: %s',thisCondName),...
                'interpreter','none')
        end
        axis tight equal
        
    end
    
    if batchExportFlag

        [condSubFolder,exportFileName] = fileparts(listing(ff).name);
        [~,condSubFolder] = fileparts(condSubFolder);
        exportFolderPath = fullfile(exportDirectory,condSubFolder);
        if ~isfolder(exportFolderPath)
            mkdir(exportFolderPath);
        end
        exportFileName = sprintf('%s.pdf',exportFileName);
        exportFilePath = fullfile(exportFolderPath,exportFileName);
        exportgraphics(gcf,exportFilePath);
        exportFileName = replace(exportFileName,'pdf','png');
        exportFilePath = fullfile(exportFolderPath,exportFileName);
        exportgraphics(gcf,exportFilePath);

    else

        waitforbuttonpress

    end

end
