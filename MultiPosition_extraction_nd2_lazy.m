%% Extract experiment image data to MATLAB files
%
% This script extracts the specified image data from the experiment image data
% files and exports separate (multi-channel) stacks to MATLAB files for
% convenient access from subsequent MATLAB analysis scripts.
%
% The channels to be extracted can be specified. Only the required image data
% are actually read from the source file.
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

clear all

%% Process parameter section

% switch figure plots on/off (central xy section, all channels, pauses script)
plotFlag = true;

% image data source directories
sourceDirectories = { ...
    fullfile('.', 'ImageData', 'Folder_1/'), ...
    fullfile('.', 'ImageData', 'Folder_2/'), ...
    };

% condition labels (number must match source directories)
condLabels = {...
    'Example 1', ...
    'Example 2', ...
    };

% condition indexing (number must match source directories)
condInds = (1:numel(condLabels))';

% limit the number of extracted stacks per file
maxNumSeries = Inf;

% extracted data output root folder
extractTargetFolder = fullfile('.', 'ExtractedStacks');

% sprintf name pattern for output sub-folder (one per experiment condition)
% (make number pattern wide enough to accomodate for all conditions, to ensure
% correct sorting in all OS environments)
condFolderPattern = 'Cond_%02d';

% sprintf name pattern for output files (one per stack)
% (make number pattern wide enough to accomodate for all stacks per condition,
% to ensure correct sorting in all OS environments)
imageFilePattern = 'Image_%03d.mat';

% directories to skip, for example if already done
% (leave empty array [] if all directories should be processed)
skipList = [];

% channels to be extracted
useChannel_inds = [7,5,4,2]; % Hoechst, Alexa 488, Alexa 594, STAR Red

%% Main script section

numChannels = numel(useChannel_inds);
numDirs = numel(sourceDirectories);

for cc = 1:numDirs
    
	fprintf('Extracting image data from directory %d of %d\n',cc,numDirs)
	
	if ismember(cc,skipList)
		
		fprintf('On skip list, skipping to next directory.\n')
		
	else
		
		thisDir = sourceDirectories{cc};
		
		listing = rdir([thisDir,'*.nd2'],'~contains(name,''._'')');
		
		numFiles = numel(listing);

        condName = condLabels{cc};
    	condInd = condInds(cc);
        
		stackCounter = 0;
		condSubFolder = sprintf(condFolderPattern,cc);
		mkdir(sprintf('./%s/%s',extractTargetFolder,condSubFolder))
		
		for ff = 1:numFiles
			
			fprintf('Extracting images from file %d of %d\n',ff,numFiles)
			
			combined_filepath = ...
				fullfile(listing(ff).name);

            reader = OMEImageReaderLazy(combined_filepath);
            nSeries = reader.getNumSeries(); % Fast, uses metadata only
            numSeries = min(nSeries,maxNumSeries);
            
            for ss = 1:numSeries
                
                % read in the series
                imgStack = reader.getChannelStack(ss,useChannel_inds); % Loads channel 1 and 3 images only, on demand
                fprintf('Image series %d of %d extracted\n',ss,numSeries)

                % save the series to an .mat file
                stackCounter = stackCounter+1;
			 	pixelSize = reader.getPixelSizeXY(ss);
                pixelSize = pixelSize(1);
				zStepSize = reader.getZStepSize(ss);
                imgSize = reader.getStackSizeXYZ(ss);
				imgFile = sprintf(imageFilePattern,stackCounter);
				save(sprintf('./%s/%s/%s',...
					extractTargetFolder,condSubFolder,imgFile),...
					'imgStack','imgSize','pixelSize','zStepSize',...
					'condInd','condName','-v7.3')

                if plotFlag

                    numChannels = numel(imgStack);

                    for cc_2 = 1:numChannels

                        subplot(1,numChannels,cc_2)

                        imagesc([0,imgSize(1)].*pixelSize,...
                            [0,imgSize(2)].*pixelSize,...
                            imgStack{cc_2}(:,:,ceil(imgSize(3)./2)),...
                            prctile(imgStack{cc_2}(:),[0.1,99.9]))
                        axis equal tight

                        xlabel('x [\mum]')

                        colormap((gray))

                    end

                    fprintf( ...
                        'Directory number %d, Condition %s, File %d/%d, Image %d/%d.\n', ...
                        cc,condLabels{condInds(cc)},ff,numFiles,ss,numSeries);

                    waitforbuttonpress

                end

            end
            reader.close(); % Always close reader after use!

		end
		
	end
	
end
