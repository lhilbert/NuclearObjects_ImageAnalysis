clear all

plotFlag = true;

sourceDirectories = { ...
    './ImageData/Folder_2/',...
    './ImageData/Folder_1/',...
    };
extractTargetFolder = './ExtractedStacks';
condLabels = {...
    'Example 2',...
    'Example 1',...
    };
condInds = (1:numel(condLabels))';
maxNumSeries = 3;

condFolderPattern = 'Cond_%02d';     % sprintf pattern used for output subfolders
imageFilePattern = 'Image_%03d.mat'; % sprintf pattern used for output files
% Make the appended number pattern wide enough to accomodate for all conditions
% and images per condition (to ensure correct sorting in all environments)

skipList = []; % Directories to skip, for example if already done
% Leave empty array [] if all directories should be processed.
useChannel_inds = [7,5,4,2]; % Hoechst, Alexa 488, Alexa 594, STAR Red

numChannels = numel(useChannel_inds);
scaleChannels = {[-Inf,Inf],[-Inf,Inf],[-Inf,Inf],[-Inf,Inf]};

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

                    disp(sprintf( ...
                        'Directory number %d, Condition %s, File %d/%d, Image %d/%d.', ...
                        cc,condLabels{condInds(cc)},ff,numFiles,ss,numSeries))

                    waitforbuttonpress

                end

            end
            reader.close(); % Always close reader after use!

		end
		
	end
	
end
