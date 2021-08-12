clear all

plotFlag = false;

sourceDirectories = { ...
	'./LIF/LIF_Ctrl24/',...
	'./LIF/LIF_Ctrl48/',...
	'./LIF/LIF_3/',...
	'./LIF/LIF_6/',...
	'./LIF/LIF_12/',...
	'./LIF/LIF_24/',...
	'./LIF/LIF_48/',...
	'./RHB/RHB_Ctrl24/'...
	'./RHB/RHB_Ctrl48/'...
	'./RHB/RHB_3/'...
	'./RHB/RHB_6/'...
	'./RHB/RHB_12/'...
	'./RHB/RHB_24/'...
	'./RHB/RHB_48/'...
	};
extractTargetFolder = 'ExtractedStacks';
condInds = (1:14)';
condLabels = {...
	'-LiF, Ctrl 24 h','-LiF, Ctrl 48 h',...
	'-LIF, 3 h','-LIF, 6 h','-LIF, 12 h','-LIF, 24 h','-LIF, 48 h',...
	'RHB, Ctrl 24 h','RHB, Ctrl 48 h',...
	'RHB, 3 h','RHB, 6 h','RHB, 12 h','RHB, 24 h','RHB, 48 h',...
	};
maxNumSeries = Inf;

skipList = []; % Directories to skip, for example if already done
% Leave empty array [] if all directories shoudl be processed.

useChannel_inds = [4,1];
numChannels = numel(useChannel_inds);
scaleChannels = {[-Inf,Inf],[-Inf,Inf]};

numDirs = numel(sourceDirectories);

for cc = 1:numDirs
    
	fprintf('Extracting image data from directory %d of %d\n',cc,numDirs)
	
	if ismember(cc,skipList)
		
		fprintf('On skip list, skipping to next directory.\n')
		
	else
		
		thisDir = sourceDirectories{cc};
		
		listing = rdir([thisDir,'*.nd2'],'~contains(name,''._'')');
		
		numFiles = numel(listing);
		
		imgCell = cell(numFiles,1);
		imgSizeCell = cell(numFiles,1);
		pixelSizeVec = cell(numFiles,1);
		zzStepVec = cell(numFiles,1);
		
		stackCounter = 0;
		mkdir(sprintf('./%s/Cond_%d',extractTargetFolder,cc))
		
		for ff = 1:numFiles
			
			fprintf('Extracting images from file %d of %d\n',ff,numFiles)
			
			combined_filepath = ...
				fullfile(listing(ff).name);
			
			[imgCell,imgSizeCell,...
				this_voxelSize_cell,seriesName_cell] = ...
				nd2read(combined_filepath,maxNumSeries);
			
			pixelSizeVec = cellfun(@(elmt)elmt(1),this_voxelSize_cell);
			zzStepVec = cellfun(@(elmt)elmt(3),this_voxelSize_cell);
			
			numImgs = numel(imgCell);
			for kk = 1:numImgs
				imgCell{kk} = imgCell{kk}(useChannel_inds);
			end
			
			if plotFlag
				
				for kk = 1:numImgs
					
					thisPixelSize = pixelSizeVec(kk);
					thisSize = imgSizeCell{kk};
					
					numChannels = numel(imgCell{kk});
					
					for cc_2 = 1:numChannels
						
						subplot(1,numChannels,cc_2)
						
						imagesc([0,thisSize(2)].*thisPixelSize,...
							[0,thisSize(1)].*thisPixelSize,...
							imgCell{kk}{cc_2}(:,:,ceil(thisSize(3)./2)),...
							scaleChannels{1})
						axis equal tight
						
						xlabel('x [\mum]')
						
						colormap((gray))
						
					end
					
					disp(sprintf( ...
						'Directory number %d, Condition %s, File %d/%d, Image %d/%d.', ...
						cc,condLabels{condInds(cc)},ff,numFiles,kk,numImgs))
					
					waitforbuttonpress
					
				end
				
			end
			
			% Apply exclusion indices
			keepInds = 1:numImgs;
			imgCell = imgCell(keepInds);
			imgSizeCell = imgSizeCell(keepInds);
			pixelSizeVec = pixelSizeVec(keepInds);
			zzStepVec = zzStepVec(keepInds);
			
			% --- saving of the retained images
			
			condName = condLabels{cc};
			condInd = condInds(cc);
			
			% Save all images into a separate file, but contained in a folder that
			% contains all images for this condition
			
			numImgs = numel(imgCell);
			for kk = 1:numImgs
				stackCounter = stackCounter+1;
				imgStack = imgCell{kk};
				imgSize = imgSizeCell{kk};
				pixelSize = pixelSizeVec(kk);
				zStepSize = zzStepVec(kk);
				save(sprintf('./%s/Cond_%d/Image_%d.mat',...
					extractTargetFolder,cc,stackCounter),...
					'imgStack','imgSize','pixelSize','zStepSize',...
					'condInd','condName')
			end
			
		end
		
	end
	
end