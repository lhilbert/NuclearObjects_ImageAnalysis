clear all

% Specify the directory that contains the extraced files from the last
% step, where you extracted from the raw files obtained from the microscope
sourceDirectory = './ExtractedStacks/Cond_1/';
sourceDirectory = './ExtractedStacks/**/';

% Channels for segmentation
NucSegChannel = 1; % Channel used to detect nuclei
S5P_SegChannel = 1; % Channel used to detect Pol II S5P clusters
S2P_SegChannel = 2; % Channel used to detect Pool II S2P clusters

% Save images of the clusters
ImgSquareSize = 0; % pixels for cut-out images, set 0 fr no images

% Target channels for intensity quantification, applied for all objects
quantChannels = [1,2];
quantBlurSigma = [0,0];
numQuantChannels = numel(quantChannels);

% Which image channels to store in example images
storeImgChannels = [];%[1,2];
numStoreChannels = numel(storeImgChannels);

nuc_segBlurSigma_small = 1.0; % in microns
nuc_segBlurSigma_large = 10; % in microns

S5P_segBlurSigma_small = 0.01;
S5P_segBlurSigma_large = 3.0;
S5P_seg_numStdDev = 2.0;

S2P_segBlurSigma_small = 0.01; % in microns
S2P_segBlurSigma_large = 5; % in microns
S2P_seg_numStdDev = 6; % number of standard deviations in robust threshold

% Minimum volume of nuclei, typical ranges for a full nucieus 10-100 cubic
% microns of volume, so set a cut-off oof 10 or 30 or so
Nuc_min_vol = 10; % cubic microns
Nuc_min_sol = 0.7; % to ensure round nuclei

% Inner and outer extension of a nuclear masks to cover cytoplasm
cytoMask_extension = 1.5; % in microns
cytoMask_distance = 1.0; % in microns

% Minimum volumes for objects inside the nuclei
S5P_minVol = 0.03; % cubic microns
S2P_minVol = 0.005; % cubic microns

% end of analysis parameter section, do not change anything else in
% this section, all necessary parameters are listed above

% --- analysis procedure begins here

numQuantChannels = numel(quantChannels);
listing = rdir([sourceDirectory,'*Image*.mat']);
numFiles = numel(listing);

% Condition index retrieval
condInds = [];
condNames = {};
for ff = 1:numFiles
	thisFilePath = listing(ff).name;
	thisCondInd = load(thisFilePath,'condInd');
	thisCondInd = thisCondInd.condInd;
	condInds = [condInds,thisCondInd];
	thisCondName = load(thisFilePath,'condName');
	thisCondName = thisCondName.condName;
	condNames = [condNames,thisCondName];
end

% Flag to log files that were successfully analyzed and contained nuclei
validFileFlag = false(1,numFiles);

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

parfor ff = 1:numFiles
	
	fprintf('Processing file %d of %d\n',ff,numFiles)
	
	thisCondInd = condInds(ff);	
	thisFilePath = listing(ff).name;
	
	loadStruct = load(thisFilePath,...
		'imgStack','imgSize','pixelSize','zStepSize');
	imgStack = loadStruct.imgStack;
	imgSize = loadStruct.imgSize;
	pixelSize = loadStruct.pixelSize;
	zStepSize = loadStruct.zStepSize;

	% Nuclei segmentation
	segImg = imgStack{NucSegChannel};
	segImg = ...
		+ imgaussfilt(segImg,nuc_segBlurSigma_small./pixelSize) ...
		- imgaussfilt(segImg,nuc_segBlurSigma_large./pixelSize);

	[bin_counts,bin_centers] = hist(segImg(:),1000);
	[nuc_seg_thresh,~] = otsuLimit(bin_centers,bin_counts,[0,Inf]);
	NucSegMask = segImg>1.0.*nuc_seg_thresh;
	
	subplot(1,3,1)
	imagesc(squeeze(imgStack{NucSegChannel}(:,:,ceil(imgSize(3)./2))))
	axis tight equal
	
	subplot(1,3,2)
	imagesc(squeeze(segImg(:,:,ceil(imgSize(3)./2))))
	axis tight equal

	subplot(1,3,3)
	imagesc(squeeze(NucSegMask(:,:,ceil(imgSize(3)./2))))
	axis tight equal
	
% Uncomment the following two lines, and remove the par in parfor above, if
% you want to check the extracted images one by one
%   fprintf('File name: %s',thisFilePath)
% 	waitforbuttonpress
	
	% --- Connected component segmentation of nuclei
	comps = bwconncomp(NucSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	minPixels = Nuc_min_vol./(pixelSize.^2)./zStepSize;
	comps.NumObjects = sum(numPxls>=minPixels);
	comps.PixelIdxList = comps.PixelIdxList(numPxls>=minPixels);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	
	props = regionprops3(comps,imgStack{NucSegChannel},...
		'Solidity');
	
	Solidity_array = [props.Solidity];
	
	comps.NumObjects = sum(Solidity_array>=Nuc_min_sol);
	comps.PixelIdxList = comps.PixelIdxList(Solidity_array>=Nuc_min_sol);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	
	numNuclei = comps.NumObjects;
	numNuclei_vec(ff) = numNuclei;
	
	nuc_intCell{ff} = cell(1,numQuantChannels);
	cyto_intCell{ff} = cell(1,numQuantChannels);
	for qq = 1:numQuantChannels
		quantImg = imgStack{quantChannels(qq)};
		quantProps = regionprops3(comps,quantImg,...
			'MeanIntensity','VoxelIdxList');
		nuc_intCell{ff}{qq} = [quantProps.MeanIntensity];
		cyto_intCell{ff}{qq} = zeros(1,numNuclei);
		for nn = 1:numNuclei
			cytoMask = false(size(quantImg));
			cytoMask(quantProps.VoxelIdxList{nn}) = true;
			coreMask = cytoMask;
			se = strel('disk',round(cytoMask_extension./pixelSize));
			cytoMask = imdilate(cytoMask,se);
			se = strel('disk',round(cytoMask_distance./pixelSize));
			coreMask = imdilate(coreMask,se);
			cytoMask(coreMask) = false;
			cytoMask = cytoMask & ~NucSegMask;
			cyto_intCell{ff}{qq}(nn) = mean(quantImg(cytoMask));
		end
	end
	
	if comps.NumObjects>0
		
		validFileFlag(ff) = true;
		
		props = regionprops3(comps,imgStack{NucSegChannel},...
			'Volume','VoxelValues','Solidity','VoxelIdxList',...
			'BoundingBox');
		
		Volume_array = [props.Volume].*pixelSize.^2.*zStepSize;
		Intensity_array = cellfun(@(vals)median(vals),props.VoxelValues);
		Solidity_array = [props.Solidity];
		
		
		% --- For each nucleus, get objects from the different channels
		
		nuc_medianVolCell{ff} = cell(1,2);
		perNuc_countCell{ff} = cell(1,2);
		for qq = 1:2
			nuc_medianVolCell{ff}{qq} = zeros(numNuclei,1);
			perNuc_countCell{ff}{qq} = zeros(numNuclei,1);
		end
		
		S5P_volume = cell(1,numNuclei);
		S5P_solidity = cell(1,numNuclei);
		S5P_elongation = cell(1,numNuclei);
		S5P_centralSlices_store = cell(1,numNuclei);
		S5P_intensity = cell(numQuantChannels,numNuclei);
		S5P_cent_store = cell(1,numNuclei);
		
		S2P_volume = cell(1,numNuclei);
		S2P_solidity = cell(1,numNuclei);
		S2P_elongation = cell(1,numNuclei);
		S2P_centralSlices_store = cell(1,numNuclei);
		S2P_intensity = cell(numQuantChannels,numNuclei);
		S2P_cent_store = cell(1,numNuclei);
		
		for nn = 1:numNuclei
			
			boxArray = props.BoundingBox(nn,:);
			
			S5P_subImage = imgStack{S5P_SegChannel}(...
				boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
				boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
				boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
			
			S2P_subImage = imgStack{S2P_SegChannel}(...
				boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
				boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
				boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
			
			S5P_subImage = ...
				+ imgaussfilt(S5P_subImage,S5P_segBlurSigma_small./pixelSize) ...
				- imgaussfilt(S5P_subImage,S5P_segBlurSigma_large./pixelSize);
			S2P_subImage = ...
				+ imgaussfilt(S2P_subImage,S2P_segBlurSigma_small./pixelSize) ...
				- imgaussfilt(S2P_subImage,S2P_segBlurSigma_large./pixelSize);
			
			NucMask = false(imgSize);
			NucMask(props.VoxelIdxList{nn}) = true;
			NucMask_subImage = NucMask(...
				boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
				boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
				boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
			
			
			seg_intensities = S5P_subImage(NucMask_subImage);
			seg_mean = mean(seg_intensities);
			seg_std = std(seg_intensities);
			S5P_mask = (S5P_subImage.*NucMask_subImage)...
				>(seg_mean+S5P_seg_numStdDev.*seg_std);
			
			seg_intensities = S2P_subImage(NucMask_subImage);
			seg_mean = mean(seg_intensities);
			seg_std = std(seg_intensities);
			S2P_mask = (S2P_subImage.*NucMask_subImage)...
				>(seg_mean+S2P_seg_numStdDev.*seg_std);
			
			subImgSize = size(S5P_subImage);
			if numel(subImgSize)==2
				subImgSize(3)=1;
			end
			
			subplot(2,2,1)
			imagesc(squeeze(max(S5P_subImage,[],3)))
			axis tight equal
			subplot(2,2,2)
			imagesc(squeeze(max(S2P_subImage,[],3)))
			axis tight equal
			subplot(2,2,3)
			imagesc(squeeze(max(S5P_mask,[],3)))
			axis tight equal
			subplot(2,2,4)
			imagesc(squeeze(max(S2P_mask,[],3)))
			axis tight equal
			% 		waitforbuttonpress
			
			% For storage of example images
			if numStoreChannels > 0
				store_subImages = cell(1,numStoreChannels);
				for color = 1:numStoreChannels
					store_subImages{color} = imgStack{storeImgChannels(color)}(...
						boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
						boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
						boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
				end
			end
			
			S5P_comps = bwconncomp(S5P_mask,18);
			S5P_numPxls = cellfun(@(elmt)numel(elmt),S5P_comps.PixelIdxList);
			minPixels = S5P_minVol./(pixelSize.^2)./zStepSize;
			S5P_comps.NumObjects = sum(S5P_numPxls>minPixels);
			S5P_comps.PixelIdxList = S5P_comps.PixelIdxList(S5P_numPxls>minPixels);
			S5P_numPxls = cellfun(@(elmt)numel(elmt),S5P_comps.PixelIdxList);
			
			S2P_comps = bwconncomp(S2P_mask,18);
			S2P_numPxls = cellfun(@(elmt)numel(elmt),S2P_comps.PixelIdxList);
			minPixels = S2P_minVol./(pixelSize.^2)./zStepSize;
			S2P_comps.NumObjects = sum(S2P_numPxls>minPixels);
			S2P_comps.PixelIdxList = S2P_comps.PixelIdxList(S2P_numPxls>minPixels);
			S2P_numPxls = cellfun(@(elmt)numel(elmt),S2P_comps.PixelIdxList);
			
			
			if S5P_comps.NumObjects>0 && S2P_comps.NumObjects>0
				
				S2P_props = regionprops3(S2P_comps,S2P_subImage,...
					'Volume','Solidity',...
					'Centroid','Image','BoundingBox');
				
				S5P_props = regionprops3(S5P_comps,S5P_subImage,...
					'Volume','Solidity',...
					'Centroid','Image','BoundingBox');
				
				S5P_Volume_array = ...
					[S5P_props.Volume].*pixelSize.^2.*zStepSize;
				S5P_Solidity_array = [S5P_props.Solidity];
				S5P_Centroid_array = ...
					S5P_props.Centroid.*[pixelSize,pixelSize,zStepSize];
				
				S2P_Volume_array = ...
					[S2P_props.Volume].*pixelSize.^2.*zStepSize;
				S2P_Solidity_array = [S2P_props.Solidity];
				S2P_Centroid_array = ...
					S2P_props.Centroid.*[pixelSize,pixelSize,zStepSize];
				
				perNuc_countCell{ff}{1}(nn) = numel(S5P_Centroid_array);
				perNuc_countCell{ff}{2}(nn) = numel(S2P_Centroid_array);
				
				nuc_medianVolCell{ff}{1}(nn) = median(S5P_Volume_array);
				nuc_medianVolCell{ff}{2}(nn) = median(S2P_Volume_array);
				
				
				% --- get cluster central plane and elongation in-plane
				S5P_Elongation_array = ...
					zeros(size(S5P_Solidity_array));
				S5P_Slices_cell = ...
					cell(size(S5P_Solidity_array));
				
				for object_nn = 1:numel(S5P_Volume_array)
					
					boundingBox = S5P_props.BoundingBox(object_nn,:);
					thisImage = squeeze(S5P_subImage(...
						boundingBox(2)+0.5:boundingBox(2)+boundingBox(5)-0.5,...
						boundingBox(1)+0.5:boundingBox(1)+boundingBox(4)-0.5,...
						ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
					thisImage = thisImage-min(thisImage(:));
					thisImage = thisImage./max(thisImage(:));
					
					thisMask = S5P_props.Image(object_nn);
					if iscell(thisMask)
						% There is some weird inconsistent behavior here,
						% sometimes .Image(object_nn) results in a cell
						% output, sometimes in a matrix output. It is not
						% clear if there is a system to it.
						thisMask = thisMask{1};
					else
						thisMask = S5P_props.Image;
					end
					centerInd = ceil(size(thisMask,3)./2);
					thisMask = squeeze(thisMask(:,:,centerInd));
					thisImage((bwperim(thisMask))) = 0;
					S5P_Slices_cell{object_nn} = cell(1,numStoreChannels+2);
					for color = 1:numStoreChannels
						S5P_Slices_cell{object_nn}{color} = ...
							squeeze(store_subImages{color}(:,:,...
							ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
					end
					S5P_Slices_cell{object_nn}{numStoreChannels+1} = ...
						squeeze(S5P_mask(:,:,...
						ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
					S5P_Slices_cell{object_nn}{numStoreChannels+2} = ...
						squeeze(S2P_mask(:,:,...
						ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
					
					thisProps = regionprops(uint8(thisMask),...
						'MajorAxisLength','MinorAxisLength');
					S5P_Elongation_array(object_nn) = ...
						thisProps.MajorAxisLength./thisProps.MinorAxisLength;
				end
				S5P_volume{nn} = S5P_Volume_array;
				S5P_solidity{nn} = S5P_Solidity_array;
				S5P_elongation{nn} = S5P_Elongation_array;
				S5P_cent_store{nn} = S5P_Centroid_array;
				S5P_centralSlices_store{nn} = S5P_Slices_cell;
				
				
				
				S2P_Elongation_array = ...
					zeros(size(S2P_Solidity_array));
				S2P_Slices_cell = ...
					cell(size(S2P_Solidity_array));
				
				for object_nn = 1:numel(S2P_Volume_array)
					
					boundingBox = S2P_props.BoundingBox(object_nn,:);
					thisImage = squeeze(S2P_subImage(...
						boundingBox(2)+0.5:boundingBox(2)+boundingBox(5)-0.5,...
						boundingBox(1)+0.5:boundingBox(1)+boundingBox(4)-0.5,...
						ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
					thisImage = thisImage-min(thisImage(:));
					thisImage = thisImage./max(thisImage(:));
					
					thisMask = S2P_props.Image(object_nn);
					if iscell(thisMask)
						% There is some weird inconsistent behavior here,
						% sometimes .Image(object_nn) results in a cell
						% output, sometimes in a matrix output. It is not
						% clear if there is a system to it.
						thisMask = thisMask{1};
					else
						thisMask = S2P_props.Image;
					end
					centerInd = ceil(size(thisMask,3)./2);
					thisMask = squeeze(thisMask(:,:,centerInd));
					thisImage((bwperim(thisMask))) = 0;
					S2P_Slices_cell{object_nn} = cell(1,numStoreChannels+2);
					for color = 1:numStoreChannels
						S2P_Slices_cell{object_nn}{color} = ...
							squeeze(store_subImages{color}(:,:,...
							ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
					end
					S2P_Slices_cell{object_nn}{numStoreChannels+1} = ...
						squeeze(S2P_mask(:,:,...
						ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
					S2P_Slices_cell{object_nn}{numStoreChannels+2} = ...
						squeeze(S5P_mask(:,:,...
						ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
					
					thisProps = regionprops(uint8(thisMask),...
						'MajorAxisLength','MinorAxisLength');
					S2P_Elongation_array(object_nn) = ...
						thisProps.MajorAxisLength./thisProps.MinorAxisLength;
				end
				S2P_volume{nn} = S2P_Volume_array;
				S2P_solidity{nn} = S2P_Solidity_array;
				S2P_elongation{nn} = S2P_Elongation_array;
				S2P_cent_store{nn} = S2P_Centroid_array;
				S2P_centralSlices_store{nn} = S2P_Slices_cell;
				
				
				% --- quantification for all target channels
				for qq = 1:numQuantChannels
					
					channelInd = quantChannels(qq);
					quant_subImage = imgStack{channelInd}(...
						boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
						boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
						boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
					
					if quantBlurSigma(qq)>0
						quant_subImage = imgaussfilt(quant_subImage,...
							quantBlurSigma(qq)./pixelSize);
					end
					
					Quant_nucleusMedian = ...
						median(quant_subImage(NucMask_subImage));
					
					S5P_quant_props = regionprops3(...
						S5P_comps,quant_subImage,'VoxelValues');
					Quant_ClusterMedian = cellfun(...
						@(vals)median(vals),S5P_quant_props.VoxelValues);
					S5P_intensity{qq,nn} = ...
						Quant_ClusterMedian./Quant_nucleusMedian;
					
					S2P_quant_props = regionprops3(...
						S2P_comps,quant_subImage,'VoxelValues');
					Quant_ClusterMedian = cellfun(...
						@(vals)median(vals),S2P_quant_props.VoxelValues);
					S2P_intensity{qq,nn} = ...
						Quant_ClusterMedian./Quant_nucleusMedian;
					
				end
				
			else
				
				S5P_volume{nn} = [];
				S5P_solidity{nn} = [];
				S5P_elongation{nn} = [];
				S5P_centralSlices_store{nn} = {};
				S5P_cent_store{nn} = [];
				
				S2P_volume{nn} = [];
				S2P_solidity{nn} = [];
				S2P_elongation{nn} = [];
				S2P_centralSlices_store{nn} = {};
				S2P_cent_store{nn} = [];
				
				for qq = 1:numQuantChannels
					S5P_intensity{qq,nn} = [];
					S2P_intensity{qq,nn} = [];
				end
				
			end
			
		end
		
		S5P_volCell{ff} = vertcat(S5P_volume{:});
		S5P_solCell{ff} = vertcat(S5P_solidity{:});
		S5P_eloCell{ff} = vertcat(S5P_elongation{:});
		S5P_imgCell{ff} = vertcat(S5P_centralSlices_store{:});
		S5P_centCell{ff} = vertcat(S5P_cent_store{:});
		S5P_intCell{ff} = cell(1,numQuantChannels);
		for qq = 1:numQuantChannels
			S5P_intCell{ff}{qq} = vertcat(S5P_intensity{qq,:});
		end
		
		S2P_volCell{ff} = vertcat(S2P_volume{:});
		S2P_solCell{ff} = vertcat(S2P_solidity{:});
		S2P_eloCell{ff} = vertcat(S2P_elongation{:});
		S2P_imgCell{ff} = vertcat(S2P_centralSlices_store{:});
		S2P_centCell{ff} = vertcat(S2P_cent_store{:});
		S2P_intCell{ff} = cell(1,numQuantChannels);
		for qq = 1:numQuantChannels
			S2P_intCell{ff}{qq} = vertcat(S2P_intensity{qq,:});
		end	
	end
	
end

% Retain only files that returned nuclei

condNames = condNames(validFileFlag);
nuc_intCell = nuc_intCell(validFileFlag);
cyto_intCell = cyto_intCell(validFileFlag);
nuc_medianVolCell = nuc_medianVolCell(validFileFlag);
perNuc_countCell = perNuc_countCell(validFileFlag);

S5P_volCell = S5P_volCell(validFileFlag);
S5P_solCell = S5P_solCell(validFileFlag);
S5P_eloCell = S5P_eloCell(validFileFlag);
S5P_imgCell = S5P_imgCell(validFileFlag);
S5P_centCell = S5P_centCell(validFileFlag);
S5P_intCell = S5P_intCell(validFileFlag);

S2P_volCell = S2P_volCell(validFileFlag);
S2P_solCell = S2P_solCell(validFileFlag);
S2P_eloCell = S2P_eloCell(validFileFlag);
S2P_imgCell = S2P_imgCell(validFileFlag);
S2P_centCell = S2P_centCell(validFileFlag);
S2P_intCell = S2P_intCell(validFileFlag);


%% Sort into conditions

uniqueCondNames = unique(condNames);
numConds = numel(uniqueCondNames);
fileIndsCell = cell(1,numConds);
numFiles_perCond = zeros(1,numConds);
for cc = 1:numConds
	fileIndsCell{cc} = cellfun(...
		@(elmt)strcmp(elmt,uniqueCondNames{cc}),condNames);
	numFiles_perCond(cc) = sum(fileIndsCell{cc});
end

sortedCondNames = cell(1,numConds);
sortedNumNuclei = zeros(1,numConds);
sortedNumFiles = zeros(1,numConds);

sortedNucIntCell = cell(1,numQuantChannels);
sortedCytoIntCell = cell(1,numQuantChannels);

sortedS5PNumCell = cell(1,numConds);
sortedS5PVolCell = cell(1,numConds);
sortedS5PSolCell = cell(1,numConds);
sortedS5PEloCell = cell(1,numConds);
sortedS5PCentralSliceCell = cell(1,numConds);
sortedS5PCentroidsCell = cell(1,numConds);
sortedS5PIntCell = cell(1,numQuantChannels);

sortedS2PNumCell = cell(1,numConds);
sortedS2PVolCell = cell(1,numConds);
sortedS2PSolCell = cell(1,numConds);
sortedS2PEloCell = cell(1,numConds);
sortedS2PCentralSliceCell = cell(1,numConds);
sortedS2PCentroidsCell = cell(1,numConds);
sortedS2PIntCell = cell(1,numQuantChannels);

for qq = 1:numQuantChannels
	sortedNucIntCell{qq} = cell(1,numConds);
	sortedCytoIntCell{qq} = cell(1,numConds);
	sortedS5PIntCell{qq} = cell(1,numConds);
	sortedS2PIntCell{qq} = cell(1,numConds);
end


for cc = 1:numConds
	
	sortedCondNames{cc} = ...
		condNames(fileIndsCell{cc});
	sortedCondNames{cc} = sortedCondNames{cc}{1};
	
	sortedNumNuclei(cc) = ...
		sum(numNuclei_vec(fileIndsCell{cc}));
	sortedNumFiles(cc) = sum(fileIndsCell{cc});
	
	S5P_nums = vertcat(arrayfun(...
		@(val)perNuc_countCell{val}{1},find(fileIndsCell{cc}),...
		'UniformOutput',false));
	S5P_nums = vertcat(S5P_nums{:});
	S5P_vols = vertcat(S5P_volCell{fileIndsCell{cc}});
	S5P_sols = vertcat(S5P_solCell{fileIndsCell{cc}});
	S5P_elos = vertcat(S5P_eloCell{fileIndsCell{cc}});
	S5P_slices = vertcat(S5P_imgCell{fileIndsCell{cc}});
	S5P_centroids = vertcat(S5P_centCell{fileIndsCell{cc}});
	S5P_ints = S5P_intCell(fileIndsCell{cc});
	
	sortedS5PNumCell{cc} = S5P_nums;
	sortedS5PVolCell{cc} = S5P_vols;
	sortedS5PSolCell{cc} = S5P_sols;
	sortedS5PEloCell{cc} = S5P_elos;
	sortedS5PCentralSliceCell{cc} = S5P_slices;
	sortedS5PCentroidsCell{cc} = S5P_centroids;


	S2P_nums = vertcat(arrayfun(...
		@(val)perNuc_countCell{val}{2},find(fileIndsCell{cc}),...
		'UniformOutput',false));
	S2P_nums = vertcat(S2P_nums{:});
	S2P_vols = vertcat(S2P_volCell{fileIndsCell{cc}});
	S2P_sols = vertcat(S2P_solCell{fileIndsCell{cc}});
	S2P_elos = vertcat(S2P_eloCell{fileIndsCell{cc}});
	S2P_slices = vertcat(S2P_imgCell{fileIndsCell{cc}});
	S2P_centroids = vertcat(S2P_centCell{fileIndsCell{cc}});
	S2P_ints = S2P_intCell(fileIndsCell{cc});

	sortedS2PNumCell{cc} = S2P_nums;
	sortedS2PVolCell{cc} = S2P_vols;
	sortedS2PSolCell{cc} = S2P_sols;
	sortedS2PEloCell{cc} = S2P_elos;
	sortedS2PCentralSliceCell{cc} = S2P_slices;
	sortedS2PCentroidsCell{cc} = S2P_centroids;

	
	for qq = 1:numQuantChannels
		
		sortedNucIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)nuc_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedNucIntCell{qq}{cc} = horzcat(sortedNucIntCell{qq}{cc}{:})';
		
		sortedCytoIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)cyto_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedCytoIntCell{qq}{cc} = horzcat(sortedCytoIntCell{qq}{cc}{:})';
		
		sortedS5PIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)S5P_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedS5PIntCell{qq}{cc} = vertcat(sortedS5PIntCell{qq}{cc}{:});
		
		sortedS2PIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)S2P_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedS2PIntCell{qq}{cc} = vertcat(sortedS2PIntCell{qq}{cc}{:});
		
	end
	
end


%% Plotting of result overview

sorted_central_slices = cell(1,numConds);

Num_median = zeros(1,numConds);
Num_CI = zeros(2,numConds);
S5P_median = zeros(1,numConds);
S5P_CI = zeros(2,numConds);
S2P_median = zeros(1,numConds);
S2P_CI = zeros(2,numConds);
Vol_median = zeros(1,numConds);
Vol_CI = zeros(2,numConds);
Elo_median = zeros(1,numConds);
Elo_CI = zeros(2,numConds);
Sol_median = zeros(1,numConds);
Sol_CI = zeros(2,numConds);

for cc = 1:numConds
	
	Vol_threshold = 0.01;0.25;	
	S5P_threshold = 0;1.2;
	inclInds = ...
		sortedS5PVolCell{cc}>=Vol_threshold ...
		& sortedS5PIntCell{1}{cc}>=S5P_threshold;
	
	S5P_Num_vals = [sortedS5PNumCell{cc}];
	S2P_Num_vals = [sortedS2PNumCell{cc}];
	Nuc_S5P_vals = [sortedNucIntCell{1}{cc}-sortedCytoIntCell{1}{cc}];
	Nuc_S2P_vals = [sortedNucIntCell{2}{cc}-sortedCytoIntCell{1}{cc}];
	
	S5P_S5P_vals = [sortedS5PIntCell{1}{cc}(inclInds)];
	S5P_S2P_vals = [sortedS5PIntCell{2}{cc}(inclInds)];
	S5P_Vol_vals = [sortedS5PVolCell{cc}(inclInds)];
    S5P_Elo_vals = [sortedS5PEloCell{cc}(inclInds)];
	S5P_Sol_vals = [sortedS5PSolCell{cc}(inclInds)];
% 	S5P_slices = [sortedS5PCentralSliceCell{cc}(inclInds)];
	
	S2P_S5P_vals = [sortedS2PIntCell{1}{cc}];
	S2P_S2P_vals = [sortedS2PIntCell{2}{cc}];
	S2P_Vol_vals = [sortedS2PVolCell{cc}];
    S2P_Elo_vals = [sortedS2PEloCell{cc}];
	S2P_Sol_vals = [sortedS2PSolCell{cc}];
% 	S2P_slices = [sortedS2PCentralSliceCell{cc}];
		
	numPoints = numel(S5P_S2P_vals);
	
	n_boot = 200;
	S5P_median(cc) = median(Nuc_S5P_vals);
	S5P_CI(:,cc) = bootci(n_boot,@median,Nuc_S5P_vals);
	S2P_median(cc) = median(Nuc_S2P_vals);
	S2P_CI(:,cc) = bootci(n_boot,@median,Nuc_S2P_vals);
% 	S5P_median(cc) = median(S5P_S5P_vals);
% 	S5P_CI(:,cc) = bootci(n_boot,@median,S5P_S5P_vals);
% 	S2P_median(cc) = median(S5P_S2P_vals);
% 	S2P_CI(:,cc) = bootci(n_boot,@median,S5P_S2P_vals);
% 	S5P_median(cc) = median(S2P_S5P_vals);
% 	S5P_CI(:,cc) = bootci(n_boot,@median,S2P_S5P_vals);
% 	S2P_median(cc) = median(S2P_S2P_vals);
% 	S2P_CI(:,cc) = bootci(n_boot,@median,S2P_S2P_vals);
	Vol_median(cc) = median(S5P_Vol_vals);
	Vol_CI(:,cc) = bootci(n_boot,@median,S5P_Vol_vals);
	Elo_median(cc) = median(S5P_Elo_vals);
	Elo_CI(:,cc) = bootci(n_boot,@median,S5P_Elo_vals);
	Sol_median(cc) = median(S5P_Sol_vals);
	Sol_CI(:,cc) = bootci(n_boot,@median,S5P_Sol_vals);
	Num_median(cc) = median(S5P_Num_vals);
	Num_CI(:,cc) = bootci(n_boot,@median,S5P_Num_vals);
	
% 	subplot(subplot_layout(1),subplot_layout(2),...
% 		plot_order(cc))
% 	
% 	plot([1,1].*0.5,[0,1],'k-','Color',[0.6,0.6,0.6],'LineWidth',1)
% 	hold on
% 	plot([0,4],[1,1].*0.35,'k-','Color',[0.6,0.6,0.6],'LineWidth',1)
% 	plot((S5P_Vol_vals),S5P_Sol_vals,'k.')
% 	set(gca,'XLim',[0,4],'YLim',[0,1])
% 	title(sortedCondNames{cc})
% 	xlabel('Volume [\mum^3]')
% 	ylabel('Solidity')
% % 	waitforbuttonpress
	
end

%% plot overview

figure(1)
clf

datasetInds = {[1,2],[3,4]};
xaxisInds = {[1,3],...
	[2,4]};
numPlotSets = numel(xaxisInds);
setSymbols = {'ko','rs'};
setFaceColors = {[0,0,0],[1,0,0]};
setNames = {'AB combination 1','AB combination 2'};

numLabels = sum(cellfun(@(elmt)numel(elmt),xaxisInds));
labelCallInds = [xaxisInds{:}];

for pp = 1:numPlotSets
	
	subplot(5,1,1)
	errorbar(xaxisInds{pp},S5P_median(datasetInds{pp}),...
		S5P_CI(1,(datasetInds{pp}))-S5P_median(datasetInds{pp}),...
		S5P_median(datasetInds{pp})-S5P_CI(2,(datasetInds{pp})),...
		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
	hold on
	
	set(gca,'XTick',1:numLabels,...
		'XTickLabels',sortedCondNames(labelCallInds))
	xlabel('')
	ylabel('Pol II Ser5P Int.')
	xtickangle(30)
	
	subplot(5,1,2)
	errorbar(xaxisInds{pp},S2P_median(datasetInds{pp}),...
		S2P_CI(1,(datasetInds{pp}))-S2P_median(datasetInds{pp}),...
		S2P_median(datasetInds{pp})-S2P_CI(2,(datasetInds{pp})),...
		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
	hold on
	set(gca,'XTick',1:numLabels,...
		'XTickLabels',sortedCondNames(labelCallInds))
	xlabel('')
	ylabel('Pol II Ser2P Int.')
	xtickangle(30)
	
	subplot(5,1,3)
	errorbar(xaxisInds{pp},Num_median(datasetInds{pp}),...
		Num_CI(1,(datasetInds{pp}))-Num_median(datasetInds{pp}),...
		Num_median(datasetInds{pp})-Num_CI(2,(datasetInds{pp})),...
		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
	hold on
	set(gca,'XTick',1:numLabels,...
		'XTickLabels',sortedCondNames(labelCallInds))
	xlabel('')
	ylabel('No. clusters')
	xtickangle(30)
	
	subplot(5,1,4)
	errorbar(xaxisInds{pp},Vol_median(datasetInds{pp}),...
		Vol_CI(1,(datasetInds{pp}))-Vol_median(datasetInds{pp}),...
		Vol_median(datasetInds{pp})-Vol_CI(2,(datasetInds{pp})),...
		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
	hold on
	set(gca,'XTick',1:numLabels,...
		'XTickLabels',sortedCondNames(labelCallInds))
	xlabel('')
	ylabel('Volume [\mum^3]')
	xtickangle(30)
	
% 	subplot(5,1,4)
% 	errorbar(xaxisInds{pp},Elo_median(datasetInds{pp}),...
% 		Elo_CI(1,(datasetInds{pp}))-Elo_median(datasetInds{pp}),...
% 		Elo_median(datasetInds{pp})-Elo_CI(2,(datasetInds{pp})),...
% 		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
% 	hold on
% 	set(gca,'XTick',1:numConds,...
% 		'XTickLabels',sortedCondNames(xaxisInds{pp}))
% 	xlabel('')
% 	ylabel('Elongation')
% 	xtickangle(30)
	
	subplot(5,1,5)
	errorbar(xaxisInds{pp},Sol_median(datasetInds{pp}),...
		Sol_CI(1,(datasetInds{pp}))-Sol_median(datasetInds{pp}),...
		Sol_median(datasetInds{pp})-Sol_CI(2,(datasetInds{pp})),...
		setSymbols{pp},'MarkerFaceColor',setFaceColors{pp})
	hold on
	set(gca,'XTick',1:numLabels,...
		'XTickLabels',sortedCondNames(labelCallInds))
	xlabel('')
	ylabel('Solidity')
	xtickangle(30)
	
end

subplot(5,1,1)

legend(setNames,'Location','Northwest')


% %% --- plot the cluster masks
% 
% figure(2)
% clf
% 
% loadStruct = load(listing(1).name,'pixelSize');
% pixelSize = loadStruct.pixelSize;
% 
% % Condition index, Cluster index, zoom switch
% exampleImgs = {...
% 	[1,10],[1,11],[1,12],[1,13],[1,14],...
% 	[1,15],[1,16],[1,17],[1,18],[1,19]};
% 
% numExamples = numel(exampleImgs);
% 
% rescaleFac = 2.5;
% windowSize = 2.2; %microns
% 
% for ee = 1:numExamples
% 	
% 	OP_img = imresize(...
% 		sortedS5PCentralSliceCell{exampleImgs{ee}(1)}{exampleImgs{ee}(2)}{1},...
% 		rescaleFac);
% 	S5P_img = imresize(...
% 		sortedS5PCentralSliceCell{exampleImgs{ee}(1)}{exampleImgs{ee}(2)}{2},...
% 		rescaleFac);
% 	OP_img = OP_img - prctile(OP_img(:),5);
% 	OP_img(OP_img<0) = 0;
% 	OP_img = OP_img./prctile(OP_img(:),99.999);
% 	OP_img(OP_img>1) = 1;
% 	S5P_img = S5P_img - prctile(S5P_img(:),50);
% 	S5P_img(S5P_img<0) = 0;
% 	S5P_img = S5P_img./prctile(S5P_img(:),99.99);
% 	S5P_img(S5P_img>1) = 1;
% 	S2P_img = imresize(...
% 		sortedS5PCentralSliceCell{exampleImgs{ee}(1)}{exampleImgs{ee}(2)}{3},...
% 		rescaleFac);
% 	S2P_img = S2P_img - prctile(S2P_img(:),20);
% 	S2P_img(S2P_img<0) = 0;
% 	S2P_img = S2P_img./prctile(S2P_img(:),99.9);
% 	S2P_img(S2P_img>1) = 1;
% 	S5P_mask = imresize(...
% 		sortedS5PCentralSliceCell{exampleImgs{ee}(1)}{exampleImgs{ee}(2)}{4},...
% 		rescaleFac);
% 	OP_mask = imresize(...
% 		sortedS5PCentralSliceCell{exampleImgs{ee}(1)}{exampleImgs{ee}(2)}{5},...
% 		rescaleFac);
% 	
% 	
% 	OP_img(bwperim(OP_mask)) = max(OP_img(:));
% 	S5P_img(bwperim(OP_mask)) = max(S5P_img(:));
% 	S2P_img(bwperim(OP_mask)) = max(S2P_img(:));
% 
% 	
% 	Centroid = sortedS5PCentroidsCell{exampleImgs{ee}(1)}(exampleImgs{ee}(2),:);
% 	
% 	imgSize = pixelSize.*size(OP_img)./rescaleFac;
% 	
% 	subplot(numExamples,4,(ee-1).*4+1)
% 	imagesc([0,imgSize(2)],[0,imgSize(1)],OP_img)
% 	axis equal tight
% 	colormap(gray)
% 	if ee == 1
% 		hold on
% 		plot(Centroid(1)-[1,0],Centroid(2)-[1,1],'w-','LineWidth',2)
% 		title('Oligopaint','FontWeight','normal')
% 	end
% 	set(gca,'XLim',Centroid(1)+[-0.5,+0.5].*windowSize,...
% 		'YLim',Centroid(2)+[-0.5,+0.5].*windowSize,...
% 		'XTick',[],'YTick',[])
% 	ylabel(sortedCondNames{exampleImgs{ee}(1)})
% 	
% 	subplot(numExamples,4,(ee-1).*4+2)
% 	imagesc([0,imgSize(2)],[0,imgSize(1)],S5P_img)
% 	axis equal tight
% 	colormap(gray)
% 	if ee == 1
% 		title('Pol II Ser5P','FontWeight','normal')
% 	end
% 	hold on
% 	set(gca,'XLim',Centroid(1)+[-0.5,+0.5].*windowSize,...
% 		'YLim',Centroid(2)+[-0.5,+0.5].*windowSize,...
% 		'XTick',[],'YTick',[])
% 
% 
% 	subplot(numExamples,4,(ee-1).*4+3)
% 	imagesc([0,imgSize(2)],[0,imgSize(1)],S2P_img)
% 	axis equal tight
% 	colormap(gray)
% 	if ee == 1
% 		title('Pol II Ser2P','FontWeight','normal')
% 	end
% 	set(gca,'XLim',Centroid(1)+[-0.5,+0.5].*windowSize,...
% 		'YLim',Centroid(2)+[-0.5,+0.5].*windowSize,...
% 		'XTick',[],'YTick',[])
% 
% 	subplot(numExamples,4,(ee-1).*4+4)
% 	greenChannel = S2P_img-min(S2P_img(:));
% 	greenChannel = greenChannel./max(greenChannel(:));
% 	magentaChannel = S5P_img-min(S5P_img(:));
% 	magentaChannel = magentaChannel./max(magentaChannel(:));
% 	rgb_img = zeros(size(S2P_img,1),size(S2P_img,2),3);
% 	rgb_img(:,:,1) = magentaChannel;
% 	rgb_img(:,:,2) = greenChannel;
% 	rgb_img(:,:,3) = magentaChannel;
% 	image([0,imgSize(2)],[0,imgSize(1)],rgb_img)
% 	axis equal tight
% 	if ee == 1
% 		title('Merge','FontWeight','normal')
% 	end
% 	set(gca,'XLim',Centroid(1)+[-0.5,+0.5].*windowSize,...
% 		'YLim',Centroid(2)+[-0.5,+0.5].*windowSize,...
% 		'XTick',[],'YTick',[])
% 
% 	
% end