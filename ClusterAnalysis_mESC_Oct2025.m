%% Extract RNA polymerase II cluster features
%
% This script analyses the experiment image data for RNA polymerase II (Pol II)
% clusters and calculates some basic features for each cluster. Two types of Pol
% II clusters are extracted:
% - Pol II S5P clusters (Serine 5-phosphorylated Pol II, "recruited Pol II") and
% - Pol II S2P clusters (Serine 2-phosphorylated Pol II, "elongating Pol II").
% The Pol II types are labeled in separate image channels, a third channel
% contains nucleus information. The extracted clusters of both types are limited
% to the nuclei regions.
%
% The segmentation of nuclei and Pol II clusters can/must be finely tuned by the
% available control parameters in the process parameters section.
%
% Notably, the script provides an easy way to split the processing into batches
% (for example for separate processing in compute cluster jobs). Batch
% processing is controlled by two parameters, that specify the total batch count
% and the current batch index. If batch processing is used, the batch count
% parameter must remain unchanged for each single batch (because it is used to
% implicitly calculate the files to be processed in each batch). The default
% setting is a value of 1 for both parameters (indicating that all files are
% processed in a single batch):
%    batch_count = 1;
%    batch_index = 1;
% Care must be taken that each batch process writes its results in a separate
% result file; a script to recombine the batch results is provided separately.

clear all

%% Process parameter section

% switch figure plots on/off (central xy section, pauses script)
plotFlag = false;

% source directory containing the extraced single-stack MATLAB files
sourceDirectory = fullfile('.', 'ExtractedStacks', '**');

% Channels for segmentation
NucSegChannel = 1; % Channel used to detect nuclei
S5P_SegChannel = 4; % Channel used to detect Pol II S5P clusters
S2P_SegChannel = 3; % Channel used to detect Pol II S2P clusters

% image channels to store (for debugging or later visualization)
storeImgChannels = [];
% radius of area around clusters to store
ImgSquareExtension = 0; % pixels for cut-out image extension (0 for no images)

% Target channels for intensity quantification, applied for all objects
quantChannels = [2,3,4]; % Hoechst, A488 (K27ac), A594 (Ser2P), A647 (Ser5P)
quantBlurSigma = [0,0,0.15];

nuc_segBlurSigma_nucleus = 1.0; % in microns
nuc_segBlurSigma_BG_removal = 10; % in microns
nuc_segErosion = 0.5; % range of erosion (in microns) to avoid margin effects
% Use topological operation to fill holes in the nuclei segmentation masks?
% Default: 3D hole-filing, set flag to value 1
% 2D hole-filling, usefil if the stack cuts nuclei on top or bottom, so
% that 3D hole-filling does not work, set flag value to 2
% To deactivate hole-filling, set flag to any other number
fillHolesFlag = 1;

% Minimum volume of nuclei, typical ranges for a full nucieus 10-100 cubic
% microns of volume, so set a cut-off oof 10 or 30 or so
Nuc_min_vol = 10; % cubic microns
Nuc_min_sol = 0.8; % to ensure round nuclei
Nuc_min_CoV = 0.0; % to ensure transcriptionally active foci

% Inner and outer extension of a nuclear masks to cover cytoplasm
cytoMask_extension = 1.5; % in microns
cytoMask_distance = 1.0; % in microns

S5P_segBlurSigma_object = 0.06; % in microns
S5P_segBlurSigma_BG_removal = 0.1; % in microns
S5P_seg_numStdDev = 2.0;

% Cluster connection range:
S5P_DBSCAN_epsilon = 0.5; % in microns, choose 0 for no clustering

S2P_segBlurSigma_object = 0.03; % in microns
S2P_segBlurSigma_BG_removal = 0.1; % in microns
S2P_seg_numStdDev = 2.25; % number of standard deviations in robust threshold

% Minimum volumes for objects inside the nuclei
S5P_minVol = 0.005; % cubic microns
S2P_minVol = 0.005; % cubic microns

% total count of batches (use 1 to process entire data in one run)
batch_count = 1;
% current batch index
batch_index = 1;

% extracted cluster data output file
save_file = fullfile(".", "NuclearClusterFeatures", "AfterObjectAnalysis.mat");

% end of analysis parameter section, do not change anything else, all necessary
% parameters are listed above

%% Main script section

numStoreChannels = numel(storeImgChannels);
numQuantChannels = numel(quantChannels);

listing = rdir(fullfile(sourceDirectory,'*Image*.mat'));
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

clear thisFilePath thisCondInd thisCondName

% Flag to log files that were successfully analyzed and contained nuclei
validFileFlag = false(1,numFiles);

% --- analyze image stacks one by one

% Variables to store properties of nuclei
numNuclei_vec = zeros(1,numFiles);
nuc_volCell = cell(1,numFiles);
nuc_intCell = cell(1,numFiles);
cyto_intCell = cell(1,numFiles);
nuc_stdCell = cell(1,numFiles);
nuc_medianVolCell = cell(1,numFiles);
perNuc_countCell = cell(1,numFiles);
perNuc_volCell = cell(1,numFiles);

% Variable to store the pixel sizes
S5P_xyVoxelSizeCell = cell(1,numFiles);
S5P_zVoxelSizeCell = cell(1,numFiles);
S2P_xyVoxelSizeCell = cell(1,numFiles);
S2P_zVoxelSizeCell = cell(1,numFiles);

% Variables to store properties of S5P objects inside nuclei
S5P_volCell = cell(1,numFiles);
S5P_solCell = cell(1,numFiles);
S5P_eloCell = cell(1,numFiles);
S5P_intCell = cell(1,numFiles);
S5P_centCell = cell(1,numFiles);
S5P_imgCell = cell(1,numFiles);
S5P_nucIntCell = cell(1,numFiles);
S5P_nucVolCell = cell(1,numFiles);
S5P_nucClustVolCell = cell(1,numFiles);

% Variables to store properties of S2P objects inside nuclei
S2P_volCell = cell(1,numFiles);
S2P_solCell = cell(1,numFiles);
S2P_eloCell = cell(1,numFiles);
S2P_intCell = cell(1,numFiles);
S2P_centCell = cell(1,numFiles);
S2P_imgCell = cell(1,numFiles);
S2P_nucIntCell = cell(1,numFiles);
S2P_nucVolCell = cell(1,numFiles);
S2P_nucClustVolCell = cell(1,numFiles);

% File index batches
batch_size = numFiles / batch_count;
batch_ff1 = round((batch_index - 1) * batch_size) + 1;
batch_ff2 = round((batch_index) * batch_size);


t1 = tic;

% parfor ff = batch_ff1:batch_ff2
for ff = batch_ff1:batch_ff2

	t2 = tic;

	fprintf('Processing file %d of %d\n',ff,numFiles);

	% thisCondInd = condInds(ff);
	thisFilePath = listing(ff).name;

	loadStruct = load(thisFilePath,...
		'imgStack','imgSize','pixelSize','zStepSize');
	imgStack = loadStruct.imgStack;
	% imgSize = loadStruct.imgSize;
	pixelSize = loadStruct.pixelSize;
	zStepSize = loadStruct.zStepSize;

	% Nuclei segmentation

    % Note: Several functions (e.g. histogram calculation with the old hist
    % function, but also the calculation of the standard deviation, as used on
    % the voxel values of the segmented nucleus candidates) don't work on
    % integer type data. The image data are therefore cast to floating-point
    % format. In order to retain the interpretability of intensity threshold
    % values, a basic cast is used instead of the im2double function.

	segImg = double(imgStack{NucSegChannel});
    if nuc_segBlurSigma_nucleus>0
        segImg = ...
            + imgaussfilt(segImg,nuc_segBlurSigma_nucleus./pixelSize) ...
            - imgaussfilt(segImg,nuc_segBlurSigma_BG_removal./pixelSize);
    else
        segImg = ...
            + segImg ...
            - imgaussfilt(segImg,nuc_segBlurSigma_BG_removal./pixelSize);
    end

    % Note: The background removal yields negative values when used on formats
    % that support it. This significantly influences the binarization
    % thresholding later on. If results similar to unsigned format background
    % removal are preferrable, negative values can be clipped to 0.

    % segImg = clip(segImg, 0, Inf);

    % Note: It has been recommended for some time to replace the hist function
    % with the newer histogram (for graphical output) or histcounts (for pure
    % calculations) functions. The results are not numerically identical, but
    % very similar.

	% [bin_counts,bin_centers] = hist(segImg(:),1000);
    [bin_counts,bin_edges] = histcounts(segImg(:),1000);
    bin_centers = mean([bin_edges(1:end-1);bin_edges(2:end)]);
	[nuc_seg_thresh,~] = otsuLimit(bin_centers,bin_counts,[0,Inf]);
	NucSegMask = segImg>1.0.*nuc_seg_thresh;
    if fillHolesFlag == 1
        % 3D hole-filling, default
        NucSegMask = imfill(NucSegMask,18,'holes');
    elseif fillHolesFlag == 2
        % 2D hole-filling, useful if the stack cuts most nuclei on top or
        % bottom
        NucSegMask = imfill(NucSegMask,8,'holes');
    end
    if nuc_segErosion>0
        se = strel('disk', ...
            round(nuc_segErosion./pixelSize)); % Two-dimensional erosion disk
        NucSegMask = imerode(NucSegMask,se);
    end

    if plotFlag

        subplot(1,3,1) %#ok<UNRCH>
	    imagesc(squeeze(imgStack{NucSegChannel}(:,:,ceil(size(segImg,3)./2))))
	    axis tight equal

	    subplot(1,3,2)
	    imagesc(squeeze(segImg(:,:,ceil(size(segImg,3)./2))))
	    axis tight equal

	    subplot(1,3,3)
	    imagesc(squeeze(NucSegMask(:,:,ceil(size(segImg,3)./2))))
	    axis tight equal

        % Uncomment the following two lines, and make sure you're in serial loop
        % processing mode, if you want to check the extracted images one by one

        % fprintf('File name: %s\n',thisFilePath)
        % waitforbuttonpress

    end % plotFlag

	% --- Connected component segmentation of nuclei
	comps = bwconncomp(NucSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	minPixels = Nuc_min_vol./(pixelSize.^2)./zStepSize;
	comps.NumObjects = sum(numPxls>=minPixels);
	comps.PixelIdxList = comps.PixelIdxList(numPxls>=minPixels);
	% numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);

	props = regionprops3(comps,double(imgStack{NucSegChannel}),...
		'Solidity','VoxelValues');

    if comps.NumObjects>0
        % Note: The current implementation of filtering only considers the
        % solidity feature. However, it seems that at some point, it was
        % supposed to consider the coefficient of variance as well. We should
        % discuss this.
        Solidity_array = [props.Solidity];
        % CoV_array = ...
        %     cellfun(@(vals)std(vals(:))./mean(vals(:)),...
        %     props.VoxelValues);
        % inclNucInds = Solidity_array>=Nuc_min_sol ...
        %     & CoV_array>=Nuc_min_CoV;
        comps.NumObjects = sum(Solidity_array>=Nuc_min_sol);
        comps.PixelIdxList = comps.PixelIdxList(Solidity_array>=Nuc_min_sol);
        % numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
    end

	numNuclei = comps.NumObjects;
	numNuclei_vec(ff) = numNuclei;

    if comps.NumObjects>0

        validFileFlag(ff) = true;

        nuc_intCell{ff} = cell(1,numQuantChannels);
        cyto_intCell{ff} = cell(1,numQuantChannels);
        nuc_stdCell{ff} = cell(1,numQuantChannels);
        for qq = 1:numQuantChannels
            quantImg = double(imgStack{quantChannels(qq)});
            quantProps = regionprops3(comps,quantImg,...
                'MeanIntensity','VoxelIdxList','VoxelValues');
            nuc_intCell{ff}{qq} = [quantProps.MeanIntensity];
            cyto_intCell{ff}{qq} = zeros(numNuclei,1);
            nuc_stdCell{ff}{qq} = cellfun(@std,quantProps.VoxelValues);

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

        props = regionprops3(comps,double(imgStack{NucSegChannel}),...
            'Volume','VoxelValues','Solidity','VoxelIdxList',...
            'BoundingBox');

        Volume_array = [props.Volume].*pixelSize.^2.*zStepSize;
        % Intensity_array = cellfun(@(vals)median(vals),props.VoxelValues);
        % Solidity_array = [props.Solidity];

        nuc_volCell{ff} = Volume_array;


        % --- For each nucleus, get objects from the different channels

        nuc_medianVolCell{ff} = cell(1,2);
        perNuc_countCell{ff} = cell(1,2);
        perNuc_volCell{ff} = cell(1,2);
        for qq = 1:2
            nuc_medianVolCell{ff}{qq} = zeros(numNuclei,1);
            perNuc_countCell{ff}{qq} = zeros(numNuclei,1);
            perNuc_volCell{ff}{qq} = cell(numNuclei,1);
        end

        S5P_xyVoxelSize = cell(1,numNuclei);
        S5P_zVoxelSize = cell(1,numNuclei);

        S2P_xyVoxelSize = cell(1,numNuclei);
        S2P_zVoxelSize = cell(1,numNuclei);

        S5P_volume = cell(1,numNuclei);
        S5P_solidity = cell(1,numNuclei);
        S5P_elongation = cell(1,numNuclei);
        S5P_centralSlices_store = cell(1,numNuclei);
        S5P_intensity = cell(numQuantChannels,numNuclei);
        S5P_cent_store = cell(1,numNuclei);
        S5P_nucIntensity = cell(numQuantChannels,numNuclei);
        S5P_nucVolume = cell(1,numNuclei);
        S5P_nucClustVolume = cell(1,numNuclei);

        S2P_volume = cell(1,numNuclei);
        S2P_solidity = cell(1,numNuclei);
        S2P_elongation = cell(1,numNuclei);
        S2P_centralSlices_store = cell(1,numNuclei);
        S2P_intensity = cell(numQuantChannels,numNuclei);
        S2P_cent_store = cell(1,numNuclei);
        S2P_nucIntensity = cell(numQuantChannels,numNuclei);
        S2P_nucVolume = cell(1,numNuclei);
        S2P_nucClustVolume = cell(1,numNuclei);

        for nn = 1:numNuclei

            boxArray = props.BoundingBox(nn,:);

            % Note: As with the nuclei channel above, several functions (e.g.
            % the calculation of the standard deviation) don't work on integer
            % type data. The image data are therefore cast to floating-point
            % format. In order to retain the interpretability of intensity
            % threshold values, a basic cast is used instead of the im2double
            % function.

            S5P_subImage = double(imgStack{S5P_SegChannel}(...
                boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
                boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
                boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5));

            S2P_subImage = double(imgStack{S2P_SegChannel}(...
                boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
                boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
                boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5));

            S5P_subImage = ...
                + imgaussfilt(S5P_subImage,S5P_segBlurSigma_object./pixelSize) ...
                - imgaussfilt(S5P_subImage,S5P_segBlurSigma_BG_removal./pixelSize);
            S2P_subImage = ...
                + imgaussfilt(S2P_subImage,S2P_segBlurSigma_object./pixelSize) ...
                - imgaussfilt(S2P_subImage,S2P_segBlurSigma_BG_removal./pixelSize);

            % Note: As with the nuclei channel above, the background removal
            % yields negative values when used on formats that support it. This
            % significantly influences the binarization thresholding. If results
            % similar to unsigned format background removal are preferrable,
            % negative values can be clipped to 0.

            % S5P_subImage = clip(S5P_subImage, 0, Inf);
            % S2P_subImage = clip(S2P_subImage, 0, Inf);

            % Note: imgSize was used here to create the mask image, but the
            % sizes in imgSize are in XYZ order
            % NucMask = false(imgSize);
            NucMask = false(size(NucSegMask));
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

            % For storage of example images
            if numStoreChannels > 0
                store_subImages = cell(1,numStoreChannels);
                for color = 1:numStoreChannels
                    store_subImages{color} = imgStack{storeImgChannels(color)}(...
                        boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
                        boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
                        boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);
                    nuc_intensities = ...
                        store_subImages{color}(NucMask_subImage);
                    store_subImages{color} = ...
                        store_subImages{color}./median(nuc_intensities(:));
                end
            end

            S5P_comps = bwconncomp(S5P_mask,18);
            S5P_numPxls = cellfun(@(elmt)numel(elmt),S5P_comps.PixelIdxList);
            minPixels = S5P_minVol./(pixelSize.^2)./zStepSize;
            S5P_comps.NumObjects = sum(S5P_numPxls>minPixels);
            S5P_comps.PixelIdxList = S5P_comps.PixelIdxList(S5P_numPxls>minPixels);
            % S5P_numPxls = cellfun(@(elmt)numel(elmt),S5P_comps.PixelIdxList);

            S2P_comps = bwconncomp(S2P_mask,18);
            S2P_numPxls = cellfun(@(elmt)numel(elmt),S2P_comps.PixelIdxList);
            minPixels = S2P_minVol./(pixelSize.^2)./zStepSize;
            S2P_comps.NumObjects = sum(S2P_numPxls>minPixels);
            S2P_comps.PixelIdxList = S2P_comps.PixelIdxList(S2P_numPxls>minPixels);
            % S2P_numPxls = cellfun(@(elmt)numel(elmt),S2P_comps.PixelIdxList);

            if S5P_comps.NumObjects>0 && S2P_comps.NumObjects>0
                % This if-clause still needs to be split into the two
                % separate object classes. There is no reason to not
                % analyze one object class because the other one is not
                % detectable

                % DBSCAN clustering of labeled regions
                if S5P_DBSCAN_epsilon > 0
                    S5P_props = regionprops3(S5P_comps,S5P_subImage,...
                        'Centroid');
                    centroid_coords = ...
                        S5P_props.Centroid.*[pixelSize,pixelSize,zStepSize];
                    dbscan_inds = ...
                        dbscan(centroid_coords,S5P_DBSCAN_epsilon,1);

                    unique_inds = unique(dbscan_inds);
                    num_inds = numel(unique_inds);
                    updated_comps = S5P_comps;
                    updated_comps.NumObjects = num_inds;
                    updated_comps.PixelIdxList = cell(1,num_inds);
                    for ii = 1:num_inds
                        updated_comps.PixelIdxList{ii} = ...
                            sort(vertcat(S5P_comps.PixelIdxList{...
                            dbscan_inds==unique_inds(ii)} ...
                            ));
                    end
                    S5P_comps = updated_comps;
                end

                LL = labelmatrix(S5P_comps);

                if plotFlag

                    subplot(2,2,1) %#ok<UNRCH>
                    centerPlaneInd = round(boxArray(6).*0.5);
                    % 			imagesc(squeeze(max(S5P_subImage,[],3)))
                    imagesc(squeeze(S5P_subImage(:,:,centerPlaneInd)))
                    axis tight equal
                    set(gca,'Colormap',gray)

                    subplot(2,2,2)
                    % 			imagesc(squeeze(max(S2P_subImage,[],3)))
                    imagesc(squeeze(S2P_subImage(:,:,centerPlaneInd)))
                    axis tight equal
                    set(gca,'Colormap',gray)

                    subplot(2,2,3)
                    % 			imagesc(squeeze(max(S5P_mask,[],3)))
                    imagesc(squeeze(LL(:,:,centerPlaneInd)))
                    axis tight equal
                    set(gca,'Colormap',lines)

                    subplot(2,2,4)
                    % 			imagesc(squeeze(max(S2P_mask,[],3)))
                    imagesc(squeeze(S2P_mask(:,:,centerPlaneInd)))
                    axis tight equal
                    set(gca,'Colormap',gray)

                    % Uncomment the following line, and make sure you're in
                    % serial loop processing mode, if you want to see the
                    % segmentation results for the two types of foci

                    % waitforbuttonpress

                end % plotFlag

                S5P_props = regionprops3(S5P_comps,S5P_subImage,...
                    'Volume','Solidity',...
                    'Centroid','Image','BoundingBox');

                S5P_Volume_array = ...
                    [S5P_props.Volume].*pixelSize.^2.*zStepSize;
                S5P_Solidity_array = [S5P_props.Solidity];
                S5P_Centroid_array = ...
                    S5P_props.Centroid.*[pixelSize,pixelSize,zStepSize];

                S2P_props = regionprops3(S2P_comps,S2P_subImage,...
                    'Volume','Solidity',...
                    'Centroid','Image','BoundingBox');

                S2P_Volume_array = ...
                    [S2P_props.Volume].*pixelSize.^2.*zStepSize;
                S2P_Solidity_array = [S2P_props.Solidity];
                S2P_Centroid_array = ...
                    S2P_props.Centroid.*[pixelSize,pixelSize,zStepSize];

                perNuc_countCell{ff}{1}(nn) = numel(S5P_Volume_array);
                perNuc_countCell{ff}{2}(nn) = numel(S2P_Volume_array);

                perNuc_volCell{ff}{1}{nn} = S5P_Volume_array;
                perNuc_volCell{ff}{2}{nn} = S2P_Volume_array;

                nuc_medianVolCell{ff}{1}(nn) = median(S5P_Volume_array);
                nuc_medianVolCell{ff}{2}(nn) = median(S2P_Volume_array);

                S5P_xyVoxelSize{nn} = ...
                    pixelSize.*ones(size(S5P_Volume_array));
                S5P_zVoxelSize{nn} = ...
                    zStepSize.*ones(size(S5P_Volume_array));

                S2P_xyVoxelSize{nn} = ...
                    pixelSize.*ones(size(S2P_Volume_array));
                S2P_zVoxelSize{nn} = ...
                    zStepSize.*ones(size(S2P_Volume_array));


                % --- get cluster central plane and elongation in-plane
                S5P_Elongation_array = ...
                    zeros(size(S5P_Solidity_array));
                S5P_Slices_cell = ...
                    cell(size(S5P_Solidity_array));

                for object_nn = 1:numel(S5P_Volume_array)

                    % Note: This bounding box had apparently been used at some
                    % point in the upcoming loop to determine S5P_Slices_cell.
                    % The current implementation uses the image limits instead,
                    % which are calculated hereafter. As long as the bounding
                    % box is no longer needed, I'll comment these calculations
                    % out.

                    % boundingBox = S5P_props.BoundingBox(object_nn,:);
                    % thisImage = squeeze(S5P_subImage(...
                    %     boundingBox(2)+0.5:boundingBox(2)+boundingBox(5)-0.5,...
                    %     boundingBox(1)+0.5:boundingBox(1)+boundingBox(4)-0.5,...
                    %     ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
                    % thisImage = thisImage-min(thisImage(:));
                    % thisImage = thisImage./max(thisImage(:));
                    
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
                    % thisImage((bwperim(thisMask))) = 0;

                    S5P_Slices_cell{object_nn} = cell(1,numStoreChannels+2);
                    centroid_1 = round(S5P_props.Centroid(object_nn,2));
                    centroid_2 = round(S5P_props.Centroid(object_nn,1));
                    center_z = round(S5P_props.Centroid(object_nn,3));
                    img_limits = [...
                        centroid_1-ImgSquareExtension,...
                        centroid_1+ImgSquareExtension,...
                        centroid_2-ImgSquareExtension,...
                        centroid_2+ImgSquareExtension];
                    img_limits = [...
                        max(1,img_limits(1)),...
                        min(subImgSize(1),img_limits(2)),...
                        max(1,img_limits(3)),...
                        min(subImgSize(2),img_limits(4))];

                    for color = 1:numStoreChannels
                        % S5P_Slices_cell{object_nn}{color} = ...
                        %     squeeze(store_subImages{color}(:,:,...
                        %     ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));

                        S5P_Slices_cell{object_nn}{color} = ...
                            squeeze(store_subImages{color}(...
                            img_limits(1):img_limits(2),...
                            img_limits(3):img_limits(4),...
                            center_z));
                    end

                    S5P_Slices_cell{object_nn}{numStoreChannels+1} = ...
                        squeeze(LL(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));
                    S5P_Slices_cell{object_nn}{numStoreChannels+2} = ...
                        squeeze(S2P_mask(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));

                    if sum(sum(uint8(thisMask)))>0
                        thisProps = regionprops(uint8(thisMask),...
                            'MajorAxisLength','MinorAxisLength');
                        S5P_Elongation_array(object_nn) = ...
                            thisProps.MajorAxisLength./thisProps.MinorAxisLength;
                    else
                        S5P_Elongation_array(object_nn) = NaN;
                    end
                end
                S5P_volume{nn} = S5P_Volume_array;
                S5P_solidity{nn} = S5P_Solidity_array;
                S5P_elongation{nn} = S5P_Elongation_array;
                S5P_cent_store{nn} = S5P_Centroid_array;
                S5P_centralSlices_store{nn} = S5P_Slices_cell;

                S5P_nucVolume{nn} = ...
                        ones(size(S5P_Volume_array)) ...
                        .*Volume_array(nn);

                % Note: This looks a bit weird; the same double vector
                % S5P_Volume_array of length N is reproduced N times to be
                % stored in a cell array of length N. It is not touched later on
                % until it is stored in the result structure. I guess what is
                % intended is to simply store the double vector once, as in:
                %   S5P_nucClustVolume{nn} = S5P_Volume_array;
                % I'll do this change, but we have to discuss this later
                S5P_nucClustVolume{nn} = S5P_Volume_array;
                % S5P_nucClustVolume{nn} = ...
                %     cell(size(S5P_Volume_array));
                % for object_nn = 1:numel(S5P_Volume_array)
                %     S5P_nucClustVolume{nn}{object_nn} = ...
                %         S5P_Volume_array;
                % end

                S2P_Elongation_array = ...
                    zeros(size(S2P_Solidity_array));
                S2P_Slices_cell = ...
                    cell(size(S2P_Solidity_array));

                for object_nn = 1:numel(S2P_Volume_array)

                    % boundingBox = S2P_props.BoundingBox(object_nn,:);
                    % thisImage = squeeze(S2P_subImage(...
                    %     boundingBox(2)+0.5:boundingBox(2)+boundingBox(5)-0.5,...
                    %     boundingBox(1)+0.5:boundingBox(1)+boundingBox(4)-0.5,...
                    %     ceil(boundingBox(3)+0.5+0.5.*(boundingBox(6)-1))));
                    % thisImage = thisImage-min(thisImage(:));
                    % thisImage = thisImage./max(thisImage(:));

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
                    % thisImage((bwperim(thisMask))) = 0;

                    centroid_1 = round(S2P_props.Centroid(object_nn,2));
                    centroid_2 = round(S2P_props.Centroid(object_nn,1));
                    center_z = round(S2P_props.Centroid(object_nn,3));
                    img_limits = [...
                        centroid_1-ImgSquareExtension,...
                        centroid_1+ImgSquareExtension,...
                        centroid_2-ImgSquareExtension,...
                        centroid_2+ImgSquareExtension];
                    img_limits = [...
                        max(1,img_limits(1)),...
                        min(subImgSize(1),img_limits(2)),...
                        max(1,img_limits(3)),...
                        min(subImgSize(2),img_limits(4))];

                    for color = 1:numStoreChannels
                        S2P_Slices_cell{object_nn}{color} = ...
                            squeeze(store_subImages{color}(...
                            img_limits(1):img_limits(2),...
                            img_limits(3):img_limits(4),...
                            center_z));
                    end
                    S2P_Slices_cell{object_nn}{numStoreChannels+1} = ...
                        squeeze(S5P_mask(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));
                    S2P_Slices_cell{object_nn}{numStoreChannels+2} = ...
                        squeeze(S2P_mask(...
                        img_limits(1):img_limits(2),...
                        img_limits(3):img_limits(4),...
                        center_z));

                    if sum(sum(uint8(thisMask)))>0
                        thisProps = regionprops(uint8(thisMask),...
                            'MajorAxisLength','MinorAxisLength');
                        S2P_Elongation_array(object_nn) = ...
                            thisProps.MajorAxisLength./thisProps.MinorAxisLength;
                    else
                        S2P_Elongation_array(object_nn) = NaN;
                    end

                end
                S2P_volume{nn} = S2P_Volume_array;
                S2P_solidity{nn} = S2P_Solidity_array;
                S2P_elongation{nn} = S2P_Elongation_array;
                S2P_cent_store{nn} = S2P_Centroid_array;
                S2P_centralSlices_store{nn} = S2P_Slices_cell;

                S2P_nucVolume{nn} = ...
                        ones(size(S2P_Volume_array)) ...
                        .*Volume_array(nn);

                % Note: This looks a bit weird; the same double vector
                % S5P_Volume_array of length N is reproduced N times to be
                % stored in a cell array of length N. It is not touched later on
                % until it is stored in the result structure. I guess what is
                % intended is to simply store the double vector once, as in:
                %   S2P_nucClustVolume{nn} = S2P_Volume_array;
                % I'll do this change, but we have to discuss this later
                S2P_nucClustVolume{nn} = S2P_Volume_array;
                % S2P_nucClustVolume{nn} = ...
                %     cell(size(S2P_Volume_array));
                % for object_nn = 1:numel(S2P_Volume_array)
                %     S2P_nucClustVolume{nn}{object_nn} = ...
                %         S2P_Volume_array;
                % end

                % --- quantification for all target channels
                for qq = 1:numQuantChannels

                    channelInd = quantChannels(qq);
                    quant_subImage = double(imgStack{channelInd}(...
                        boxArray(2)+0.5:boxArray(2)+boxArray(5)-0.5,...
                        boxArray(1)+0.5:boxArray(1)+boxArray(4)-0.5,...
                        boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5));

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
                    S5P_nucIntensity{qq,nn} = ...
                        ones(size(S5P_intensity{qq,nn})) ...
                        .*nuc_intCell{ff}{qq}(nn);

                    S2P_quant_props = regionprops3(...
                        S2P_comps,quant_subImage,'VoxelValues');
                    Quant_ClusterMedian = cellfun(...
                        @(vals)median(vals),S2P_quant_props.VoxelValues);
                    S2P_intensity{qq,nn} = ...
                        Quant_ClusterMedian./Quant_nucleusMedian;
                    S2P_nucIntensity{qq,nn} = ...
                        ones(size(S2P_intensity{qq,nn})) ...
                        .*nuc_intCell{ff}{qq}(nn);

                end

            else

                S5P_xyVoxelSize{nn} = [];
                S5P_zVoxelSize{nn} = [];
                S2P_xyVoxelSize{nn} = [];
                S2P_zVoxelSize{nn} = [];

                S5P_volume{nn} = [];
                S5P_solidity{nn} = [];
                S5P_elongation{nn} = [];
                S5P_centralSlices_store{nn} = {};
                S5P_cent_store{nn} = [];
                S5P_nucVolume{nn} = [];
                S5P_nucClustVolume{nn} = [];

                S2P_volume{nn} = [];
                S2P_solidity{nn} = [];
                S2P_elongation{nn} = [];
                S2P_centralSlices_store{nn} = {};
                S2P_cent_store{nn} = [];
                S2P_nucVolume{nn} = [];
                S2P_nucClustVolume{nn} = [];

                for qq = 1:numQuantChannels
                    S5P_intensity{qq,nn} = [];
                    S5P_nucIntensity{qq,nn} = [];
                    S2P_intensity{qq,nn} = [];
                    S2P_nucIntensity{qq,nn} = [];
                end

            end

        end

        S5P_xyVoxelSizeCell{ff} = vertcat(S5P_xyVoxelSize{:});
        S5P_zVoxelSizeCell{ff} = vertcat(S5P_zVoxelSize{:});
        S2P_xyVoxelSizeCell{ff} = vertcat(S2P_xyVoxelSize{:});
        S2P_zVoxelSizeCell{ff} = vertcat(S2P_zVoxelSize{:});

        S5P_volCell{ff} = vertcat(S5P_volume{:});
        S5P_solCell{ff} = vertcat(S5P_solidity{:});
        S5P_eloCell{ff} = vertcat(S5P_elongation{:});
        S5P_imgCell{ff} = vertcat(S5P_centralSlices_store{:});
        S5P_centCell{ff} = vertcat(S5P_cent_store{:});
        S5P_intCell{ff} = cell(1,numQuantChannels);
        S5P_nucIntCell{ff} = cell(1,numQuantChannels);
        for qq = 1:numQuantChannels
            S5P_intCell{ff}{qq} = vertcat(S5P_intensity{qq,:});
            S5P_nucIntCell{ff}{qq} = vertcat(S5P_nucIntensity{qq,:});
        end
        % Note: I strongly assume that these two variables should also be
        % ocmbined to a single large double vector, just as all the others
        % before. I'll do this, but we have to discuss this later.
        S5P_nucVolCell{ff} = vertcat(S5P_nucVolume{:});
        S5P_nucClustVolCell{ff} = vertcat(S5P_nucClustVolume{:});

        S2P_volCell{ff} = vertcat(S2P_volume{:});
        S2P_solCell{ff} = vertcat(S2P_solidity{:});
        S2P_eloCell{ff} = vertcat(S2P_elongation{:});
        S2P_imgCell{ff} = vertcat(S2P_centralSlices_store{:});
        S2P_centCell{ff} = vertcat(S2P_cent_store{:});
        S2P_intCell{ff} = cell(1,numQuantChannels);
        S2P_nucIntCell{ff} = cell(1,numQuantChannels);
        for qq = 1:numQuantChannels
            S2P_intCell{ff}{qq} = vertcat(S2P_intensity{qq,:});
            S2P_nucIntCell{ff}{qq} = vertcat(S2P_nucIntensity{qq,:});
        end
        % Note: I strongly assume that these two variables should also be
        % ocmbined to a single large double vector, just as all the others
        % before. I'll do this, but we have to discuss this later.
        S2P_nucVolCell{ff} = vertcat(S2P_nucVolume{:});
        S2P_nucClustVolCell{ff} = vertcat(S2P_nucClustVolume{:});

    end

	fprintf('Processing file %d of %d finished\n',ff,numFiles);

    toc(t2)

end

toc(t1)

%% Retain only files that returned nuclei

condInds = condInds(validFileFlag);
condNames = condNames(validFileFlag);
numNuclei_vec = numNuclei_vec(validFileFlag);
nuc_volCell = nuc_volCell(validFileFlag);
nuc_intCell = nuc_intCell(validFileFlag);
cyto_intCell = cyto_intCell(validFileFlag);
nuc_stdCell = nuc_stdCell(validFileFlag);
nuc_medianVolCell = nuc_medianVolCell(validFileFlag);
perNuc_countCell = perNuc_countCell(validFileFlag);
perNuc_volCell = perNuc_volCell(validFileFlag);

S5P_xyVoxelSizeCell = S5P_xyVoxelSizeCell(validFileFlag);
S5P_zVoxelSizeCell = S5P_zVoxelSizeCell(validFileFlag);
S2P_xyVoxelSizeCell = S2P_xyVoxelSizeCell(validFileFlag);
S2P_zVoxelSizeCell = S2P_zVoxelSizeCell(validFileFlag);

S5P_volCell = S5P_volCell(validFileFlag);
S5P_solCell = S5P_solCell(validFileFlag);
S5P_eloCell = S5P_eloCell(validFileFlag);
S5P_imgCell = S5P_imgCell(validFileFlag);
S5P_centCell = S5P_centCell(validFileFlag);
S5P_intCell = S5P_intCell(validFileFlag);
S5P_nucIntCell = S5P_nucIntCell(validFileFlag);
S5P_nucClustVolCell = S5P_nucClustVolCell(validFileFlag);
S5P_nucVolCell = S5P_nucVolCell(validFileFlag);

S2P_volCell = S2P_volCell(validFileFlag);
S2P_solCell = S2P_solCell(validFileFlag);
S2P_eloCell = S2P_eloCell(validFileFlag);
S2P_imgCell = S2P_imgCell(validFileFlag);
S2P_centCell = S2P_centCell(validFileFlag);
S2P_intCell = S2P_intCell(validFileFlag);
S2P_nucIntCell = S2P_nucIntCell(validFileFlag);
S2P_nucClustVolCell = S2P_nucClustVolCell(validFileFlag);
S2P_nucVolCell = S2P_nucVolCell(validFileFlag);

%% Save results on disk to make available for future analysis

clear( ...
    "CoV_array", "Intensity_array", "LL", "NucMask", "NucMask_subImage", ...
    "NucSegMask", "Quant_ClusterMedian", "Quant_nucleusMedian", ...
    "S2P_Centroid_array", "S2P_Elongation_array", "S2P_Slices_cell", ...
    "S2P_Solidity_array", "S2P_Volume_array", "S2P_cent_store", ...
    "S2P_centralSlices_store", "S2P_comps", "S2P_elongation", "S2P_intensity", ...
    "S2P_mask", "S2P_nucClustVolume", "S2P_nucIntensity", "S2P_nucVolume", ...
    "S2P_numPxls", "S2P_props", "S2P_quant_props", "S2P_solidity", ...
    "S2P_subImage", "S2P_volume", "S2P_xyVoxelSize", "S2P_zVoxelSize", ...
    "S5P_Centroid_array", "S5P_Elongation_array", "S5P_Solidity_array", ...
    "S5P_Volume_array", "S5P_Slices_cell", "S5P_cent_store", ...
    "S5P_centralSlices_store", "S5P_comps", "S5P_elongation", "S5P_intensity", ...
    "S5P_mask", "S5P_nucClustVolume", "S5P_nucIntensity", "S5P_nucVolume", ...
    "S5P_numPxls", "S5P_props", "S5P_quant_props", "S5P_solidity", ...
    "S5P_subImage", "S5P_volume", "S5P_xyVoxelSize", "S5P_zVoxelSize", ...
    "Solidity_array", "Volume_array", "bin_centers", "bin_counts", "bin_edges", ...
    "boundingBox", "boxArray", "centerInd", "center_z", "centroid_1", ...
    "centroid_2", "centroid_coords", "channelInd", "color", "comps", "coreMask", ...
    "cytoMask", "dbscan_inds", "ff", "ii", "imgSize", "imgStack", "img_limits", ...
    "inclNucInds", "loadStruct", "minPixels", "nn", "nuc_seg_thresh", ...
    "numNuclei", "numPxls", "num_inds", "object_nn", "pixelSize", "props", "qq", ...
    "quantImg", "quantProps", "quant_subImage", "se", "segImg", ...
    "seg_intensities", "seg_mean", "seg_std", "subImgSize", "t1", "t2", ...
    "thisCondInd", "thisCondName", "thisFilePath", "thisImage", "thisMask", ...
    "thisProps", "unique_inds", "updated_comps", "zStepSize");

folder = fileparts(save_file);
if ~isfolder(folder)
    mkdir(folder);
end
clear folder;

% make sure not to use the v7.3 file format, it would become GIGANTIC because of
% the nested cell structures
save(save_file)

%% Sort into conditions

uniqueCondNames = unique(condNames);
numPlotSets = numel(uniqueCondNames);
fileIndsCell = cell(1,numPlotSets);
numFiles_perCond = zeros(1,numPlotSets);
for cc = 1:numPlotSets
	fileIndsCell{cc} = cellfun(...
		@(elmt)strcmp(elmt,uniqueCondNames{cc}),condNames);
	numFiles_perCond(cc) = sum(fileIndsCell{cc});
end

sortedCondNames = cell(1,numPlotSets);
sortedNumFiles = zeros(1,numPlotSets);

sortedNucVolCell = cell(1,numQuantChannels);
sortedNucIntCell = cell(1,numQuantChannels);
sortedCytoIntCell = cell(1,numQuantChannels);
sortedNucStdCell = cell(1,numQuantChannels);

sortedS5PPixelSize_xy = cell(1,numPlotSets);
sortedS5PPixelSize_z = cell(1,numPlotSets);
sortedS2PPixelSize_xy = cell(1,numPlotSets);
sortedS2PPixelSize_z = cell(1,numPlotSets);

sortedS5PNumCell = cell(1,numPlotSets);
sortedS5PVolCell = cell(1,numPlotSets);
sortedS5PVolPerNucCell = cell(1,numPlotSets);
sortedS5PSolCell = cell(1,numPlotSets);
sortedS5PEloCell = cell(1,numPlotSets);
sortedS5PCentralSliceCell = cell(1,numPlotSets);
sortedS5PCentroidsCell = cell(1,numPlotSets);
sortedS5PIntCell = cell(1,numQuantChannels);
sortedS5PNucIntCell = cell(1,numQuantChannels);
sortedS5PNucVolCell = cell(1,numPlotSets);  
sortedS5PNucClustVolCell = cell(1,numPlotSets);

sortedS2PNumCell = cell(1,numPlotSets);
sortedS2PVolCell = cell(1,numPlotSets);
sortedS2PVolPerNucCell = cell(1,numPlotSets);
sortedS2PSolCell = cell(1,numPlotSets);
sortedS2PEloCell = cell(1,numPlotSets);
sortedS2PCentralSliceCell = cell(1,numPlotSets);
sortedS2PCentroidsCell = cell(1,numPlotSets);
sortedS2PIntCell = cell(1,numQuantChannels);
sortedS2PNucIntCell = cell(1,numQuantChannels);
sortedS2PNucVolCell = cell(1,numPlotSets);  
sortedS2PNucClustVolCell = cell(1,numPlotSets);

for qq = 1:numQuantChannels
	sortedNucIntCell{qq} = cell(1,numPlotSets);
	sortedCytoIntCell{qq} = cell(1,numPlotSets);
    sortedNucStdCell{qq} = cell(1,numPlotSets);
	sortedS5PIntCell{qq} = cell(1,numPlotSets);
    sortedS5PNucIntCell{qq} = cell(1,numPlotSets);
	sortedS2PIntCell{qq} = cell(1,numPlotSets);
    sortedS2PNucIntCell{qq} = cell(1,numPlotSets);
end


for cc = 1:numPlotSets
	
	sortedCondNames{cc} = ...
		condNames(fileIndsCell{cc});
	sortedCondNames{cc} = sortedCondNames{cc}{1};
	
	sortedNumFiles(cc) = sum(fileIndsCell{cc});
	
    sortedNucVolCell{cc} = ...
        vertcat(nuc_volCell{fileIndsCell{cc}});

    sortedS5PPixelSize_xy{cc} = ...
        vertcat(S5P_xyVoxelSizeCell{fileIndsCell{cc}});
    sortedS5PPixelSize_z = ...
        vertcat(S5P_zVoxelSizeCell{fileIndsCell{cc}});
    sortedS2PPixelSize_xy = ...
        vertcat(S2P_xyVoxelSizeCell{fileIndsCell{cc}});
    sortedS2PPixelSize_z = ...
        vertcat(S2P_zVoxelSizeCell{fileIndsCell{cc}});

	S5P_nums = vertcat(arrayfun(...
		@(val)perNuc_countCell{val}{1},find(fileIndsCell{cc}),...
		'UniformOutput',false));
	S5P_nums = vertcat(S5P_nums{:});
	S5P_vols = vertcat(S5P_volCell{fileIndsCell{cc}});
	S5P_perNucVols = vertcat(arrayfun(...
        @(val)perNuc_volumeCell{val}{1},find(fileIndsCell{cc}),...
        'UniformOutput',false));
    S5P_perNucVols = vertcat(S5P_perNucVols{:});
    S5P_sols = vertcat(S5P_solCell{fileIndsCell{cc}});
	S5P_elos = vertcat(S5P_eloCell{fileIndsCell{cc}});
	S5P_slices = vertcat(S5P_imgCell{fileIndsCell{cc}});
	S5P_centroids = vertcat(S5P_centCell{fileIndsCell{cc}});
	S5P_ints = S5P_intCell(fileIndsCell{cc});
    S5P_nucInts = S5P_nucIntCell(fileIndsCell{cc});
    S5P_nucVols = horzcat(S5P_nucVolCell{fileIndsCell{cc}})';
    S5P_nucClustVols = S5P_nucClustVolCell(fileIndsCell{cc});
    
	sortedS5PNumCell{cc} = S5P_nums;
    sortedS5PVolCell{cc} = S5P_vols;
    sortedS5PVolPerNucCell{cc} = S5P_perNucVols;
    sortedS5PSolCell{cc} = S5P_sols;
	sortedS5PEloCell{cc} = S5P_elos;
	sortedS5PCentralSliceCell{cc} = S5P_slices;
	sortedS5PCentroidsCell{cc} = S5P_centroids;
    sortedS5PNucVolCell{cc} = vertcat(S5P_nucVols{:});
    sortedS5PNucClustVolCell{cc} = horzcat(S5P_nucClustVols{:})';
	
                    
	S2P_nums = vertcat(arrayfun(...
		@(val)perNuc_countCell{val}{2},find(fileIndsCell{cc}),...
		'UniformOutput',false));
	S2P_nums = vertcat(S2P_nums{:});
	S2P_vols = vertcat(S2P_volCell{fileIndsCell{cc}});
	S2P_perNucVols = vertcat(arrayfun(...
        @(val)perNuc_volumeCell{val}{2},find(fileIndsCell{cc}),...
        'UniformOutput',false));
    S2P_sols = vertcat(S2P_solCell{fileIndsCell{cc}});
	S2P_elos = vertcat(S2P_eloCell{fileIndsCell{cc}});
	S2P_slices = vertcat(S2P_imgCell{fileIndsCell{cc}});
	S2P_centroids = vertcat(S2P_centCell{fileIndsCell{cc}});
	S2P_ints = S2P_intCell(fileIndsCell{cc});
    S2P_nucInts = S2P_nucIntCell(fileIndsCell{cc});
    S2P_nucVols = horzcat(S2P_nucVolCell{fileIndsCell{cc}})';
   

	sortedS2PNumCell{cc} = S2P_nums;
	sortedS2PVolCell{cc} = S2P_vols;
    sortedS2PVolPerNucCell{cc} = S2P_perNucVols;
    sortedS2PSolCell{cc} = S2P_sols;
	sortedS2PEloCell{cc} = S2P_elos;
	sortedS2PCentralSliceCell{cc} = S2P_slices;
	sortedS2PCentroidsCell{cc} = S2P_centroids;
    sortedS2PNucVolCell{cc} = vertcat(S2P_nucVols{:});
    sortedS2PNucClustVolCell{cc} = horzcat(S2P_nucClustVolCell{:})';
	
	for qq = 1:numQuantChannels
		
		sortedNucIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)nuc_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedNucIntCell{qq}{cc} = vertcat(sortedNucIntCell{qq}{cc}{:})';
		
		sortedCytoIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)cyto_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedCytoIntCell{qq}{cc} = vertcat(sortedCytoIntCell{qq}{cc}{:})';
		
        sortedNucStdCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)nuc_stdCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedNucStdCell{qq}{cc} = vertcat(sortedNucStdCell{qq}{cc}{:})';
		
		sortedS5PIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)S5P_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedS5PIntCell{qq}{cc} = vertcat(sortedS5PIntCell{qq}{cc}{:});

        sortedS5PNucIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)S5P_nucIntCell{ind}{qq},....
            find(fileIndsCell{cc}),...
            'UniformOutput',false));
        sortedS5PNucIntCell{qq}{cc} = vertcat(sortedS5PNucIntCell{qq}{cc}{:});

        
		sortedS2PIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)S2P_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedS2PIntCell{qq}{cc} = vertcat(sortedS2PIntCell{qq}{cc}{:});
		
        sortedS2PNucIntCell{qq}{cc} = vertcat(arrayfun(...
            @(ind)S2P_nucIntCell{ind}{qq},....
            find(fileIndsCell{cc}),...
            'UniformOutput',false));
        sortedS2PNucIntCell{qq}{cc} = vertcat(sortedS2PNucIntCell{qq}{cc}{:});

        
	end
	
end

%% Show reasoning for different volume cutoffs

Vol_threshold = 0.1;

% Define which number of data set corresponds to which condition
datasetInds = {[11,5,9,1,3,7,11+1,5+1,9+1,1+1,3+1,7+1],...
    [11+14,5+14,9+14,1+14,3+14,7+14,11+15,5+15,9+15,1+15,3+15,7+15]};

datasetNames = {'-LIF','RHB'};

refDatasetInds = [...
    1,2,1,2,1,2,1,2,1,2,1,2,...
    15,16,15,16,15,16,15,16,15,16,15,16];

targets = {...
    [1,2,3,4,5,6,1,2,3,4,5,6],...
    [1,2,3,4,5,6,1,2,3,4,5,6]};
targetNames = {'Ctrl 24 h','3 h','6 h','12 h','24 h','48 h'};

numPlots = numel(datasetInds);

% Make density plots of cluster volume

figure(2)
clf

for pp = 1:numPlots
    
    numTargets = numel(unique(targets{pp}));

    % Collect all data needed to analyze this data set
    inclDatasets = datasetInds{pp};
    numDatasets = numel(inclDatasets);
    
    % Plot for different differentiation time points

    collect_S5P = [];
    collect_S2P = [];
    collect_Vol = [];
    collect_S2P_S5P = [];
    collect_S2P_S2P = [];
    collect_S2P_Vol = [];
    
    for kk = 1:numDatasets
        
        cc = inclDatasets(kk);
        tt = targets{pp}(kk);
        
        inclInds = ...
            sortedS5PVolCell{cc}>=0;
        
        S5P_Num_vals = [sortedS5PNumCell{cc}];
        S2P_Num_vals = [sortedS2PNumCell{cc}];
        Nuc_S5P_vals = [sortedNucIntCell{1}{cc}-sortedCytoIntCell{1}{cc}]';
        Nuc_S2P_vals = [sortedNucIntCell{2}{cc}-sortedCytoIntCell{2}{cc}]';
        
        S5P_S5P_vals = [sortedS5PIntCell{1}{cc}(inclInds)];
        S5P_S2P_vals = [sortedS5PIntCell{2}{cc}(inclInds)];
        S5P_NucS5P_vals = [sortedS5PNucIntCell{1}{cc}(inclInds)];
        S5P_NucS2P_vals = [sortedS5PNucIntCell{2}{cc}(inclInds)];
        S5P_Vol_vals = [sortedS5PVolCell{cc}(inclInds)];
        S5P_Elo_vals = [sortedS5PEloCell{cc}(inclInds)];
        S5P_Sol_vals = [sortedS5PSolCell{cc}(inclInds)];
        
        S2P_S5P_vals = [sortedS2PIntCell{1}{cc}];
        S2P_S2P_vals = [sortedS2PIntCell{2}{cc}];
        S2P_NucS5P_vals = [sortedS2PNucIntCell{1}{cc}];
        S2P_NucS2P_vals = [sortedS2PNucIntCell{2}{cc}];
        S2P_Vol_vals = [sortedS2PVolCell{cc}];
        S2P_Elo_vals = [sortedS2PEloCell{cc}];
        S2P_Sol_vals = [sortedS2PSolCell{cc}];

        collect_S5P = vertcat(collect_S5P,S5P_S5P_vals);
        collect_S2P = vertcat(collect_S2P,S5P_S2P_vals);
        collect_Vol = vertcat(collect_Vol,S5P_Vol_vals);
        collect_S2P_S5P = vertcat(collect_S2P_S5P,S2P_S5P_vals);
        collect_S2P_S2P = vertcat(collect_S2P_S2P,S2P_S2P_vals);
        collect_S2P_Vol = vertcat(collect_S2P_Vol,S2P_Vol_vals);

    end

    VV_sample_points = linspace(0,0.5,200);

    subplot(numPlots,4,4.*(pp-1)+1)
    cla
    [VV_prob,VV_supp] = ksdensity(collect_Vol,VV_sample_points,...
        'Support','unbounded','Bandwidth',0.005);
    plot(VV_supp,VV_prob,'k-','LineWidth',1)
    hold on
    
    [VV_prob,VV_supp] = ksdensity(collect_S2P_Vol,VV_sample_points,...
        'Support','unbounded','Bandwidth',0.005);
    plot(VV_supp,VV_prob,'r--','LineWidth',1)
    
    set(gca,'YScale','log','YLim',[0.01,100],'XLim',[0,0.5])
    title(sprintf('%s',datasetNames{pp}),...
        'FontWeight','normal')

    xlabel('Volume [\mum^3]')
    ylabel('Prob. density')
    
    plot([1,1].*Vol_threshold,[0.001,1000],'b:',...
        'LineWidth',1.5)

    legend('Pol II Ser5P clusters','Pol II Ser2P spots',...
        sprintf('V_{threshold}=%2.2f \\mum^3',Vol_threshold))

    
    subplot(numPlots,4,4.*(pp-1)+2)
    
%     scatter1 = scatter(...
%         collect_S5P(collect_Vol>=Vol_threshold),...
%         collect_S2P(collect_Vol>=Vol_threshold),2,...
%         'MarkerFaceColor',[1.0,0,0],...
%         'MarkerEdgeColor','none');
%     alpha(scatter1,.05)
%     hold on
    scatter1 = scatter(...
        collect_S5P(collect_Vol<Vol_threshold),...
        collect_S2P(collect_Vol<Vol_threshold),2,...
        'MarkerFaceColor',[0.0,0.0,0.0],...
        'MarkerEdgeColor','none');
    alpha(scatter1,.05)
    hold off

    legend('V<V_{threshold}')
    xlabel('Pol II Ser5P Int.')
    ylabel('Pol II Ser2P Int.')
    set(gca,'XLim',[0.75,4],'YLim',[0,6],'Box','on')

    title(sprintf('%s, small Pol II Ser5P clusters',datasetNames{pp}),...
            'FontWeight','normal')

    subplot(numPlots,4,4.*(pp-1)+3)

    scatter1 = scatter(...
        collect_S5P(collect_Vol<Vol_threshold),...
        collect_S2P(collect_Vol<Vol_threshold),2,...
        'MarkerFaceColor',[0.5,0.5,0.5],...
        'MarkerEdgeColor','none');
    alpha(scatter1,.05)
    hold on
    scatter1 = scatter(...
        collect_S5P(collect_Vol>=Vol_threshold),...
        collect_S2P(collect_Vol>=Vol_threshold),2,...
        'MarkerFaceColor',[1,0,0],...
        'MarkerEdgeColor','none');
    alpha(scatter1,0.2)
    hold off

    legend('V<V_{threshold}','V\geqV_{threshold}')
    xlabel('Pol II Ser5P Int.')
    ylabel('Pol II Ser2P Int.')

    set(gca,'XLim',[0.75,4],'YLim',[0,6],'Box','on')

    title(sprintf('%s, large Pol II Ser5P clusters',datasetNames{pp}),...
            'FontWeight','normal')


  
    subplot(numPlots,4,4.*(pp-1)+4)
    plot(collect_S5P,collect_S2P,'ko',...
        'MarkerFaceColor',[1,0.5,0.5],...
        'MarkerEdgeColor','none',...
        'MarkerSize',2)
    hold on
    plot(collect_S2P_S5P,collect_S2P_S2P,'ro',...
        'MarkerFaceColor',[0,0,0],...
        'MarkerEdgeColor','none',...
        'MarkerSize',2)
    hold on
    legend('Ser5P clusters','Ser2P spots')
    xlabel('Pol II Ser5P Int.')
    ylabel('Pol II Ser2P Int.')
    set(gca,'XLim',[0,5],'YLim',[0,30])

    title(sprintf('%s, Pol II Ser2P spots',datasetNames{pp}),...
            'FontWeight','normal')    



end


%% Analyze the results pooled from both repeats

Vol_threshold = 0.1;

minNucVol = 0; % Miminal nuclear volume
maxNucVol = Inf; % Maxinal nuclear volume

XLimVec = [0,30];

datasetNames = {'-LIF','RHB'};

datasetInds = {[11,5,9,1,3,7,11+1,5+1,9+1,1+1,3+1,7+1],...
    [11+14,5+14,9+14,1+14,3+14,7+14,11+15,5+15,9+15,1+15,3+15,7+15]};

refDatasetInds = [...
    1,2,1,2,1,2,1,2,1,2,1,2,...
    15,16,15,16,15,16,15,16,15,16,15,16];

targets = {...
    [1,2,3,4,5,6,1,2,3,4,5,6],...
    [1,2,3,4,5,6,1,2,3,4,5,6]};
targetNames = {'Ctrl 24 h','3 h','6 h','12 h','24 h','48 h'};

mapLimits = {[0,6,0,3],[0,8,0,6]};

S5P_threshold = 0.0;

n_boot = 50;%1000;


numPlots = numel(datasetInds);

plotStyles = {'k-','r--'};

figure(1)
clf

figure(2)
clf

for pp = 1:numPlots
    
    numTargets = numel(unique(targets{pp}));

    % Collect all data needed to analyze this data set
    inclDatasets = datasetInds{pp};
    numDatasets = numel(inclDatasets);
    
    % Collect for distribution plots and normalization
    
    Nuc_group_vec = [];
    Nuc_Vol_vec = [];
    Nuc_S5P_vec = [];
    Nuc_S2P_vec = [];
    Nuc_ClusterNum_vec_unscaled = [];
    Nuc_ClusterNum_vec = [];


    Num_S5P_vec = [];
    Num_S2P_vec = [];
        
    S5P_group_vec = [];
    S5P_S5P_vec = [];
    S5P_S2P_vec = [];
    S5P_NucClusterNum_vec = [];
    S5P_NucS5P_vec = [];
    S5P_NucS2P_vec = [];
    S5P_Vol_vec = [];
    S5P_Sol_vec = [];
    S5P_Elo_vec = [];

    S5P_all_group_vec = [];
    S5P_allVol_vec = [];
    
    S2P_group_vec = [];
    S2P_S5P_vec = [];
    S2P_S2P_vec = [];
    S2P_NucS5P_vec = [];
    S2P_NucS2P_vec = [];

    for kk = 1:numDatasets
        
        cc = inclDatasets(kk);
        tt = targets{pp}(kk);
        
        pooledNucVols = [sortedS5PNucVolCell{cc}];
%        pooledNucVols = vertcat(pooledNucVols{:});

        inclInds = ...
            sortedS5PVolCell{cc}>=Vol_threshold ...
            & sortedS5PIntCell{1}{cc}>=S5P_threshold ...
            & pooledNucVols>=minNucVol ...
            & pooledNucVols<=maxNucVol;
        
        Nuc_Vol_vals = [sortedNucVolCell{cc}];
        S5P_Num_vals = [sortedS5PNumCell{cc}];
        S2P_Num_vals = [sortedS2PNumCell{cc}];
        Nuc_S5P_vals = [sortedNucIntCell{1}{cc}-sortedCytoIntCell{1}{cc}]';
        Nuc_S2P_vals = [sortedNucIntCell{2}{cc}-sortedCytoIntCell{2}{cc}]';

        NucClusterNum_vals = cellfun( ...
            @(elmt)sum(elmt>=Vol_threshold),...
            sortedS5PVolPerNucCell{cc});
        NucClusterNum_vals_unscaled = NucClusterNum_vals;
        NucClusterNum_vals = NucClusterNum_vals./Nuc_Vol_vals;

        % ---- Calculate clusters per volume over the whole nucleus, but
        % saved for every individual S5P cluster

        % Volume of all clusters, per S5P cluster
        ClustVols = vertcat(sortedS5PNucClustVolCell{cc}{:});
        %ClustVols = vertcat(ClustVols{:});
        S5P_NucClusterNum_vals = cellfun(...
            @(elmt)sum(elmt>=Vol_threshold),ClustVols);
        S5P_NucClusterNum_vals = S5P_NucClusterNum_vals(inclInds);

        % Volume of the nuclei, per S5P cluster
        NucVol = pooledNucVols(inclInds);

        S5P_NucClusterNum_vals = ...
            S5P_NucClusterNum_vals./NucVol;
        % ----




        S5P_S5P_vals = [sortedS5PIntCell{1}{cc}(inclInds)];
        S5P_S2P_vals = [sortedS5PIntCell{2}{cc}(inclInds)];
        S5P_NucS5P_vals = [sortedS5PNucIntCell{1}{cc}(inclInds)];
        S5P_NucS2P_vals = [sortedS5PNucIntCell{2}{cc}(inclInds)];
        S5P_allVol_vals = [sortedS5PVolCell{cc}];
        S5P_Vol_vals = [sortedS5PVolCell{cc}(inclInds)];
        S5P_Elo_vals = [sortedS5PEloCell{cc}(inclInds)];
        S5P_Sol_vals = [sortedS5PSolCell{cc}(inclInds)];
        
        S2P_S5P_vals = [sortedS2PIntCell{1}{cc}];
        S2P_S2P_vals = [sortedS2PIntCell{2}{cc}];
        S2P_NucS5P_vals = [sortedS2PNucIntCell{1}{cc}];
        S2P_NucS2P_vals = [sortedS2PNucIntCell{2}{cc}];
        S2P_Vol_vals = [sortedS2PVolCell{cc}];
        S2P_Elo_vals = [sortedS2PEloCell{cc}];
        S2P_Sol_vals = [sortedS2PSolCell{cc}];
        
        Nuc_group_vec = vertcat(Nuc_group_vec,tt.*ones(size(Nuc_S5P_vals)));
        Nuc_Vol_vec = vertcat(Nuc_Vol_vec,Nuc_Vol_vals);
        Nuc_S5P_vec = vertcat(Nuc_S5P_vec,Nuc_S5P_vals);
        Nuc_S2P_vec = vertcat(Nuc_S2P_vec,Nuc_S2P_vals);
        Num_S5P_vec = vertcat(Num_S5P_vec,S5P_Num_vals);
        Num_S2P_vec = vertcat(Num_S2P_vec,S5P_Num_vals);
        Nuc_ClusterNum_vec_unscaled = ...
            vertcat(Nuc_ClusterNum_vec_unscaled,...
            NucClusterNum_vals_unscaled);
        Nuc_ClusterNum_vec = ...
            vertcat(Nuc_ClusterNum_vec,NucClusterNum_vals);

        S5P_group_vec = vertcat(S5P_group_vec,tt.*ones(size(S5P_S5P_vals)));
        S5P_S5P_vec = vertcat(S5P_S5P_vec,S5P_S5P_vals);
        S5P_S2P_vec = vertcat(S5P_S2P_vec,S5P_S2P_vals);
        S5P_NucClusterNum_vec = vertcat(S5P_NucClusterNum_vec,...
            S5P_NucClusterNum_vals);
        S5P_NucS5P_vec = vertcat(S5P_NucS5P_vec,S5P_NucS5P_vals);
        S5P_NucS2P_vec = vertcat(S5P_NucS2P_vec,S5P_NucS2P_vals);
        S5P_Vol_vec = vertcat(S5P_Vol_vec,S5P_Vol_vals);
        S5P_Elo_vec = vertcat(S5P_Elo_vec,S5P_Elo_vals);
        S5P_Sol_vec = vertcat(S5P_Sol_vec,S5P_Sol_vals);
        
        S5P_all_group_vec = vertcat(S5P_all_group_vec,...
            tt.*ones(size(S5P_allVol_vals)));
        S5P_allVol_vec = vertcat(S5P_allVol_vec,S5P_allVol_vals);

        S2P_group_vec = vertcat(S2P_group_vec,tt.*ones(size(S2P_S5P_vals)));
        S2P_S5P_vec = vertcat(S2P_S5P_vec,S2P_S5P_vals);
        S2P_S2P_vec = vertcat(S2P_S2P_vec,S2P_S2P_vals);
        S2P_NucS5P_vec = vertcat(S2P_NucS5P_vec,S2P_NucS5P_vals);
        S2P_NucS2P_vec = vertcat(S2P_NucS2P_vec,S2P_NucS2P_vals);
                
    end

    Nuc_ClusterNum_vec = Nuc_ClusterNum_vec.*median(Nuc_Vol_vec);
    S5P_NucClusterNum_vec = S5P_NucClusterNum_vec.*median(Nuc_Vol_vec);
    
    Num_S5P_mean = zeros(1,numTargets);
    Num_S5P_CI = zeros(2,numTargets);
    Nuc_S5P_median = zeros(1,numTargets);
    Nuc_S5P_CI = zeros(2,numTargets);
    Nuc_S2P_median = zeros(1,numTargets);
    Nuc_S2P_CI = zeros(2,numTargets);
    Nuc_ClusterNum_mean = zeros(1,numTargets);
    Nuc_ClusterNum_CI = zeros(2,numTargets);
    Nuc_Vol_median = zeros(1,numTargets);
    Nuc_Vol_CI = zeros(2,numTargets);
    Nuc_Sol_median = zeros(1,numTargets);
    Nuc_Sol_CI = zeros(2,numTargets);

    S5P_allVol_mean = zeros(1,numTargets);
    S5P_allVol_CI = zeros(2,numTargets);
     
    for tt = 1:numTargets
        
        targetInclInds = Nuc_group_vec==tt;
        
        Num_S5P_mean(tt) = mean(Num_S5P_vec(targetInclInds));
        Num_S5P_CI(:,tt) = ...
            bootci(n_boot,@mean,Num_S5P_vec(targetInclInds));
        
        Nuc_S5P_median(tt) = median(Nuc_S5P_vec(targetInclInds));
        Nuc_S5P_CI(:,tt) = ...
            bootci(n_boot,@median,Nuc_S5P_vec(targetInclInds));
        
        Nuc_S2P_median(tt) = median(Nuc_S2P_vec(targetInclInds));
        Nuc_S2P_CI(:,tt) = ...
            bootci(n_boot,@median,Nuc_S2P_vec(targetInclInds));
                
        Nuc_ClusterNum_mean(tt) = mean(Nuc_ClusterNum_vec_unscaled(targetInclInds));
        Nuc_ClusterNum_CI(:,tt) = ...
           bootci(n_boot,@mean,Nuc_ClusterNum_vec_unscaled(targetInclInds));

        targetInclInds = S5P_group_vec==tt;
        Nuc_Vol_median(tt) = median(S5P_Vol_vec(targetInclInds));
        Nuc_Vol_CI(:,tt) = ...
            bootci(n_boot,@median,S5P_Vol_vec(targetInclInds));

        Nuc_Sol_median(tt) = median(S5P_Sol_vec(targetInclInds));
        Nuc_Sol_CI(:,tt) = ...
            bootci(n_boot,@median,S5P_Sol_vec(targetInclInds));


    end
        
    Nuc_S5P_median = Nuc_S5P_median./median(Nuc_S5P_vec);
    Nuc_S5P_CI = Nuc_S5P_CI./median(Nuc_S5P_vec);
    Nuc_S2P_median = Nuc_S2P_median./median(Nuc_S2P_vec);
    Nuc_S2P_CI = Nuc_S2P_CI./median(Nuc_S2P_vec);
    
    Nuc_S5P_vec = Nuc_S5P_vec./median(Nuc_S5P_vec);
    Nuc_S2P_vec = Nuc_S2P_vec./median(Nuc_S2P_vec);
    
    S5P_S5P_vec = S5P_S5P_vec./median(S5P_S5P_vec);
    S5P_S2P_vec = S5P_S2P_vec./median(S5P_S2P_vec);
    S5P_NucS5P_vec = S5P_NucS5P_vec./median(S5P_NucS5P_vec);
    S5P_NucS2P_vec = S5P_NucS2P_vec./median(S5P_NucS2P_vec);
    
    S2P_S5P_vec = S2P_S5P_vec./median(S2P_S5P_vec);
    S2P_S2P_vec = S2P_S2P_vec./median(S2P_S2P_vec);
    S2P_NucS5P_vec = S2P_NucS5P_vec./median(S5P_NucS5P_vec);
    S2P_NucS2P_vec = S2P_NucS2P_vec./median(S5P_NucS2P_vec);
    
    figure(1)
    numDisplays = 6;

    % Nucleus-level results

    subplot(numPlots.*2,numDisplays,(pp-1).*numDisplays.*2+1)

    % distributionPlot(Nuc_S5P_vec,'groups',Nuc_group_vec,...
    %     'xNames',targetNames,'color',[0.6,0.6,0.6],...
    %     'showMM',0,'addSpread',0)
    plotSpread(Nuc_S5P_vec,'distributionIdx',Nuc_group_vec,...
        'xNames',targetNames,'distributionColors',[0.6,0.6,0.6])
    hold on
    
    errorbar(1:numTargets,Nuc_S5P_median,...
        Nuc_S5P_CI(1,:)-Nuc_S5P_median,...
        Nuc_S5P_median-Nuc_S5P_CI(2,:),...
        'k-o','MarkerFaceColor',[0,0,0],...
        'LineWidth',1);

    set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
        'XTick',1:numTargets,'XTickLabel',targetNames,...
        'XLim',[0.3,numTargets+0.8])

    ylabel('Pol II Ser5P Int.')

    title(datasetNames{pp})



    subplot(numPlots.*2,numDisplays,(pp-1).*numDisplays.*2+2)

    % distributionPlot(Nuc_S2P_vec,'groups',Nuc_group_vec,...
    %     'xNames',targetNames,'color',[0.6,0.6,0.6],...
    %     'showMM',0,'addSpread',0)
    plotSpread(Nuc_S2P_vec,'distributionIdx',Nuc_group_vec,...
        'xNames',targetNames,'distributionColors',[0.6,0.6,0.6])
    hold on
    
    errorbar(1:numTargets,Nuc_S2P_median,...
        Nuc_S2P_CI(1,:)-Nuc_S2P_median,...
        Nuc_S2P_median-Nuc_S2P_CI(2,:),...
        'k-o','MarkerFaceColor',[0,0,0],...
        'LineWidth',1);

    set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
        'XTick',1:numTargets,'XTickLabel',targetNames,...
        'XLim',[0.3,numTargets+0.8])

    ylabel('Pol II Ser2P Int.')



    subplot(numPlots.*2,numDisplays,(pp-1).*numDisplays.*2+3)

    % distributionPlot(Nuc_ClusterNum_vec_unscaled,'groups',Nuc_group_vec,...
    %     'xNames',targetNames,'color',[0.6,0.6,0.6],...
    %     'showMM',0,'addSpread',0)
    plotSpread(Nuc_ClusterNum_vec,'distributionIdx',Nuc_group_vec,...
        'xNames',targetNames,'distributionColors',[0.6,0.6,0.6])
    hold on
    
    errorbar(1:numTargets,Nuc_ClusterNum_mean,...
        Nuc_ClusterNum_CI(1,:)-Nuc_ClusterNum_mean,...
        Nuc_ClusterNum_mean-Nuc_ClusterNum_CI(2,:),...
        'k-o','MarkerFaceColor',[0,0,0],...
        'LineWidth',1);

    set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
        'XTick',1:numTargets,'XTickLabel',targetNames,...
        'XLim',[0.3,numTargets+0.8])

    ylabel('Number large clusters')


    subplot(numPlots.*2,numDisplays,(pp-1).*numDisplays.*2+4)

    % distributionPlot(S5P_Vol_vec,'groups',S5P_group_vec,...
    %     'xNames',targetNames,'color',[0.6,0.6,0.6],...
    %     'showMM',0,'addSpread',0)
    % hold on
    
    errorbar(1:numTargets,Nuc_Vol_median,...
        Nuc_Vol_CI(1,:)-Nuc_Vol_median,...
        Nuc_Vol_median-Nuc_Vol_CI(2,:),...
        'k-o','MarkerFaceColor',[1,1,1],...
        'LineWidth',1);

    set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
        'XTick',1:numTargets,'XTickLabel',targetNames,...
        'XLim',[0.3,numTargets+0.8])

    ylabel('Cluster Volume [\mum^3]')

    subplot(numPlots.*2,numDisplays,(pp-1).*numDisplays.*2+5)
    % 
    % distributionPlot(S5P_Sol_vec,'groups',S5P_group_vec,...
    %     'xNames',targetNames,'color',[0.6,0.6,0.6],...
    %     'showMM',0,'addSpread',0)
    % hold on
    
    errorbar(1:numTargets,Nuc_Sol_median,...
        Nuc_Sol_CI(1,:)-Nuc_Sol_median,...
        Nuc_Sol_median-Nuc_Sol_CI(2,:),...
        'k-o','MarkerFaceColor',[1,1,1],...
        'LineWidth',1);

    set(gca,'XLim',[0.5,numTargets+0.5],'Box','on','YLim',[-Inf,+Inf],...
        'XTick',1:numTargets,'XTickLabel',targetNames,...
        'XLim',[0.3,numTargets+0.8])

    ylabel('Cluster Solidity')



    % --- running window plots of sorted per-nucleus results

    Nuc_Num = Nuc_ClusterNum_vec;
    Nuc_S5P = Nuc_S5P_vec;
    Nuc_S2P = Nuc_S2P_vec;
    
    numWindows = 60;
    minNucCount = 12;

    lowLim = 1;
    highLim = 30;
    numWindows = (highLim-lowLim+1).*4;
    windowCenters = linspace(lowLim,highLim,numWindows);
    windowWidth = 8.0;

    inWindowCount = zeros(1,numWindows);
    windowIndsCell = cell(1,numWindows);
    for ww = 1:numWindows
        windowEdges = windowCenters(ww)+[-0.5,+0.5].*windowWidth;
        windowInds = find(Nuc_Num>=windowEdges(1) & Nuc_Num<windowEdges(2));
        windowIndsCell{ww}=windowInds;
        inWindowCount(ww) = numel(windowInds);
    end


    subplot(numPlots.*2,numDisplays,((pp-1).*2+1).*numDisplays+1)
    
    mean_S5P = zeros(1,numWindows);
    CI_S5P = zeros(2,numWindows);
    for ww = 1:numWindows
        windowInds = windowIndsCell{ww};
        if inWindowCount(ww)>=minNucCount
            mean_S5P(ww) = median(Nuc_S5P(windowInds));
            CI_S5P(:,ww) = bootci(n_boot,@median,Nuc_S5P(windowInds));
            CI_S5P(:,ww) = CI_S5P(:,ww)-mean_S5P(ww);
        else
            mean_S5P(ww) = NaN;
            CI_S5P(:,ww) = NaN;
        end
    end
    
    err_h = errorbar(...
        windowCenters,mean_S5P,...
        CI_S5P(1,:),CI_S5P(2,:),...
        'ko','LineWidth',1,'MarkerSize',3,...
        'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    err_h.CapSize = 0;
    err_h.Color = [0.7,0.7,0.7];

    xlabel('Number large clusters')
    ylabel('Nucleus S5P')
    set(gca,'XDir','reverse','YLim',[0,4],'XLim',XLimVec)



    subplot(numPlots.*2,numDisplays,((pp-1).*2+1).*numDisplays+2)
    
    mean_S2P = zeros(1,numWindows);
    CI_S2P = zeros(2,numWindows);
    for ww = 1:numWindows
        windowInds = windowIndsCell{ww};
        if inWindowCount(ww)>=minNucCount
            mean_S2P(ww) = median(Nuc_S2P(windowInds));
            CI_S2P(:,ww) = bootci(n_boot,@median,Nuc_S2P(windowInds));
            CI_S2P(:,ww) = CI_S2P(:,ww)-mean_S2P(ww);
        else
            mean_S2P(ww) = NaN;
            CI_S2P(:,ww) = NaN;
        end
    end
    
    err_h = errorbar(...
        windowCenters,mean_S2P,...
        CI_S2P(1,:),CI_S2P(2,:),...
        'ko','LineWidth',1,'MarkerSize',3,...
        'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    err_h.CapSize = 0;
    err_h.Color = [0.7,0.7,0.7];

    xlabel('Number large clusters')
    ylabel('Nucleus S2P')
    set(gca,'XDir','reverse','YLim',[0,2],'XLim',XLimVec)

    figure(2)
    
    for tt = 1:numTargets

        subplot(numTargets,numPlots, ...
            numPlots.*(tt-1)+pp)

        xx_vec = linspace(0,50,400);
        [pp_val_ctrl,xx_val_ctrl] = ...
            ksdensity(Nuc_ClusterNum_vec(Nuc_group_vec==1),xx_vec,...
            'Support','nonnegative','BoundaryCorrection','reflection',...
            'Bandwidth',2.0);
        fill([0,xx_val_ctrl,0],...
            [0,pp_val_ctrl,0],...
            [0.65,0.65,0.65],'EdgeColor','none')
        hold on
        [pp_val,xx_val] = ...
            ksdensity(Nuc_ClusterNum_vec(Nuc_group_vec==tt),xx_vec,...
            'Support','nonnegative','BoundaryCorrection','reflection',...
            'Bandwidth',2.0);
        legend_handle = plot(xx_val,pp_val,'k-','LineWidth',1.5);
        plot(Nuc_ClusterNum_mean(tt).*[1,1],[0,0.1],'k:',...
            'LineWidth',1)
        text(8,0.08,...
            sprintf('%1.1f',Nuc_ClusterNum_mean(tt)))
        set(gca,'XLim',[0,35],'XDir','normal','YLim',[0,0.1])
        legend(legend_handle,targetNames(tt),'Location','Northeast')

        if tt == 1
            title(datasetNames{pp})
        end

        if tt == numTargets
            xlabel('Number large clusters')
        else
            set(gca,'XTickLabel',[])
        end

        ylabel('Probability density')

    end

    figure(1)

    % --- running window plots of sorted per-cluster results



    S5P_Num = S5P_NucClusterNum_vec(S5P_group_vec>0);
    S5P_S5P = S5P_S5P_vec(S5P_group_vec>0);
    S5P_S2P = S5P_S2P_vec(S5P_group_vec>0);
    S5P_Vol = S5P_Vol_vec(S5P_group_vec>0);
    S5P_Sol = S5P_Sol_vec(S5P_group_vec>0);

    % lowLim = round(prctile(S5P_Num,1));
    % highLim = round(prctile(S5P_Num,99));
    % windowCenters = linspace(lowLim,highLim,numWindows);
    % windowWidth = (highLim-lowLim)./2.5;

    minClustCount = 12;

    lowLim = round(prctile(S5P_Num,1));
    highLim = round(prctile(S5P_Num,99));
    numWindows = (highLim-lowLim+1).*4;
    windowCenters = linspace(lowLim,highLim,numWindows);
    windowWidth = 4.0;

    XLimVec = [0,highLim];

    inWindowCount = zeros(1,numWindows);
    for ww = 1:numWindows
        windowEdges = windowCenters(ww)+[-0.5,+0.5].*windowWidth;
        windowInds = find(S5P_Num>=windowEdges(1) & S5P_Num<windowEdges(2));
        windowIndsCell{ww}=windowInds;
        inWindowCount(ww) = numel(windowInds);
    end

    subplot(numPlots.*2,numDisplays,((pp-1).*2+1).*numDisplays+3)
    
    mean_S5P = zeros(1,numWindows);
    CI_S5P = zeros(2,numWindows);
    for ww = 1:numWindows
        windowInds = windowIndsCell{ww};
        if inWindowCount(ww)>=minClustCount
            mean_S5P(ww) = median(S5P_S5P(windowInds));
            CI_S5P(:,ww) = bootci(n_boot,@median,S5P_S5P(windowInds));
            CI_S5P(:,ww) = CI_S5P(:,ww)-mean_S5P(ww);
        else
            mean_S5P(ww) = NaN;
            CI_S5P(:,ww) = NaN;
        end
    end
    
    err_h = errorbar(...
        windowCenters,mean_S5P,...
        CI_S5P(1,:),CI_S5P(2,:),...
        'ko','LineWidth',1,'MarkerSize',3,...
        'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    err_h.CapSize = 0;
    err_h.Color = [0.7,0.7,0.7];

    xlabel('Number large clusters')
    ylabel('Cluster S5P')
    set(gca,'XDir','reverse','YLim',[0.8,1.3],'XLim',XLimVec)

    
    
    
    subplot(numPlots.*2,numDisplays,((pp-1).*2+1).*numDisplays+4)

    mean_S2P = zeros(1,numWindows);
    CI_S2P = zeros(2,numWindows);
    for ww = 1:numWindows
        windowInds = windowIndsCell{ww};
        if inWindowCount(ww)>=minClustCount
            mean_S2P(ww) = median(S5P_S2P(windowInds));
            CI_S2P(:,ww) = bootci(n_boot,@median,S5P_S2P(windowInds));
            CI_S2P(:,ww) = CI_S2P(:,ww)-mean_S2P(ww);
        else
            mean_S2P(ww) = NaN;
            CI_S2P(:,ww) = NaN;
        end
        
    end
    
    err_h = errorbar(...
        windowCenters,mean_S2P,CI_S2P(1,:),CI_S2P(2,:),...
        'ko','LineWidth',1,'MarkerSize',3,...
        'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    err_h.CapSize = 0;
    err_h.Color = [0.7,0.7,0.7];

    xlabel('Number large clusters')
    ylabel('Cluster S2P')
    set(gca,'XDir','reverse','YLim',[0.9,1.2],'XLim',XLimVec)

    
    subplot(numPlots.*2,numDisplays,((pp-1).*2+1).*numDisplays+5)

    mean_Vol = zeros(1,numWindows);
    CI_Vol = zeros(2,numWindows);
    for ww = 1:numWindows
        windowInds = windowIndsCell{ww};
        if inWindowCount(ww)>=minClustCount
            mean_Vol(ww) = median(S5P_Vol(windowInds));
            CI_Vol(:,ww) = bootci(n_boot,@median,S5P_Vol(windowInds));
            CI_Vol(:,ww) = CI_Vol(:,ww)-mean_Vol(ww);
        else
            mean_Vol(ww) = NaN;
            CI_Vol(:,ww) = NaN;
        end
        
    end
    
    err_h = errorbar(...
        windowCenters,mean_Vol,CI_Vol(1,:),CI_Vol(2,:),...
        'ko','LineWidth',1,'MarkerSize',3,...
        'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    err_h.CapSize = 0;
    err_h.Color = [0.7,0.7,0.7];

    xlabel('Number large clusters')
    ylabel('Volume [\mum^3]')
    set(gca,'XDir','reverse','YLim',[0.125,0.15],'XLim',XLimVec)

    
    subplot(numPlots.*2,numDisplays,((pp-1).*2+1).*numDisplays+6)

    mean_Sol = zeros(1,numWindows);
    CI_Sol = zeros(2,numWindows);
    for ww = 1:numWindows
        windowInds = windowIndsCell{ww};
        if inWindowCount(ww)>=minClustCount
            mean_Sol(ww) = median(S5P_Sol(windowInds));
            CI_Sol(:,ww) = bootci(n_boot,@median,S5P_Sol(windowInds));
            CI_Sol(:,ww) = CI_Sol(:,ww)-mean_Sol(ww);
        else
            mean_Sol(ww) = NaN;
            CI_Sol(:,ww) = NaN;
        end
    end
    
    err_h = errorbar(...
        windowCenters,mean_Sol,CI_Sol(1,:),CI_Sol(2,:),...
        'ko','LineWidth',1,'MarkerSize',3,...
        'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    err_h.CapSize = 0;
    err_h.Color = [0.7,0.7,0.7];

    set(gca,'XDir','reverse','YLim',[-Inf,Inf],'XLim',XLimVec)

    ylabel('Solidity')
    xlabel('Number large clusters')

end



