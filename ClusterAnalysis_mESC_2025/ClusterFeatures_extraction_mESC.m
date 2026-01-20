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
save_file = fullfile(".", "NuclearClusterFeatures", "ClusterFeatures.mat");

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
