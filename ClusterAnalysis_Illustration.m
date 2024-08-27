clear all

% Specify the directory that contains the extraced files from the last
% step, where you extracted from the raw files obtained from the microscope
sourceDirectory = './ExtractedStacks/**/';

Vol_threshold = 0.1;

% Example nuclei
exampleFileInds =    [73,131,160,169,203,237]; 
exampleNucleusInds = [3, 2,  1,  1,  3,  1];
exampleClusterInds = [15,16, 2,  6,  6,  3];

exampleFileInds = exampleFileInds([1,2,3,4,6]);
exampleNucleusInds = exampleNucleusInds([1,2,3,4,6]);
exampleClusterInds = exampleClusterInds([1,2,3,4,6]);

% exampleFileInds = 237; 231; 229;
% exampleNucleusInds = 1; 3; 3;
% exampleClusterInds = 3; 2; 4;

% Targets:
% Ctrl 24 h, 25 clusters: 73, 3, 15
% -LIF 3 h, 25 clusters; 131, 2, 16; 126, 1,14
% -LIF 6 h, 19/20 clusters, high Ser2P; 160, 1, 2
% -LIF 12 h, 15 clusters, high solidity; 169, 1, 6
% -LIF 24 h, 10 clusters, normal solidity; 203, 3, 6; 203, 1, 7 (exclude)
% -LIF 48 h, 5 clusters, normal solidity; 229, 1, 4

S5P_lims = [90,545];
S2P_lims = [160,1200];
S5P_zoom_lims = [0.6,3.1];
S2P_zoom_lims = [0.6,3.0];

numExamples = numel(exampleFileInds);

% Channels for segmentation
NucSegChannel = 1; % Channel used to detect nuclei
S5P_SegChannel = 1; % Channel used to detect Pol II S5P clusters
S2P_SegChannel = 2; % Channel used to detect Pol II S2P clusters

% Save images of the clusters
ImgSquareExtension = 14; % pixels for cut-out image extension, set 0 for no images
% Which image channels to store in example images
storeImgChannels = [1,2];
numStoreChannels = numel(storeImgChannels);

% Target channels for intensity quantification, applied for all objects
quantChannels = [1,2];
quantBlurSigma = [0,0.15];
numQuantChannels = numel(quantChannels);

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


figure(1)
clf

for ff = 1:numExamples
	
	fprintf('Processing example %d of %d\n',ff,numExamples)
	
    exampleInd = exampleFileInds(ff);

	thisCondInd = condInds(exampleInd);	
	thisFilePath = listing(exampleInd).name;
    thisCondName = condNames{exampleInd};
	
	loadStruct = load(thisFilePath,...
		'imgStack','imgSize','pixelSize','zStepSize');
	imgStack = loadStruct.imgStack;
	imgSize = loadStruct.imgSize;
	pixelSize = loadStruct.pixelSize;
	zStepSize = loadStruct.zStepSize;

	% Nuclei segmentation
	segImg = imgStack{NucSegChannel};
    if nuc_segBlurSigma_nucleus>0
        segImg = ...
            + imgaussfilt(segImg,nuc_segBlurSigma_nucleus./pixelSize) ...
            - imgaussfilt(segImg,nuc_segBlurSigma_BG_removal./pixelSize);
    else
        segImg = ...
            + segImg ...
            - imgaussfilt(segImg,nuc_segBlurSigma_BG_removal./pixelSize);
    end

	[bin_counts,bin_centers] = hist(segImg(:),1000);
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
    
		
	% --- Connected component segmentation of nuclei
	comps = bwconncomp(NucSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	minPixels = Nuc_min_vol./(pixelSize.^2)./zStepSize;
	comps.NumObjects = sum(numPxls>=minPixels);
	comps.PixelIdxList = comps.PixelIdxList(numPxls>=minPixels);
	numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
	
	props = regionprops3(comps,imgStack{NucSegChannel},...
		'Solidity','VoxelValues');
	
    if comps.NumObjects>0
        Solidity_array = [props.Solidity];
        CoV_array = ...
            cellfun(@(vals)std(vals(:))./mean(vals(:)),...
            props.VoxelValues);
        inclNucInds = Solidity_array>=Nuc_min_sol ...
            & CoV_array>=Nuc_min_CoV;
        comps.NumObjects = sum(Solidity_array>=Nuc_min_sol);
        comps.PixelIdxList = comps.PixelIdxList(Solidity_array>=Nuc_min_sol);
        numPxls = cellfun(@(elmt)numel(elmt),comps.PixelIdxList);
    end
    
	numNuclei = comps.NumObjects;
	
    if comps.NumObjects>0
                
        nuc_intCell = cell(1,numQuantChannels);
        cyto_intCell = cell(1,numQuantChannels);
        nuc_stdCell = cell(1,numQuantChannels);
        for qq = 1:numQuantChannels
            quantImg = imgStack{quantChannels(qq)};
            quantProps = regionprops3(comps,quantImg,...
                'MeanIntensity','VoxelIdxList','VoxelValues');
            nuc_intCell{qq} = [quantProps.MeanIntensity];
            cyto_intCell{qq} = zeros(numNuclei,1);
            nuc_stdCell{qq} = cellfun(@std,quantProps.VoxelValues);
            
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
                cyto_intCell{qq}(nn) = mean(quantImg(cytoMask));
            end
        end
        
        
        props = regionprops3(comps,imgStack{NucSegChannel},...
            'Volume','VoxelValues','Solidity','VoxelIdxList',...
            'BoundingBox');
        
        Volume_array = [props.Volume].*pixelSize.^2.*zStepSize;
        Intensity_array = cellfun(@(vals)median(vals),props.VoxelValues);
        Solidity_array = [props.Solidity];

        nuc_volCell = Volume_array;
        
        
        % --- For each nucleus, get objects from the different channels
        
        nn = exampleNucleusInds(ff);

        boxArray = props.BoundingBox(nn,:);

        S5P_subImage = imgStack{S5P_SegChannel}(...
            boxArray(2)-20+0.5:boxArray(2)+boxArray(5)+20-0.5,...
            boxArray(1)-20+0.5:boxArray(1)+boxArray(4)+20-0.5,...
            boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);

        S2P_subImage = imgStack{S2P_SegChannel}(...
            boxArray(2)-20+0.5:boxArray(2)+boxArray(5)+20-0.5,...
            boxArray(1)-20+0.5:boxArray(1)+boxArray(4)+20-0.5,...
            boxArray(3)+0.5:boxArray(3)+boxArray(6)-0.5);

        S5P_subImage_raw = imgaussfilt(S5P_subImage,...
            S5P_segBlurSigma_object./pixelSize);
        S2P_subImage_raw = imgaussfilt(S2P_subImage,...
            S2P_segBlurSigma_object./pixelSize);

        S5P_subImage = ...
            + imgaussfilt(S5P_subImage,S5P_segBlurSigma_object./pixelSize) ...
            - imgaussfilt(S5P_subImage,S5P_segBlurSigma_BG_removal./pixelSize);
        S2P_subImage = ...
            + imgaussfilt(S2P_subImage,S2P_segBlurSigma_object./pixelSize) ...
            - imgaussfilt(S2P_subImage,S2P_segBlurSigma_BG_removal./pixelSize);

        NucMask = false(imgSize);
        NucMask(props.VoxelIdxList{nn}) = true;
        NucMask_subImage = NucMask(...
            boxArray(2)-20+0.5:boxArray(2)+boxArray(5)+20-0.5,...
            boxArray(1)-20+0.5:boxArray(1)+boxArray(4)+20-0.5,...
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
                    boxArray(2)-20+0.5:boxArray(2)+boxArray(5)+20-0.5,...
                    boxArray(1)-20+0.5:boxArray(1)+boxArray(4)+20-0.5,...
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
        S5P_numPxls = cellfun(@(elmt)numel(elmt),S5P_comps.PixelIdxList);

        S2P_comps = bwconncomp(S2P_mask,18);
        S2P_numPxls = cellfun(@(elmt)numel(elmt),S2P_comps.PixelIdxList);
        minPixels = S2P_minVol./(pixelSize.^2)./zStepSize;
        S2P_comps.NumObjects = sum(S2P_numPxls>minPixels);
        S2P_comps.PixelIdxList = S2P_comps.PixelIdxList(S2P_numPxls>minPixels);
        S2P_numPxls = cellfun(@(elmt)numel(elmt),S2P_comps.PixelIdxList);

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

            S5P_numPxls = cellfun(@(elmt)numel(elmt),...
                S5P_comps.PixelIdxList);

            S5P_comps_small_inds = ...
                S5P_numPxls.*pixelSize.^2.*zStepSize<Vol_threshold;
            S5P_comps_small = S5P_comps;
            S5P_comps_small.PixelIdxList = ...
                S5P_comps_small.PixelIdxList(S5P_comps_small_inds);
            S5P_comps_small.NumObjects = sum(S5P_comps_small_inds);

            S5P_comps_large_inds = ...
                S5P_numPxls.*pixelSize.^2.*zStepSize>=Vol_threshold;
            S5P_comps_large = S5P_comps;
            S5P_comps_large.PixelIdxList = ...
                S5P_comps_large.PixelIdxList(S5P_comps_large_inds);
            S5P_comps_large.NumObjects = sum(S5P_comps_large_inds);

            NumClusters = S5P_comps_large.NumObjects;

            LL_small = labelmatrix(S5P_comps_small);
            LL_large = labelmatrix(S5P_comps_large);
            maskImage = zeros(size(LL_small));
            maskImage(LL_small>0) = 1;
            maskImage(LL_large>0) = 2;

            S5P_comps = S5P_comps_large;

            S5P_props = regionprops3(S5P_comps,S5P_subImage,...
                'Volume','Solidity',...
                'Centroid','Image','BoundingBox');

            S5P_Volume_array = ...
                [S5P_props.Volume].*pixelSize.^2.*zStepSize;
            S5P_Solidity_array = [S5P_props.Solidity];
            S5P_Centroid_array = S5P_props.Centroid;
            S5P_Centroid_array(:,1) = S5P_Centroid_array(:,1)-1;
            S5P_Centroid_array(:,2) = S5P_Centroid_array(:,2)-1;
            S5P_Centroid_array = S5P_Centroid_array...
                .*[pixelSize,pixelSize,zStepSize];

            Volumes = S5P_Volume_array'
            Solidities = S5P_Solidity_array'

            clusterInd = exampleClusterInds(ff);
            S5P_Volume = S5P_Volume_array(clusterInd);
            S5P_Solidity = S5P_Solidity_array(clusterInd);
            S5P_Centroid = S5P_Centroid_array(clusterInd,:);
            zPlaneInd = round(S5P_props.Centroid(clusterInd,3));

            % --- get cluster central plane and elongation in-plane

            % --- analyze and plot target cluster
            object_nn = clusterInd;

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

            S5P_Slices = cell(1,numStoreChannels+2);
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
                S5P_Slices{color} = ...
                    squeeze(store_subImages{color}(...
                    img_limits(1):img_limits(2),...
                    img_limits(3):img_limits(4),...
                    center_z));
            end
            
            S5P_Slices{numStoreChannels+1} = ...
                squeeze(maskImage(...
                img_limits(1):img_limits(2),...
                img_limits(3):img_limits(4),...
                center_z));
            
            % --- quantification for all target channels
            for qq = 1:numQuantChannels

                channelInd = quantChannels(qq);
                quant_subImage = imgStack{channelInd}(...
                    boxArray(2)-20+0.5:boxArray(2)+boxArray(5)+20-0.5,...
                    boxArray(1)-20+0.5:boxArray(1)+boxArray(4)+20-0.5,...
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
                S5P_intensity{qq} = ...
                    Quant_ClusterMedian./Quant_nucleusMedian;
                S5P_nucIntensity{qq} = nuc_intCell{qq}(nn);

            end
        
        end

        subplot(6,numExamples,ff)

        x_min = (centroid_2-ImgSquareExtension).*pixelSize;
        x_max = (centroid_2+ImgSquareExtension).*pixelSize;
        y_min = (centroid_1-ImgSquareExtension).*pixelSize;
        y_max = (centroid_1+ImgSquareExtension).*pixelSize;

        boxSize = size(S5P_subImage_raw);
        centerPlaneInd = center_z;
        imagesc([1,boxSize(2)].*pixelSize,...
            [1,boxSize(1)].*pixelSize,...
            squeeze(max(S5P_subImage_raw(:,:,...
            (zPlaneInd-1):(zPlaneInd+1)),[],3)),...
            S5P_lims)
        axis tight equal
        set(gca,'Colormap',inferno,'XTick',[],'YTick',[])
        hold on
        plot([0.5,4.5],[0.5,0.5],'w-','LineWidth',4)
        text(2.5,1.5,'4 \mum','Color',[1,1,1])

        plot([x_min,x_min,x_max,x_max,x_min],...
            [y_max,y_min,y_min,y_max,y_max],...
            'w-','LineWidth',1)

        title(thisCondName,'FontWeight','normal')

        hold off


        subplot(6,numExamples,ff+1.*numExamples)
        imagesc([1,boxSize(2)].*pixelSize,...
            [1,boxSize(1)].*pixelSize,...
            squeeze(max(S2P_subImage_raw(:,:,...
            (zPlaneInd-1):(zPlaneInd+1)),[],3)),...
            S2P_lims)
        axis tight equal
        set(gca,'Colormap',inferno,'XTick',[],'YTick',[])
        hold on
        plot([0.5,4.5],[0.5,0.5],'w-','LineWidth',4)

        plot([x_min,x_min,x_max,x_max,x_min],...
            [y_max,y_min,y_min,y_max,y_max],...
            'w-','LineWidth',1)

        %title(thisCondName,'FontWeight','normal')

        hold off


        subplot(6,numExamples,ff+2.*numExamples)
        imagesc([1,boxSize(2)].*pixelSize,...
            [1,boxSize(1)].*pixelSize,...
            squeeze(max(maskImage(:,:,...
            (zPlaneInd-1):(zPlaneInd+1)),[],3)))
        axis tight equal
        set(gca,'Colormap',...
            [1,1,1;0.7,0.7,0.7;0.6,0,0],'XTick',[],'YTick',[])
        hold on
        plot([0.5,4.5],[0.5,0.5],'k-','LineWidth',4)

        plot([x_min,x_min,x_max,x_max,x_min],...
            [y_max,y_min,y_min,y_max,y_max],...
            'k-','LineWidth',1)

        hold off

        title(sprintf('%d large clusters',...
            NumClusters),...
            'Color',[1,0,0],'FontWeight','normal')

        subplot(6,numExamples,ff+3.*numExamples)
        
        interpolFac = 3;

        smallS5PImg = imresize(S5P_Slices{1},interpolFac);
        smallS2PImg = imresize(S5P_Slices{2},interpolFac);
        smallMaskImg = imresize(S5P_Slices{3}==2,interpolFac);

        smallPerimImg = bwperim(smallMaskImg);
        
        smallS5PImg_outline = smallS5PImg;
        smallS2PImg_outline = smallS2PImg;

        smallS5PImg_outline(smallPerimImg) = S5P_zoom_lims(2);
        smallS2PImg_outline(smallPerimImg) = S2P_zoom_lims(2);

        small_boxSize = size(smallS5PImg);
        imagesc([1,small_boxSize(2)].*pixelSize./interpolFac,...
            [1,small_boxSize(1)].*pixelSize./interpolFac,...
            smallS5PImg_outline,...
            S5P_zoom_lims)
        axis tight equal
        set(gca,'Colormap',inferno,'XTick',[],'YTick',[])
        hold on
        plot([0.1,0.6],[0.1,0.1],'w-','LineWidth',4)
        text(0.35,0.3,'500 nm','Color',[1,1,1])
        set(gca,'XTick',[],'YTick',[])

        hold off

        title(sprintf('I_{S5P}=%2.2f',...
            S5P_intensity{1}(object_nn)),...
            'Color',[1,0,0],'FontWeight','normal')

        subplot(6,numExamples,ff+4.*numExamples)
        
        imagesc([1,small_boxSize(2)].*pixelSize./interpolFac,...
            [1,small_boxSize(1)].*pixelSize./interpolFac,...
            smallS2PImg_outline,...
            S2P_zoom_lims)
        axis tight equal
        set(gca,'Colormap',inferno,'XTick',[],'YTick',[])
        hold on
        plot([0.1,0.6],[0.1,0.1],'w-','LineWidth',4)
        set(gca,'XTick',[])

        hold off

        title(sprintf('I_{S2P}=%2.2f',...
            S5P_intensity{2}(object_nn)),...
            'Color',[0.4,0.4,0.4],'FontWeight','normal')

        
        subplot(6,numExamples,ff+5.*numExamples)

        this_magenta_plot = ...
            (smallS5PImg-S5P_zoom_lims(1)) ...
            ./(S5P_zoom_lims(2)-S5P_zoom_lims(1));
        this_green_plot = ...
            (smallS2PImg-S2P_zoom_lims(1)) ...
            ./(S2P_zoom_lims(2)-S2P_zoom_lims(1));
        thisSize = size(this_magenta_plot);

%         this_magenta_plot = this_magenta_plot-min(this_magenta_plot(:));
%         this_magenta_plot = this_magenta_plot./max(this_magenta_plot(:));
%         this_green_plot = this_green_plot-min(this_green_plot(:));
%         this_green_plot = this_green_plot./max(this_green_plot(:));


        
    	redChannel = this_magenta_plot;
    	blueChannel = this_magenta_plot;
    	greenChannel = this_green_plot;
    	rgb_img = zeros(thisSize(1),thisSize(2),3);
    	rgb_img(:,:,1) = redChannel;
    	rgb_img(:,:,2) = greenChannel;
    	rgb_img(:,:,3) = blueChannel;

    	image([1,thisSize(2)].*pixelSize./interpolFac,...
    		[1,thisSize(1)].*pixelSize./interpolFac,...
    		rgb_img)

        axis tight equal

        hold on
        plot([0.1,0.6],[0.1,0.1],'w-','LineWidth',4)
        set(gca,'XTick',[],'YTick',[])

        hold off

        title(sprintf('V=%2.2f \\mum^3, S=%2.2f',...
            S5P_Volume_array(object_nn),S5P_Solidity_array(object_nn)),...
            'Color',[0,0,0],'FontWeight','normal')


%         subplot(6,6,22)
% 
%         small_boxSize = size(S5P_Slices{3});
%         imagesc([1,small_boxSize(2)].*pixelSize,...
%             [1,small_boxSize(1)].*pixelSize,...
%             S5P_Slices{3})
%         axis tight equal
%         set(gca,'Colormap',...
%             [1,1,1;0.7,0.7,0.7;0.6,0,0])        
%         hold on
%         plot([0.1,0.6],[0.1,0.1],'k-','LineWidth',6)
%         set(gca,'XTick',[])
% 
%         hold off
% 
%         title(sprintf('Sol=%2.2f',...
%             S5P_Solidity_array(object_nn)),...
%             'Color',[0,0,0],'FontWeight','normal')

        

    end
end