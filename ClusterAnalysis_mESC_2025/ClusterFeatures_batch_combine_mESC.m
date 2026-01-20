%% Combine RNA polymerase II cluster features from multi-batch extraction
%
% This script combines imports RNA polymerase II cluster features from a set of
% result files from multi-batch-based cluster extraction and combines them into
% a single result file.
%
% This step is required for all further analysis steps, because they work on
% combined results in a single result file.
%
% This step can be omitted if the cluster feature extraction has been performed
% in a single-batch run.

clear all

%% Process parameter section

% data containing the batch processing results from the cluster analysis step
sourceFiles = fullfile(".", "NuclearClusterFeatures", ...
    [ ...
    "ClusterFeatures_batch01.mat", ...
    "ClusterFeatures_batch02.mat", ...
    "ClusterFeatures_batch03.mat", ...
    "ClusterFeatures_batch04.mat", ...
    "ClusterFeatures_batch05.mat", ...
    "ClusterFeatures_batch06.mat", ...
    "ClusterFeatures_batch07.mat", ...
    "ClusterFeatures_batch08.mat", ...
    ]);

% combined cluster data output file
save_file = fullfile(".", "NuclearClusterFeatures", "ClusterFeatures.mat");

%% Main script section

% Load result files
data = cell(size(sourceFiles));
for ff = 1:length(sourceFiles)
    data{ff} = load(fullfile(sourceDirectory, sourceFiles(ff)));
end

% Combine process parameters
plotFlag = data{1}.plotFlag;
sourceDirectory = data{1}.sourceDirectory;
NucSegChannel = data{1}.NucSegChannel;
S5P_SegChannel = data{1}.S5P_SegChannel;
S2P_SegChannel = data{1}.S2P_SegChannel;
ImgSquareExtension = data{1}.ImgSquareExtension;
storeImgChannels = data{1}.storeImgChannels;
quantChannels = data{1}.quantChannels;
quantBlurSigma = data{1}.quantBlurSigma;
nuc_segBlurSigma_nucleus = data{1}.nuc_segBlurSigma_nucleus;
nuc_segBlurSigma_BG_removal = data{1}.nuc_segBlurSigma_BG_removal;
nuc_segErosion = data{1}.nuc_segErosion;
fillHolesFlag = data{1}.fillHolesFlag;
Nuc_min_vol = data{1}.Nuc_min_vol;
Nuc_min_sol = data{1}.Nuc_min_sol;
Nuc_min_CoV = data{1}.Nuc_min_CoV;
cytoMask_extension = data{1}.cytoMask_extension;
cytoMask_distance = data{1}.cytoMask_distance;
S5P_segBlurSigma_object = data{1}.S5P_segBlurSigma_object;
S5P_segBlurSigma_BG_removal = data{1}.S5P_segBlurSigma_BG_removal;
S5P_seg_numStdDev = data{1}.S5P_seg_numStdDev;
S5P_DBSCAN_epsilon = data{1}.S5P_DBSCAN_epsilon;
S2P_segBlurSigma_object = data{1}.S2P_segBlurSigma_object;
S2P_segBlurSigma_BG_removal = data{1}.S2P_segBlurSigma_BG_removal;
S2P_seg_numStdDev = data{1}.S2P_seg_numStdDev;
S5P_minVol = data{1}.S5P_minVol;
S2P_minVol = data{1}.S2P_minVol;
numStoreChannels = data{1}.numStoreChannels;
numQuantChannels = data{1}.numQuantChannels;
listing = data{1}.listing;
numFiles = data{1}.numFiles;
for ff = 2:length(data)
    assert(plotFlag == data{ff}.plotFlag);
    % assert(sourceDirectory == data{ff}.sourceDirectory);
    assert(NucSegChannel == data{ff}.NucSegChannel);
    assert(S5P_SegChannel == data{ff}.S5P_SegChannel);
    assert(S2P_SegChannel == data{ff}.S2P_SegChannel);
    assert(ImgSquareExtension == data{ff}.ImgSquareExtension);
    assert(all(storeImgChannels == data{ff}.storeImgChannels));
    assert(all(quantChannels == data{ff}.quantChannels));
    assert(all(quantBlurSigma == data{ff}.quantBlurSigma));
    assert(nuc_segBlurSigma_nucleus == data{ff}.nuc_segBlurSigma_nucleus);
    assert(nuc_segBlurSigma_BG_removal == data{ff}.nuc_segBlurSigma_BG_removal);
    assert(nuc_segErosion == data{ff}.nuc_segErosion);
    assert(fillHolesFlag == data{ff}.fillHolesFlag);
    assert(Nuc_min_vol == data{ff}.Nuc_min_vol);
    assert(Nuc_min_sol == data{ff}.Nuc_min_sol);
    assert(Nuc_min_CoV == data{ff}.Nuc_min_CoV);
    assert(cytoMask_extension == data{ff}.cytoMask_extension);
    assert(cytoMask_distance == data{ff}.cytoMask_distance);
    assert(S5P_segBlurSigma_object == data{ff}.S5P_segBlurSigma_object);
    assert(S5P_segBlurSigma_BG_removal == data{ff}.S5P_segBlurSigma_BG_removal);
    assert(S5P_seg_numStdDev == data{ff}.S5P_seg_numStdDev);
    assert(S5P_DBSCAN_epsilon == data{ff}.S5P_DBSCAN_epsilon);
    assert(S2P_segBlurSigma_object == data{ff}.S2P_segBlurSigma_object);
    assert(S2P_segBlurSigma_BG_removal == data{ff}.S2P_segBlurSigma_BG_removal);
    assert(S2P_seg_numStdDev == data{ff}.S2P_seg_numStdDev);
    assert(S5P_minVol == data{ff}.S5P_minVol);
    assert(S2P_minVol == data{ff}.S2P_minVol);
    assert(numStoreChannels == data{ff}.numStoreChannels);
    assert(numQuantChannels == data{ff}.numQuantChannels);
    % assert(listing == data{ff}.listing);
    assert(numFiles == data{ff}.numFiles);
end

% Combine Result variables
validFileFlag = data{1}.validFileFlag;
condInds = data{1}.condInds;
condNames = data{1}.condNames;
numNuclei_vec = data{1}.numNuclei_vec;
nuc_volCell = data{1}.nuc_volCell;
nuc_intCell = data{1}.nuc_intCell;
cyto_intCell = data{1}.cyto_intCell;
nuc_stdCell = data{1}.nuc_stdCell;
nuc_medianVolCell = data{1}.nuc_medianVolCell;
perNuc_countCell = data{1}.perNuc_countCell;
perNuc_volCell = data{1}.perNuc_volCell;
S5P_xyVoxelSizeCell = data{1}.S5P_xyVoxelSizeCell;
S5P_zVoxelSizeCell = data{1}.S5P_zVoxelSizeCell;
S2P_xyVoxelSizeCell = data{1}.S2P_xyVoxelSizeCell;
S2P_zVoxelSizeCell = data{1}.S2P_zVoxelSizeCell;
S5P_volCell = data{1}.S5P_volCell;
S5P_solCell = data{1}.S5P_solCell;
S5P_eloCell = data{1}.S5P_eloCell;
S5P_intCell = data{1}.S5P_intCell;
S5P_centCell = data{1}.S5P_centCell;
S5P_imgCell = data{1}.S5P_imgCell;
S5P_nucIntCell = data{1}.S5P_nucIntCell;
S5P_nucVolCell = data{1}.S5P_nucVolCell;
S5P_nucClustVolCell = data{1}.S5P_nucClustVolCell;
S2P_volCell = data{1}.S2P_volCell;
S2P_solCell = data{1}.S2P_solCell;
S2P_eloCell = data{1}.S2P_eloCell;
S2P_intCell = data{1}.S2P_intCell;
S2P_centCell = data{1}.S2P_centCell;
S2P_imgCell = data{1}.S2P_imgCell;
S2P_nucIntCell = data{1}.S2P_nucIntCell;
S2P_nucVolCell = data{1}.S2P_nucVolCell;
S2P_nucClustVolCell = data{1}.S2P_nucClustVolCell;
for ff = 2:length(data) %#ok<*AGROW>
    validFileFlag = validFileFlag | data{ff}.validFileFlag;
    condInds = [condInds, data{ff}.condInds];
    condNames = [condNames, data{ff}.condNames];
    numNuclei_vec = [numNuclei_vec, data{ff}.numNuclei_vec];
    nuc_volCell = [nuc_volCell, data{ff}.nuc_volCell];
    nuc_intCell = [nuc_intCell, data{ff}.nuc_intCell];
    cyto_intCell = [cyto_intCell, data{ff}.cyto_intCell];
    nuc_stdCell = [nuc_stdCell, data{ff}.nuc_stdCell];
    nuc_medianVolCell = [nuc_medianVolCell, data{ff}.nuc_medianVolCell];
    perNuc_countCell = [perNuc_countCell, data{ff}.perNuc_countCell];
    perNuc_volCell = [perNuc_volCell, data{ff}.perNuc_volCell];
    S5P_xyVoxelSizeCell = [S5P_xyVoxelSizeCell, data{ff}.S5P_xyVoxelSizeCell];
    S5P_zVoxelSizeCell = [S5P_zVoxelSizeCell, data{ff}.S5P_zVoxelSizeCell];
    S2P_xyVoxelSizeCell = [S2P_xyVoxelSizeCell, data{ff}.S2P_xyVoxelSizeCell];
    S2P_zVoxelSizeCell = [S2P_zVoxelSizeCell, data{ff}.S2P_zVoxelSizeCell];
    S5P_volCell = [S5P_volCell, data{ff}.S5P_volCell];
    S5P_solCell = [S5P_solCell, data{ff}.S5P_solCell];
    S5P_eloCell = [S5P_eloCell, data{ff}.S5P_eloCell];
    S5P_intCell = [S5P_intCell, data{ff}.S5P_intCell];
    S5P_centCell = [S5P_centCell, data{ff}.S5P_centCell];
    S5P_imgCell = [S5P_imgCell, data{ff}.S5P_imgCell];
    S5P_nucIntCell = [S5P_nucIntCell, data{ff}.S5P_nucIntCell];
    S5P_nucVolCell = [S5P_nucVolCell, data{ff}.S5P_nucVolCell];
    S5P_nucClustVolCell = [S5P_nucClustVolCell, data{ff}.S5P_nucClustVolCell];
    S2P_volCell = [S2P_volCell, data{ff}.S2P_volCell];
    S2P_solCell = [S2P_solCell, data{ff}.S2P_solCell];
    S2P_eloCell = [S2P_eloCell, data{ff}.S2P_eloCell];
    S2P_intCell = [S2P_intCell, data{ff}.S2P_intCell];
    S2P_centCell = [S2P_centCell, data{ff}.S2P_centCell];
    S2P_imgCell = [S2P_imgCell, data{ff}.S2P_imgCell];
    S2P_nucIntCell = [S2P_nucIntCell, data{ff}.S2P_nucIntCell];
    S2P_nucVolCell = [S2P_nucVolCell, data{ff}.S2P_nucVolCell];
    S2P_nucClustVolCell = [S2P_nucClustVolCell, data{ff}.S2P_nucClustVolCell];
end

% Save combined results

% make sure not to use the v7.3 file format, it would become GIGANTIC because of
% the nested cell structures

clear("data", "ff", "sourceFiles");

save(save_file)
