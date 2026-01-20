%% Analyze RNA polymerase II cluster features
%
% This script performs deeper analyses on previously extracted RNA polymerase II
% (Pol II) clusters. Two types of Pol II clusters are differentiated:
% - Pol II S5P clusters (Serine 5-phosphorylated Pol II, "recruited Pol II") and
% - Pol II S2P clusters (Serine 2-phosphorylated Pol II, "elongating Pol II").
%
% This script requires a result data file containing the extracted (and
% potentially combined, after multi-batch extraction) cluster features, as well
% as the feature table result file that must have been created previously from
% the same result data file.

clear all

%% Process parameter section

% import data file path
sourceFile_dat = fullfile(".", "NuclearClusterFeatures", "ClusterFeatures.mat");
sourceFile_tab = fullfile(".", "NuclearClusterFeatures", "ClusterFeatures_tables.mat");

%% Main script section

% ----- Create required result structures -----
%
% Note that these result structures are not already created like this in the
% extraction process on purpose. To save the data as created here in a MATLAB
% data file may require to use the file format version 7.3, as it may contain
% arrays larger than 2 GB. Unfortunately, the data structures created here
% contain massive amounts of nested cell arrays, which make the data file size
% in the file format version 7.3 explode. We are talking about a file size of
% ~60 GB, where the two source files that are imported here (and which contain
% all required data) have a combined size of ~1 GB. It is actually faster to do
% this initial data rearrangement in memory here than to import the huge data
% file that would be required to hold this data.

fprintf("Loading data ...\n");
tic

load(resfile1);
data = load(resfile2);

t = toc;
fprintf("Loading data finished (%.1f seconds)\n", t);

fprintf("Initializing new result structures ...\n");
tic

S2P_nucClustVolCell_old = S2P_nucClustVolCell;
S5P_nucClustVolCell_old = S5P_nucClustVolCell;
S2P_nucClustVolCell_new = cell(size(S2P_nucClustVolCell));
S5P_nucClustVolCell_new = cell(size(S5P_nucClustVolCell));
assert(size(S2P_nucClustVolCell_new, 2) == size(data.resTableDatasets, 1));
assert(size(S5P_nucClustVolCell_new, 2) == size(data.resTableDatasets, 1));

row = 1;
for d = 1:size(data.resTableDatasets, 1)
    S2P_nucClustVolCell_new{d} = cell(1, data.resTableDatasets.numNuclei(d));
    S5P_nucClustVolCell_new{d} = cell(1, data.resTableDatasets.numNuclei(d));
    for n = 1:data.resTableDatasets.numNuclei(d)
        S2P_nucClustVolCell_new{d}{n} = cell(data.resTableNuclei.numClustersS2P(row), 1);
        S5P_nucClustVolCell_new{d}{n} = cell(data.resTableNuclei.numClustersS5P(row), 1);
        row = row + 1;
    end
end

t = toc;
fprintf("Initializing new result structures finished (%.1f seconds)\n", t);

fprintf("Filling new result structures ...\n");
tic

dinds = find(validFileFlag);
for row = 1:size(data.resTableNuclei, 1)
    d = find(dinds == data.resTableNuclei.indDataset(row));
    n = data.resTableNuclei.indNucleus(row);
    assert(d <= length(S2P_nucClustVolCell_new));
    assert(n <= length(S2P_nucClustVolCell_new{d}));
    if n == 1
        i1s2p = 1;
        i1s5p = 1;
    else
        i1s2p = i2s2p + 1;
        i1s5p = i2s5p + 1;
    end
    i2s2p = i1s2p + data.resTableNuclei.numClustersS2P(row) - 1;
    i2s5p = i1s5p + data.resTableNuclei.numClustersS5P(row) - 1;
    % fprintf("(%5d / %5d) d = %4d, n = %3d, S2P clusters: %5d (%5d : %5d), S2P clusters: %5d (%5d : %5d)\n", row, size(data.resTableNuclei, 1), d, n, data.resTableNuclei.numClustersS2P(row), i1s2p, i2s2p, data.resTableNuclei.numClustersS5P(row), i1s5p, i2s5p);
    assert(i2s2p <= length(S2P_nucClustVolCell_old{d}));
    assert(i2s5p <= length(S5P_nucClustVolCell_old{d}));
    for c = 1:length(i1s2p:i2s2p)
        S2P_nucClustVolCell_new{d}{n}{c} = S2P_nucClustVolCell_old{d}(i1s2p:i2s2p);
    end
    for c = 1:length(i1s5p:i2s5p)
        S5P_nucClustVolCell_new{d}{n}{c} = S5P_nucClustVolCell_old{d}(i1s5p:i2s5p);
    end
end
S2P_nucClustVolCell = S2P_nucClustVolCell_new;
S5P_nucClustVolCell = S5P_nucClustVolCell_new;

t = toc;
fprintf("Filling new result structures finished (%.1f seconds)\n", t);

clear data S2P_nucClustVolCell_old S2P_nucClustVolCell_new S5P_nucClustVolCell_old S5P_nucClustVolCell_new d dinds row n i1s2p i1s5p i2s2p i2s5p c



% ----- Sort into conditions -----

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
