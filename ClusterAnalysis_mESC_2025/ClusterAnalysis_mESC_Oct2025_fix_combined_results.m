resfile1 = fullfile(".", "AfterObjectAnalysis_combined_v7.mat");
resfile2 = fullfile(".", "AfterObjectAnalysis_combined_tab_v7.mat");
% save_file_2 = fullfile(".", "AfterObjectAnalysis_combined_2_v7.mat");
save_file_3 = fullfile(".", "AfterObjectAnalysis_combined_3_v73.mat");
load(resfile1);
data = load(resfile2);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for save_file_3 only
    for n = 1:data.resTableDatasets.numNuclei(d)
        S2P_nucClustVolCell_new{d}{n} = cell(data.resTableNuclei.numClustersS2P(row), 1);
        S5P_nucClustVolCell_new{d}{n} = cell(data.resTableNuclei.numClustersS5P(row), 1);
        row = row + 1;
    end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

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
    fprintf("(%5d) d = %4d, n = %3d, S2P clusters: %5d (%5d : %5d), S2P clusters: %5d (%5d : %5d)\n", row, d, n, data.resTableNuclei.numClustersS2P(row), i1s2p, i2s2p, data.resTableNuclei.numClustersS5P(row), i1s5p, i2s5p);
    assert(i2s2p <= length(S2P_nucClustVolCell_old{d}));
    assert(i2s5p <= length(S5P_nucClustVolCell_old{d}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for save_file_2 only
    % S2P_nucClustVolCell_new{d}{n} = S2P_nucClustVolCell_old{d}(i1s2p:i2s2p);
    % S5P_nucClustVolCell_new{d}{n} = S5P_nucClustVolCell_old{d}(i1s5p:i2s5p);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for save_file_3 only
    for c = 1:length(i1s2p:i2s2p)
        S2P_nucClustVolCell_new{d}{n}{c} = S2P_nucClustVolCell_old{d}(i1s2p:i2s2p);
    end
    for c = 1:length(i1s5p:i2s5p)
        S5P_nucClustVolCell_new{d}{n}{c} = S5P_nucClustVolCell_old{d}(i1s5p:i2s5p);
    end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
S2P_nucClustVolCell = S2P_nucClustVolCell_new;
S5P_nucClustVolCell = S5P_nucClustVolCell_new;

clear resfile1 resfile2 data S2P_nucClustVolCell_old S2P_nucClustVolCell_new S5P_nucClustVolCell_old S5P_nucClustVolCell_new d dinds row n i1s2p i1s5p i2s2p i2s5p c cyto_intCell

% save(save_file_2)
save(save_file_3, "-v7.3")
