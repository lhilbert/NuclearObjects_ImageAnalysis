function [ImgSeries_cell,arrayDim_cell,voxelSize_cell,seriesName_cell] = ...
   nd2read(combined_filepath,varargin)
% Retrieves intensity value arrays from .nd2 files (saved from Nikon
% Elements). As input, specific the comined filepath, including the path of
% the directory that the source file is stored in, and the indices of the
% images that you want to read in.
%
% The output is a cell of cells of 3D stacks, containing the intensity
% values from the different images that were requested. Following are a
% cell that contains the dimensions of the array, a cell that 
% contains the physical voxel sizes of the different arrays, and a cell
% that contains the names given to the different images before saving thee
% data on the microscope

reader = bfGetReader(combined_filepath);

% --- extract stack and microscope info from meta data

numSeries = reader.getSeriesCount;

if nargin > 1
	numSeries = min(numSeries,varargin{1});
end

ImgSeries_cell = cell(1,numSeries);
arrayDim_cell = cell(1,numSeries);
voxelSize_cell = cell(1,numSeries);
seriesName_cell = cell(1,numSeries);

actual_readImg_inds = 1:numSeries;

for nn = 1:numSeries
	
	actualSeries = actual_readImg_inds(nn)-1; % 0 based indexing!
	
	reader.setSeries(actualSeries);
	
	omeMeta = reader.getMetadataStore();
    
    % --- get the voxel edge sizes
    voxelSizeX = omeMeta.getPixelsPhysicalSizeX(actual_readImg_inds(nn)-1);
    voxelSizeX = voxelSizeX.value(ome.units.UNITS.MICROM);
    rawVoxelSizeX = voxelSizeX.doubleValue();
    voxelSizeY = omeMeta.getPixelsPhysicalSizeY(actual_readImg_inds(nn)-1);
    voxelSizeY = voxelSizeY.value(ome.units.UNITS.MICROM);
    rawVoxelSizeY = voxelSizeY.doubleValue();
	try
		voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(actual_readImg_inds(nn)-1);
		voxelSizeZ = voxelSizeZ.value(ome.units.UNITS.MICROM);
		rawVoxelSizeZ = voxelSizeZ.doubleValue();
	catch
		rawVoxelSizeZ = NaN;
	end
		
    voxelSize_cell{nn} = [rawVoxelSizeX,rawVoxelSizeY,rawVoxelSizeZ];
    
    % --- get the spatial stack dimensions
    rawStackSizeX = omeMeta.getPixelsSizeX(actual_readImg_inds(nn)-1).getValue(); % image width, pixels
    rawStackSizeY = omeMeta.getPixelsSizeY(actual_readImg_inds(nn)-1).getValue(); % image height, pixels
    rawStackSizeZ = omeMeta.getPixelsSizeZ(actual_readImg_inds(nn)-1).getValue(); % image height, pixels
	
	arrayDim_cell{nn} = [rawStackSizeY,rawStackSizeX,rawStackSizeZ];
	
	numChannels = reader.getSizeC();
	
	ImgSeries_cell{nn} = cell(1,numChannels);

	for cc = 1:numChannels
		
		ImgSeries_cell{nn}{cc} = ...
			zeros(rawStackSizeY,rawStackSizeX,rawStackSizeZ);
		
		for zz = 1:rawStackSizeZ
			
			planeInd = reader.getIndex(zz-1,cc-1,0)+1;
			planeImg = bfGetPlane(reader,planeInd);
			ImgSeries_cell{nn}{cc}(:,:,zz) = planeImg;
			
		end
		
	end
	
	ImgName = omeMeta.getImageName(actualSeries);
	seriesName_cell{nn} = ImgName.toCharArray';
	
end

reader.close()