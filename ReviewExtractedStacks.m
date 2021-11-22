clear all

% Specify the directory that contains the extraced files from the last
% step, where you extracted from the raw files obtained from the microscope
sourceDirectory = './ExtractedStacks/**/';
% sourceDirectory = './ExtractedStacks/Cond_15/';

% --- file read and display procedure begins here

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

for ff = 1:numFiles
	
	fprintf('Processing file %d of %d\n',ff,numFiles)
	
	thisCondInd = condInds(ff);	
    thisCondName = condNames{ff};
	thisFilePath = listing(ff).name;
	
	loadStruct = load(thisFilePath,...
		'imgStack','imgSize','pixelSize','zStepSize');
	imgStack = loadStruct.imgStack;
	imgSize = loadStruct.imgSize;
	pixelSize = loadStruct.pixelSize;
	zStepSize = loadStruct.zStepSize;
    centralSlice = round(imgSize(3)./2);
    numChannels = numel(imgStack);
	
    for nn = 1:numChannels
        
        subplot(1,numChannels,nn)
        imagesc([0,imgSize(2)].*pixelSize,...
            [0,imgSize(1)].*pixelSize,...
            squeeze(imgStack{nn}(:,:,centralSlice)))
        if nn == 2
            title(sprintf('File name: %s',thisFilePath),...
                'interpreter','none')
        elseif nn == 1
            title(sprintf('Condition name: %s',thisCondName),...
                'interpreter','none')
        end
        axis tight equal
        
    end
    
    fprintf('File name: %s\n',thisFilePath)
    waitforbuttonpress
    
end
