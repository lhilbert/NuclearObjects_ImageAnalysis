# NuclearObjects_ImageAnalysis

# MatLab scripts for image analysis of objects inside cell nuclei

Currently, the repository contains the scripts for extracting and analyzing Nikon .nd2 files with multiple XY positions in them.

# Requirements

You will need MatLab installed on your computer. You will also need the Image Processing toolbox. The Parallel toolbox is helpful, as it can speed up your image processing.

You also first have to download bfmatlab, and add it to your MatLab path.

bfmatlab:
https://www.openmicroscopy.org/bio-formats/downloads/

Adding a directory to the MatLab path permanently:
https://www.mathworks.com/matlabcentral/answers/97215-how-do-i-automatically-add-folders-to-the-matlab-path-on-startup

# Use instructions

1. Specifiy which files shoudl be extracted, and what color channels shoudl be retained in MultiPosition_extraction_nd2.m, then run the script to oextract the raw data into .mat MatLab format data files
2. Specify analysis parameters in ClusterAnalysis.m, then run the script to actually analyze your data
3. Once you see that this is working, you can specify the cootents of the ClusterAnalysis.m file t fit the needs of your analysis and graphical outputs
4. You can also try out ExampleImages.m, so you can prepare nice plots based on the extracted data

# Acknowledgments

The scripts make use of the rdir.m scripts for recursive directory search

https://www.mathworks.com/matlabcentral/fileexchange/47125-rdir-m
