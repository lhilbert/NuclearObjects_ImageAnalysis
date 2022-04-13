# NuclearObjects_ImageAnalysis

MatLab scripts for image analysis of objects inside cell nuclei

Currently, the repository contains the scripts for extracting and analyzing Nikon .nd2 files with multiple XY positions in them.

# Requirements

You will need MatLab installed on your computer. You will also need the **Image Processing toolbox** and the **Statistic and Machine Learning toolbox**. The **Parallel toolbox** is helpful, as it can speed up your image processing.

You also first have to download bfmatlab, and add it to your MatLab path.

bfmatlab:
https://www.openmicroscopy.org/bio-formats/downloads/

Adding a directory to the MatLab path permanently:
1. Click the option "Set Path" in the Tab "Home"
2. Navigate to the place where you save bfmatlab, open the directory and add to path
3. Click the button "Preferences" in the Tab "Home"
4. Select the bfmatlab directory, then click the button "Save"

# Use instructions

1. Specifiy which files should be extracted, and what color channels should be retained in MultiPosition_extraction_nd2.m, then run the script to oextract the raw data into .mat MatLab format data files
2. Specify analysis parameters in ClusterAnalysis.m, then run the script to actually analyze your data
3. Once you see that this is working, you can specify the cootents of the ClusterAnalysis.m file t fit the needs of your analysis and graphical outputs
4. You can also try out ExampleImages.m, so you can prepare nice plots based on the extracted data

# Example data set

A useful data set for an example analysis is provided by the Hilbert lab via Zenodo:
https://zenodo.org/record/5242952#.YZkwLy1Q1pQ

First, create a folder called ImageData, which you place in the directory where you have saved the MatLab scripts. It is the easiest way to keep the scripts and data together in one and the same place.

Second, inside the ImageData folder, create two more folders, called Ctrl and FP.

In the Ctrl folder, store the following two files:

https://zenodo.org/record/5242952/files/SetC_Control_004_crop.nd2?download=1
https://zenodo.org/record/5242952/files/SetD_Control_007_crop.nd2?download=1

In the FP folder, store the following two files:

https://zenodo.org/record/5242952/files/SetC_Flavopiridol_005_crop.nd2?download=1
https://zenodo.org/record/5242952/files/SetD_Flavopiridol_008_crop.nd2?download=1

# Acknowledgments

The scripts make use of the rdir.m scripts for recursive directory search

https://www.mathworks.com/matlabcentral/fileexchange/47125-rdir-m
