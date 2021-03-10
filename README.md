# ImageAnalysisMyxo

This repository contains custom scripts for the analyses presented in the paper "PomX, a ParA/MinD ATPase activating protein, is a triple regulator of cell division in *Myxococcus xanthus*" by Dominik Schumacher, Andrea Harms, Silke Bergeler, Erwin Frey and Lotte Søgaard-Andersen (https://www.biorxiv.org/content/10.1101/2020.12.14.422651v1). 

The scripts automatically analyze *Myxococcus xanthus* cells and fluorescently labelled protein foci within the cells. It builds on the software Oufti (http://www.oufti.org) from the Jacobs-Wagner lab. Note that the scripts were written specifically for the analyses in the paper above. The scripts have only been tested on the specific images used in the above publication on a Mac with Matlab 2020b and on Windows with MatLab 2020b. You may use and modify the scripts at your own risk. 

# How to run the scripts

There are two main scripts, one for the image analysis (main_imageAnalysis.m) and a second one for the visualization of the results (main_visualization.m). First, run the main_imageAnalysis.m script to automatically obtain measurements of the cells and, if applicable, fluorescent spots within the cells. The results are stored together with the data. Then, run the main_visualization.m script to obtain customized plots of the results and create Excel files with the data. 

# Image analysis 

## Input
The following input is needed to run the image analysis script:
* .mat file from Oufti containing the results from cell segmentation (the meshes for the cells) 
* Phase contrast images of cells (same ones used in Oufti to detect cell outlines)
* Corresponding fluorescence images of the cells (if applicable)

## Types of analyses
The script allows for two different kind of analyses:
1. Analysis of snapshots
2. Analysis of time lapse images
In the latter case, the analysis distinguishes between cells that divide and those that don't divide. For cells that do not divide, the analysis can be performed with cell outlines (obtained from Oufti) for the first frame only. The position of the cells in the first frame is then used to detect clusters within the cells and these clusters are then tracked for subsequent frames. This allows for a quantification of the clusters even if cell segmentation becomes difficult as in the case of very long cells. 

## Cell tracking (for time lapse images)
We modified scripts from Oufti to optimize cell tracking for our data. Here, association of cells in consecutive frames are obtained via the Kuhn-Munkres algorithm. 

## Spot detection
Spots of fluorescently-labelled proteins within cells are detected based on a threshold intensity that is defined based on the diffuse fluorescence signal in the cell and a minimal number of connected high-intensity pixels.

## Spot tracking (for cells that do not divide)
For cells that do not divide, only the first frame needs to have cell outlines (from Oufti). Spots are then detected within the cells of the first frame. For subsequent frames, spots are detected and tracked using the Kuhn-Munkres algorithm. 

Note: this analysis is still work in progress. 

# Visualization
The analysis results are visualized and saved in tables that are suitable for further analysis and visualization e.g. in Excel. 




