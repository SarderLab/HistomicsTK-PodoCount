# PodoCount in the Cloud: A cloud-based tool for whole-slide podocyte quantification

**Version 1.0.0**

This repository contains the source codes for the publication, "PodoCount: A robust, fully auotmated whole-slide podocyte quantification tool." All algorithms were developed and written by Briana Santo. The function xmltomask has been adapted from *Lutnick* et al.'s work "An integrated iterative annotation technique for easing neural network training in medical image analysis," *Nature Machine Intelligence*, 2019. The function "TranslateXMLtoJSON" has been adapted from *Shashiprakash* et al.'s work "A distributed system improves inter-observer and AI concordance in annotating interstitial fibrosis and tubular atrophy" [2020] as well as *Lutnick* et al.'s work "User friendly, cloud based, whole slide image segmentation" [2021].

---
## A pre-print version of this work is available at: X

*Prepared by Briana Santo at SUNY Buffalo on 27July2021*

## Image Data

Whole slide images (WSIs) of murine kidney data are available at: http://bit.ly/3rdGPEd. Murine data includes whole kidney sections derived from both wild type and diseased mice across six mouse models [T2DM A, T2DM B, Aging, FSGS (SAND), HIVAN, Progeroid]. All kidney specimens were stained with p57kip2 immunohistochemistry and Periodic Acid-Schiff (without Hematoxylin counter stain) prior to digitization. 

## Requirements

This code runs using python3, and is enabled by HistomicsTK and Docker.

### Dependencies

- argparse [1.1]
- cv2 [3.2.0]
- lxml.etree [4.2.1]
- matplotlib [3.3.4]
- numpy [1.19.5]
- openslide-python [1.1.1]
- pandas [0.22.0]
- scikit-image [0.17.2]
- scipy [1.5.4]


## Usage: 
### Running PodoCount in the cloud via the Sarder Lab's Digital Slide Archive and Viewer

Access the Sarder Lab's Digital Slide Archive and Viewer at http://hermes.med.buffalo.edu:8080 

An account has been created and designated for end users who wish to experience PodoCount in the Cloud. The login credentials are:

	Username: experiencecloudpathology
	Password: sarderlab2021!
  
An instructional video has been prepared as a guide for first-time users and is available at http://bit.ly/3rdGPEd. 

Detailed instructions also follow.

### For questions or feedback, please contact:
- Briana Santo <basanto@buffalo.edu>
- Pinaki Sarder <pinakisa@buffalo.edu>

### Detailed Instructions

Accessing and logging on to the Sarder Lab's DSA:
- Navigate to http://hermes.med.buffalo.edu:8080
- Login
- In the upper right hand corner of the window select the drop down arrow to the right of your username
	- Navigate to {'Your username'} > My folders > Public
- Upload WSIs and corresponding glomerulus annotation files for analysis. 

To initiate PodoCount analysis:

- Click the line item corresponding to the name of the WSI that you would like to analyze. The WSI will open in ImageViewer. 
- In the upper right hand corner of the window select the drop down arrow to the right of "Actions". Navigate to Actions > Open in HistomicsUI. HistomicsUI will open as a new tab in your browser, and the selected WSI will be displayed. 
- To the left of your username are the Analyses options available. Using the drop down arrow navigate to Analyses > brianasanto/bas_feb2021_build9 > latest > PodoCount_human_analysis or PodoCount_mouse_analysis, depending on the kidney data that you're analyzing. 
- An analysis window will open on the left side of the screen. 
	Under the IO section, you will find the following subsections: 
-- Input Image -- should already be defined as the WSI selected for analysis.
-- Input Annotation File 1 -- you need to choose the corresponding glomerulus annotation file. Select the folder icon and navigate to {'Your Username'} > Home > Public > {'corresponding glomerulus annotation (.xml) file'}.
	   c.) Output Annotation File 1 -- you need to enter the desired name and location for the podocyte boundary annototation file that PodoCount will output. You will also need to indicate where you would like this file to be saved. 		We recommend navigating to {'Your Username'} > Home > Public, entering the desired file name, and clicking 'Save.'
   	   d.) Output Annotation File 2 -- you need to enter the desired name and location for the podocyte counter annototation file that PodoCount will output. We recommend doing as in item 4c. 
   	   e.) csvFile -- you need to enter the desired name and location for the podocyte and glomerulus feature file that PodoCount will output. We reccomend doing as in items 4c,d. 
	Under the UserParameters section, you will find the following subsection
   	   a.) Slider -- this parameter sets the threshold on the dab stain to isolate the immunohistochemically-positive podocyte nuclear regions. Optional values are in the range [0,3]. We recommend using the optimized value of 2.5 to 			start. You may adjust this value accordingly to best fit your image data. 
	Submit the job. 
5.) The UI will inform you when the job is complete. This step should take less than a minute. Once the job is complete, you may preview the podocyte boundary annotations from 4c by either downloading the .xml file to your local computer 	and opening the WSI is AperioImage Scope or converting the annotation file to a UI-friendly file format for web-based viewing. To convert the .xml to json format, use the drop down arrow to navigate to Analyses > 		brianasanto/bas_feb2021_build9 > latest > TranslateXMLToJson. 
6.) An analysis window will open on the left side of the screen. 
	Under the IO section, you will find the following subsections: 
   	   a.) Input Image -- should already be defined as the WSI selected for analysis.
   	   b.) Input Annotation File -- you need to choose the corresponding podocyte annotation (.xml) file output by PodoCount that you saved in 4c. Select the folder icon and navigate to the podocyte boundary annotation (.xml) file. 
	   c.) Output Annotation File -- you need to enter the desired name and location for the podocyte boundary annototation (json) file. 
	Submit the job. The dtected podocyte boundaries will automatically display in the UI when the job is complete. This step should take less than a minute. 



