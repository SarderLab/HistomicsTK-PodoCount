# PodoCount in the Cloud 
__A cloud-based tool for whole-slide podocyte quantification__

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
  
An instructional video has been prepared as a guide for first-time users and is available at http://bit.ly/3rdGPEd 


### For questions or feedback, please contact:
- Briana Santo <basanto@buffalo.edu>
- Pinaki Sarder <pinakisa@buffalo.edu>


