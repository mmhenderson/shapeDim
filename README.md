Code related to our preprint: 

Henderson, M.M., Serences, J.T., & Rungratsameetaweemana, N. (2023). Dynamic categorization rules alter representations in human visual cortex.
[bioRxiv](https://doi.org/10.1101/2023.09.11.557257)

See our OSF repository for data files. 

### Instructions:
1. Clone this repository: ``` git clone https://github.com/mmhenderson/shapeDim ```
	* The location you clone to (including /shapeDim/) will become your ``` root ``` directory.
	
2. Download all files from the OSF repository, and unzip the zipped files.
	* After unzipping these files, put the folders ```DataBehavior``` and ```Samples``` inside your root directory, on same level as ```Analysis```
	* Put the folders ``` decoding_results``` and ```image_similarity``` into the ```Analysis``` folder.

3. Use the notebooks inside ```Analysis/final_analyses``` to reproduce our figures.
	* At the top of each notebook, be sure to change the variable ```root``` to your local path.

4. If you want to re-run our multivariate classifier analyses from scratch, use ```Analysis/multinomial_decoding/decode_multiclass.py```. 
	* You can also skip this (slow) step and just go straight to the jupyter notebooks - the saved decoding results are available as part of the OSF repo (```decoding_results_all.zip```)

### Contents:
* Analysis: Python code & jupyter notebooks needed to reproduce all our main analyses. 
* AnalyzeSilhouetteLoc: Matlab code providing a template for analyzing our localizer using FSL FEAT.
* ExpScriptsMRI: Matlab/Psychtoolbox code used to run our experiment.
* Preprocessing: Matlab code used for preprocessing our data.
* StimulusGeneration: Matlab code used for image generation.


Other notes:

* In parts of this code we refer to the "Nonlinear" task from our paper as the "Checker" task.
* This repo includes code for running and analyzing another task called the "Repeat" or "One-Back" task. We didn't include analyses of this task in our paper.
  
