# Time Cell Analysis project: Python and pybind11/C++ files for analysis and demos.


## Overview

This repository contains the Matlab, Python and related files for 
Time Cell Analysis,

This is from a forthcoming paper: 

Synthetic Data Resource and Benchmarks for Time Cell Analysis and Detection Algorithms
*K. Ananthamurthy and U.S. Bhalla,* 
in preparation.

The current directory, **TcPy**, has the Python demos and code (including 
pybind11/C++ code) for the tc module for Python. 

## Examples:
---------
For the impatient, once you have installed the files, here are some things
to do from the TcPy directory. In the examples below, the input data files
are Matlab version 7.3 files, organized as `[dataset][frame#][trial#][cell#]`
These files can be generated from data analysis code, or from the synthetic
data generation programs developed in this project.

Here are the key script calls:

> python ti_demo.py *filename*

This will run a partial analysis for time-cells in the specified file,
and report something like this:

```
Classification of first 30 cells for Dataset 0
CellIdx sigMean    SigTI   SigBoth pkFrame   FracTrialsFired
     0       0       0       0        0       0.6000
     1       0       0       0        3       0.3833
     2       1       1       1        7       0.7167
     3       1       1       1        7       0.7167
...
Number of time cells classified by each method, for each dataset
Dataset    #SigMean    #sigBoot     #sigBoth
   0          65          67          64
   1          64          67          63
   2          66          67          64
   3          44          17          13
...
```

Similar output is obtained from *r2b_demo.py*

A few additional stats are reported by *peq_demo.py*


If your source data file is generated from the synthetic data program, you 
can run:

> python ground_truth_check.py *filename*

and it will use ground-truth data from the file to see how well the various
methods handle the data. Output is like this:

```
    |Mau Peak Sig |   Mau TI    | Mau pk&&TI  |  r2b thresh | r2b bootstr | peq score
idx | tn fn fp tp | tn fn fp tp | tn fn fp tp | tn fn fp tp | tn fn fp tp | tn fn fp tp
  0 | 68  2  0 65 | 67  1  1 66 | 68  3  0 64 | 64  1  4 66 | 68  1  0 66 | 55 32 13 35
  1 | 67  5  1 62 | 68  0  0 67 | 68  5  0 62 | 65  0  3 67 | 68  1  0 66 | 44 35 24 32
  2 | 67  2  1 65 | 68  1  0 66 | 68  3  0 64 | 61  2  7 65 | 68  2  0 65 | 59 24  9 43
...
```

If you want to run a batch analysis on a synthetic data file do the following. 
It will dump all the data into .csv files. If a file exists with the same name
then the output is appended on to it.

> python run_batch_analysis.py *filename*

This will generate the files `ti.csv`, `r2b.csv`, `peq.csv` and `groundTruth.csv`

Here are three lines from a sample `ti.csv`:

```
0,0,0,0.4535,-0.1383,-0.1337,0,0,0.6000,0
0,1,0,0.1643,-0.0821,-0.0681,0,0,0.3833,3
0,2,1,2.6316,-0.1009,-0.1667,1,1,0.7167,7
```

The comma-separated entries are

	- datasetIdx:	Index of the dataset, i.e. recording session.
	- cellIdx:	Index of the cell whose stats are reported on this line
	- groundTruth:	Flag: Is the  cell a time cell?
	- isTimeCell:	Is the cell a time-cell? 1 if true, 0 if false.
	- meanScore: 	Peak of the mean of all trials for this cell
	- baseScore:	Temporal information for cell
	- percentileScore:	TI of cell vs TI of shuffled cells
	- sigMean:	Is the mean peak significantly different from shuffled?
	- sigBootstrap:	Is the TI significantly different from shuffled?
	- fracTrialsFired:	Hit Trial Ratio, fraction of trials where cell fired.
	- meanPkIdx:	Frame number of peak of mean cell activity. Note that this method bins frames (default by 3), and this is the binned frame number.

r2b.csv is very similar. r2b stands for Ridge-to-background, a metric for
the height of the peak vs the height of the shuffled peak.

	- datasetIdx:	Index of the dataset, i.e. recording session.
	- cellIdx:	Index of the cell whose stats are reported on this line
	- groundTruth:	Flag: Is the  cell a time cell?
	- meanScore: 	r2b average for shuffled trials.
	- baseScore:	r2b value for original trial.
	- percentileScore:	where does r2bscore place among shuffled r2bs?
	- sigMean:	Is the base score above a threshold? Usually 3.0.
	- sigBootstrap:	Is the percentile score above a threshold? Usually 99.5
	- fracTrialsFired:	Not computed, set to 0.
	- meanPkIdx:	Frame number of peak of mean cell activity. Note that this method does not do frame binning.

`peq.csv` is a bit different. peq stands for parametric equation quality, a metric for
estimating the quality of the recording in terms of key parameters of the dataset.
These parameters are noise, imprecision, event width, and hit trial ratio.
Here are the comma-separated entries for the peq.csv:

	- datasetIdx:	Index of the dataset, i.e. recording session.
	- cellIdx:	Index of the cell whose stats are reported on this line
	- groundTruth:	Flag: Is the  cell a time cell?
	- sigPeq:	Is this a time cell, ie, is the peq score > threshold?
	- baseScore:	The peq score, PEQ = HTR x exp -( alpha.N/S + beta*sdevEW/meanEW + gamma*sdevImp/stimWindow )
	- noise: 	Noise estimate for dataset.
	- eventWidth:	Mean event width for dataset
	- imprecision:	Imprecision of timing of the event.
	- fracTrialsFired:	Hit trial ratio, fraction of trials where cell fired.

Finally, `groundTruth.csv` reports how well the different methods compare to
the ground truth, which is the known presence of time cells as put into
the synthetic data.

```
0,68,2,0,65,67,1,1,66,68,3,0,64,64,1,4,66,29,1,39,66
1,67,4,1,63,68,0,0,67,68,4,0,63,63,0,5,67,29,0,39,67
2,67,2,1,65,66,1,2,66,68,3,0,64,61,2,7,65,41,2,27,65
```

The comma-separated entries are:

	- datasetIdx:	Index of the dataset, i.e. recording session.
	- A block of four entries for the sigMean for the TI/Mau method:
		- True Negative
		- False Negative
		- False Positive
		- True Positive
	- A similar block of four entries for the TI for the TI/Mau method
	- A similar block of four entries where both sigMean and TI must be true
	- A similar block of four entries for r2b threshold
	- A similar block of four entries for r2b bootstrap score.


## Contents of this directory:


### Python demo scripts:
- ti_demo.py	: Runs Mau's Temporal Information set of algorithms.<br/>
- r2b_demo.py	: Runs Modi's Ridge-to-background set of algorithms.<br/> 
- peq_demo.py	: Runs a simple classifier based on noise, hit-trial ratio, event width and jitter. Also reports these useful stats of the dataset.<br/> 
- ground_truth_check.py : Uses synthetic data files to assess accuracy of
			classification by the various Mau and Modi algorithms.<br/> 
- benchmark.py	: Simple time and memory benchmarks for the Mau and Modi
			algorithms.<br/>
- run_batch_analysis.py	: Runs a batch analysis on a datafile, generates 3 csv 
			files with the output. Note that the output files are
			created in append mode, so that previous runs are
			not overwritten. 

### Matlab(R) demo scripts:

We have placed a few Matlab demo scripts
here to show how the synthetic data files are generated, and to show how
the analysis can be run using Python wrapper functions called from Matlab.

run_batch_analysis.m	: Runs a partial batch analysis on a datafile,
	output is printed on console.

generateSyntheticData.m	: Generates a MATLAB data file with synthetic data
	representing time cell activity, according to several control 
	parameters.


### Files to build the pyBind11 module

r2bScore.cpp  
peqScore.cpp  
tcBind.cpp  
tcDefaults.cpp  
timeCell.cpp  
tiScore.cpp 
tcHeader.h
Makefile

### Matlab Interface Files

pyBindMap.py	: Provides an interface to the python/C++ functions using
			two wrapper functions:
	
	runTIanalysis
	runR2Banalysis



## pybind11 module functions

This module provides an interface for the Mau algorithms, the Modi algorithms,
the peq calculations, and the data structures used to configure them.

Functions:
----------

### tiScore

> CellScore tc.tiScore( data, analysisParams, tiAnalysisParams )

This performs the analysis of Mau et al 2018. It effectively handles
three classifications of time cells, all involving bootstrapping

	- **Returns**: an array of CellScore structures. See below.
	- **Arguments**: 
		data: Python Array of doubles organized as
		data[frameIdx][trailIdx][cellIdx]
	- analysisParams: AnalysisParams data structure, see below
	- tiAnalysisParams: TiAnalysisParams data structure, see below.

### r2bScore

> CellScore tc.r2bScore( data, analysisParams, threshold, percentile )
This performs the analysis of Modi et al 2014. It handles two
classification methods, one by thresholding and the other by bootstrapping.

	- **Returns**: an array of CellScore structures. See below.
	- **Arguments**: 
		data: Python Array of doubles organized as
		data[frameIdx][trailIdx][cellIdx]

	- analysisParams: AnalysisParams data structure, see below
	- threshold: Threshold for classifying a cell as a time cell using
		the original Modi 2014 method. A good value is 3.0.
	- percentile: Percentile cutoff for r2b bootstrap. Suggest 99.5.

### peqScore

> CellScore tc.peqScore( data, analysisParams, peqAnalysisParams )
This performs the parametric equation analysis defined in Ananthamurthy and Bhalla 2022.
It isn't a very good time-cell classifier, but it returns estimates of a number of key
parameters of the dataset: the noise, the event width, the jitter, and the hit trial ratio.

	- **Returns**: an array of CellScore structures. See below.
	- **Arguments**: 
		data: Python Array of doubles organized as
		data[frameIdx][trailIdx][cellIdx]

	- analysisParams: AnalysisParams data structure, see below
	- peqAnalysisParams: PeqAnalysisParams data structure, see below.



## pybind11 module data structures

This module provides an interface for the Mau algorithms, the Modi algorithms,
the PEQ analysis, and the data structures used to configure them.

Data structures:
----------------

### CellScore

CellScore is the class of the return object from each of the analysis routines.
It reports the stats for a given cell in a given session, over multiple trials.
Its fields are:

| Field name | Type    | Meaning in ti method | Meaning in r2b method | Meaning in peq method |
|------------|---------|----------------------|-----------------------|-----------------------|
| meanScore  | double  | Peak of mean over trials | r2b average of shuffled trials | mean of all frames and trials |
| baseScore  | double  | Temporal information of cell | r2b value of original trial | PEQ score |
| percentileScore  | double  | percentile for original TI among TI of shuffled trials | percentile for original r2b among r2b of shuffled trials | Ignored |
| sdev       | double  | Ignored              |   Ignored             | Standard deviation of all frames and trials |
| eventWidthMean | double  | Ignored          |   Ignored             | Mean of event width, in frames |
| eventWidthSdev | double  | Ignored          |   Ignored             | Mean of event width, in frames |
| imprecision | double     | Ignored          |   Ignored             | Imprecision in onset of event, in frames |
| sigMean  | bool  | does mean pk differ from shuffled? | is baseScore above r2b_threshold? | Ignored |
| sigBootstrap  | bool  | does TI differ from shuffled? | is r2b percentileScore above threshold? | Ignored |
| fracTrialsFired  | double  | Hit trial ratio, fraction of trials with sig response | not computed, set to 0 | Hit trial ratio |
| meanTrace   | array of doubles | Average activity vs time over all trials | Average activity vs time over all trials | Average activity vs time over all trials |
| meanPkIdx   | int | frame # of peak of mean over trials, typically bins of 3 frames | frame # of peak of mean over trials, no binning | frame # of peak of mean over trials, no binning |


### AnalysisParams

AnalysisParams objects are passed in to the tiScore and r2b methods to 
set parameters for the analysis.

| Field name     | Type  | default |         Meaning                   |
|----------------|-------|---------|-----------------------------------|
| csOnsetFrame   | int   | 75 | First frame of conditioned stimulus      |
| usOnsetFrame   | int   | 190 | First frame of unconditioned stimulus    |
| circPad        | int   | 20 | # frames on either side of cs-us to include in shuffling |
| circShuffleFrames | int   | 155 | # frames for circular shuffling |
| binFrames | int   | 3 | # frames for binning input trace, used in TI method |
| numShuffle | int   | 1000 | Number of times to do shuffling for bootstrap |
| epsilon | double   | 1.0e-6 | Smallest number for log operations |


### TiAnalysisParams

These are some additional parameters for the TI method from Mau 2018.

| Field name     | Type  | default |         Meaning                   |
|----------------|-------|---------|-----------------------------------|
| transientThresh   | double   | 2.0 | Number of sdev over mean that an event should be to qualify as a transient |
| tiPercentile   | double   | 99.0 | Percentile rank for TI to qualify |
| fracTrialsFiredThresh | double   | 0.25 | Fraction of trials where there must be a transient for the cell to qualify. |
| frameDt | double   | 0.08 | Duration of a 2p imaging frame for the Ca signal |

### PeqAnalysisParams

These are additional parameters for the PEQ method from Ananthamurthy and Bhalla 2022.

| Field name    | Type   | default |         Meaning                   |
|---------------|--------|---------|-----------------------------------|
| alpha         | double | 10.0    | Scale factor for noise term in PEQ calculation |
| beta          | double | 1.0     | Scale factor for event width term in PEQ calculation |
| gamma         | double | 10.0    | Scale factor for imprecision term in PEQ calculation |
| transientThresh | double   | 2.0 | Number of sdev over mean for an event to qualify as a transient |
| hitWindow     | double   | 5.0 | Number of frames + or - from peak to claim event as a hit |



## Installation:
	pip install TimeCellAnalysys


or

	git clone TimeCellAnalysis
	pip install h5py
	pip install pybind11
	make
