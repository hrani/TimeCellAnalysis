#Time Cell Analysis files: Python and pybind11/C++ files for analysis and demos.


##Overview

This repository contains the Matlab, Python and related files for 
Time Cell Analysis,
Kambadur Ananthamurthy and U.S. Bhalla, in preparation.


##Description
Time cells are neurons whose activity encodes the time since a reference 
stimulus. They have been observed in the hippocamal CA1, CA3, and also in
entorhinal cortex of rodents (refs). They may encode times of the order of
100ms (Modi et al 2014 eLife ) to 20s (Mau et al 2018 Current Biology).

Several algorithms have been developed to identify time cells from amongst a
population of firing neurons. With the advent of large-scale unit recordings
using 2-photon Calcium imaging or high-density electrodes, it is important to
have reliable ways to identify time cells automatically.

This project has implemented a way to assess the performance of time-cell
algorithms. We have done two key things. First, we implemented code to generate
synthetic neuronal activity data in which we know the ground truth of which 
cells are time-cells, and we can control parameters such as noise, background 
activity, jitter, and hit trial ratio (fraction of trials in which the time 
cell was active). 
Second, we implemented and extended published time-cell analysis algorithms. 
While some of the original algorithms were in Matlab, we have re-implemented 
key ones in C++ using the pybind11 libraries to provide a simple Python 
interface. This gives us considerable improvements in speed and memory 
efficiency, at the cost of some complexity in the code.
The Python functions can also be accessed via Matlab, and we illustrate how
this is done.

## Directories:

	- TcPy: Time Cell analysis Python demos, pybind11 and example driver 
	code from Matlab.
	- rho-matlab: Time Cell analysis Matlab libraries.

## Examples:
---------
For the impatient, once you have installed the files, here are some things
to do from the TcPy directory:

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

Similarly for *r2b_demo.py*


If your source data file is generated from the synthetic data program, you 
can run:

> python ground_truth_check.py <filename>

and it will use ground-truth data from the file to see how well the various
methods handle the data. Output is like this:

```
idx  Mau Peak Sig |   Mau TI    | Mau pk&&TI  |  r2b thresh | r2b bootstr
idx | tn fn fp tp | tn fn fp tp | tn fn fp tp | tn fn fp tp | tn fn fp tp
  0 | 68  2  0 65 | 68  1  0 66 | 68  3  0 64 | 64  1  4 66 | 29  1 39 66
  1 | 67  5  1 62 | 68  0  0 67 | 68  5  0 62 | 63  0  5 67 | 29  0 39 67
  2 | 67  2  1 65 | 66  1  2 66 | 68  3  0 64 | 61  2  7 65 | 41  2 27 65
...
```

If you want to run a batch analysis on a synthetic data file do the following. 
It will dump all the data into .csv files:

> python run_batch_analysis.py *filename*

This will generate the files `ti.csv`, `r2b.csv` and `groundTruth.csv`

Here are three lines from a sample ti,csv:

```
0,0,0.4535,-0.1383,-0.1332,0,0,0.6000,0
0,1,0.1643,-0.0821,-0.0700,0,0,0.3833,3
0,2,2.6316,-0.1009,-0.1673,1,1,0.7167,7
```

The comma-separated entries are

	- datasetIdx:	Index of the dataset, i.e. recording session.
	- cellIdx:	Index of the cell whose stats are reported on this line
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
	- meanScore: 	r2b score for shuffled trials.
	- baseScore:	r2b score for original trial.
	- percentileScore:	where does r2bscore place among shuffled r2bs?
	- sigMean:	Is the base score above a threshold? Usually 3.0.
	- sigBootstrap:	Is the percentile score above a threshold? Usually 99.5
	- fracTrialsFired:	Not computed, set to 0.
	- meanPkIdx:	Frame number of peak of mean cell activity. Note that this method does not do frame binning.

finally, groundTruth.csv reports how well the different methods compare to
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
ti_demo.py	: Runs Mau's Temporal Information set of algorithms.
r2b_demo.py	: Runs Modi's Ridge-to-background set of algorithms.
ground_truth_check.py : Uses synthetic data files to assess accuracy of
			classification by the various Mau and Modi algorithms.
benchmark.py	: Simple time and memory benchmarks for the Mau and Modi
			algorithms.

run_batch_analysis.py	: Runs a batch analysis on a datafile, generates 3 csv 
			files with the output.

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
tcBind.cpp  
tcDefaults.cpp  
timeCell.cpp  
tiScore.cpp 
tcHeadher.h
Makefile

### Matlab Interface Files

pyBindMap.py	: Provides an interface to the python/C++ functions using
			two wrapper functions:
	
	runTIanalysis
	runR2Banalysis



## pybind11 module functions

This module provides an interface for the Mau algorithms, the Modi algorithms,
and the data structures used to configure them.

Functions:

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
	classification methods, one by thresholding and the other by 
	bootstrapping.
	**Returns**: an array of CellScore structures. See below.
	**Arguments**: 
		data: Python Array of doubles organized as
		data[frameIdx][trailIdx][cellIdx]

		analysisParams: AnalysisParams data structure, see below
		threshold: Threshold for classifying a cell as a time cell using
		the original Modi method. A good value is 3.0.
		percentile: Percentile cutoff for r2b bootstrap. Suggest 99.5.





Installation:
-------------
	pip install TimeCellAnalysys


or

	git clone TimeCellAnalysis
	pip install h5py
	pip install pybind11
	make
