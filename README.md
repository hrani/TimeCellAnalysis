# Time Cell Analysis project: Python and pybind11/C++ files for analysis and demos.


## Overview

This repository contains the Matlab, Python and related files for 
Time Cell Analysis,

This is from a forthcoming paper: 

Synthetic Data Resource and Benchmarks for Time Cell Analysis and Detection Algorithms
*K. Ananthamurthy and U.S. Bhalla,* 
in preparation.

## Description
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
	Please see README in TcPy for details on running demos etc.
	- rho-matlab: Time Cell analysis Matlab libraries.
	Please see README in rho-matlab for details on running demos.


## Installation:
	pip install TimeCellAnalysis

Following this command, the actual files will be placed in

	~/.local/bin/TimeCellAnalysis


or

	git clone TimeCellAnalysis
	pip install h5py
	pip install pybind11
	pip install matplotlib
	make


Following this command, the actual files will be placed in 

	./TimeCellAnalysis

