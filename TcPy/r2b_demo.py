# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
# 

'''
*******************************************************************
 * File:            r2b_demo.py
 * Description:     Example of use of timeCell python module on Matlab files
 *                  This uses the analysis from 
 *                  Modi, Dhawale and Bhalla, eLife, 2014
 *                  3:e01982. doi: 10.7554/eLife.01982.
 *                  This demo takes as input a MATLAB file in version 7.3.
 *                  It expects data in the form 
 *                  data[DATASET][CELL][TRIAL][FRAME]
 *                  It  reports the classification of cells as time/non time
 *                  cells by 2 sub-methods in the Modi et al analysis.
 * Author:          Upinder S. Bhalla
 * E-mail:          bhalla@ncbs.res.in
 * Copyright (c) Upinder S. Bhalla
 ********************************************************************/
 '''

import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
import tc       # This is the timeCell analysis code module.

DATA_LOCATION = "/sdo_batch/syntheticDATA"   # Synthetic dataset.
R2B_THRESH = 3.0
R2B_PERCENTILE = 99.5

'''
# These are the data structures. Params go into the function, and 
# CellScore comes out. Default values are indicated here.
# These are initialized in C++, shown here for clarity.
class AnalysisParams():
    def __init__( self ):
        self.csOnsetFrame = 75
        self.usOnsetFrame = 190
        self.circPad = 20
        self.circShuffleFrames = 40 + 190 - 75
        self.binFrames = 3
        self.numShuffle = 1000
        self.epsilon = 1.0e-6

class TiAnalysisParams():
    def __init__( self ):
        self.transientThresh = 2.0
        self.tiPercentile = 99.0
        self.fracTRialsFiredThresh = 0.25
        self.frameDt = 1.0 / 12.5

# Note that CellScore is read-only. Its values are filled by the tc code.
#class CellScore():
#    float self.meanScore        #Mau: pk of mean trace. r2b: shuffled mean
#    float self.baseScore        # Mau: Raw TI score, raw r2b ratio
#    float self.percentileScore  # Mau: Temporal Info. r2b: bootstrap score
#    bool self.sigMean           # Mau: Is mean sig. r2b: Is mean ratio sig?
#    bool self.sigBootstrap      # Mau and r2b: Is over bootstrap thresh.
#    float self.fracTrialsFired  # Hit trial ratio.
#    np.array meanTrace          # meanTrace[frame#]. Ave trials for a cell
#    int meanPkIdx               # Idx of peak frame in above.
'''

def analyzeDatasets( dat ):
    ap = tc.AnalysisParams()        # Use defaults for AnalysisParams
    sd0 = dat[DATA_LOCATION]        # Datasets.

    # Go through all entries in synthetic dataset. Each corresponds to
    # a recording session with different conditions of noise, background...
    for idx, ss in enumerate( sd0 ):
        # r2bScore returns an array of class CellScore.
        r2bScore = np.array( tc.r2bScore( dat[ss[0]], ap, R2B_THRESH, R2B_PERCENTILE ) )

        # Print classification of first 30 cells for first dataset
        if idx == 0:   
            print( "Classification of first 30 cells for Dataset 0" )
            print( "CellIdx sigMean    SigTI   pkFrame   ")
            for cellIdx, rr in enumerate( r2bScore ):
                print( "{:6d}{:8d}{:8d} {:8d}".format( cellIdx, rr.sigMean, rr.sigBootstrap, rr.meanPkIdx ) )
                if cellIdx >= 30:
                    break
            print( "\nNumber of time cells classified by each method, for each dataset")
            print( "Dataset    #SigMean    #sigBoot" )
        # Print number of classified cells for all the datasets.
        numMean = 0
        numBoot = 0
        for rr in r2bScore:
            numMean += rr.sigMean
            numBoot += rr.sigBootstrap

        print( "{:4d}{:12d}{:12d}".format( idx, numMean, numBoot), flush = True )

def main():
    parser = argparse.ArgumentParser( description = "Perform r2b time cell analysis on given matlab file" )
    parser.add_argument( "datafile",  type = str, help = "Required. File name to load, in matlab format" )
    args = parser.parse_args()

    dat = h5py.File( args.datafile, 'r' )
    analyzeDatasets( dat )

if __name__ == '__main__':
    main()
