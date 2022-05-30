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
 * File:            ti_demo.py
 * Description:     Example of use of timeCell python module on Matlab files
 *                  This uses the analysis from 
 *                  Mau et al, Curr Biol. 2018 
 *                  28(10):1499-1508.e4. doi: 10.1016/j.cub.2018.03.051.
 *                  This demo takes as input a MATLAB file in version 7.3.
 *                  It expects data in the form 
 *                  data[DATASET][CELL][TRIAL][FRAME]
 *                  It  reports the classification of cells as time/non time
 *                  cells by 3 sub-methods in the Mau analysis.
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
    tip = tc.TiAnalysisParams()     # Use defaults for TIAnalysisParams
    tip.frameDt = 1.0/ 12.5         # Reassign default frameDt
    sd0 = dat[DATA_LOCATION]        # Datasets.

    # Go through all entries in synthetic dataset. Each corresponds to
    # a recording session with different conditions of noise, background...
    for idx, ss in enumerate( sd0 ):
        # tiScore returns an array of class CellScore.
        tiScore = np.array(tc.tiScore( dat[ss[0]], ap, tip ) )
        # Print classification of first 30 cells for first dataset
        if idx == 0:   
            print( "Classification of first 30 cells for Dataset 0" )
            print( "CellIdx sigMean    SigTI   SigBoth pkFrame   FracTrialsFired")
            for cellIdx, ts in enumerate( tiScore ):
                print( "{:6d}{:8d}{:8d}{:8d} {:8d} {:12.4f}".format( cellIdx, ts.sigMean, ts.sigBootstrap, ts.sigMean and ts.sigBootstrap, ts.meanPkIdx, ts.fracTrialsFired ) )
                if cellIdx >= 30:
                    break
            print( "\nNumber of time cells classified by each method, for each dataset")
            print( "Dataset    #SigMean    #sigBoot     #sigBoth" )
        # Print number of classified cells for all the datasets.
        numMean = 0
        numBoot = 0
        numBoth = 0
        for ts in tiScore:
            numMean += ts.sigMean
            numBoot += ts.sigBootstrap
            numBoth += (ts.sigBootstrap and ts.sigMean)

        print( "{:4d}{:12d}{:12d}{:12d}".format( idx, numMean, numBoot, numBoth ), flush = True )

def main():
    parser = argparse.ArgumentParser( description = "Perform ti (Mau 2018) time cell analysis on given matlab file" )
    parser.add_argument( "datafile",  type = str, help = "Required. File name to load, in matlab format" )
    args = parser.parse_args()

    dat = h5py.File( args.datafile, 'r' )
    analyzeDatasets( dat )

if __name__ == '__main__':
    main()
