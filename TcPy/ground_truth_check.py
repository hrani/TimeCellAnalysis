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
 * File:            timeCellMatlabExample.py
 * Description:     Example of use of timeCell python module on Matlab files
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

R2B_THRESH = 3.0
R2B_PERCENTILE = 0.995

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

def printDatasetInfo( dat ):
    ap = tc.AnalysisParams()        # Use defaults for AnalysisParams
    tip = tc.TiAnalysisParams()     # Use defaults for TIAnalysisParams
    tip.frameDt = 1.0/ 12.5         # Reassign default frameDt
    ptc = dat["sdo_batch/ptcList"] # ptc is list of positive time cells.
    sd0 = dat["/sdo_batch/syntheticDATA"]   # Synthetic dataset.

    # Go through all entries in synthetic dataset. Each corresponds to
    # a recording session with different conditions of noise, background...
    print( "idx  Mau Peak Sig |   Mau TI    | Mau pk&&TI  |  r2b thresh | r2b bootstr" )
    for idx, ss in enumerate( sd0 ):
        #truth array indexed as: truth[scoringMethod][isCorrect][isPositive]
        truth = np.zeros( (5,2,2), dtype=int)
        # These are the calls to the analysis routines. Return is an
        # array of CellScores, see above
        tiScore = np.array(tc.tiScore( dat[ss[0]], ap, tip ) )
        r2bScore = np.array( tc.r2bScore( dat[ss[0]], ap, R2B_THRESH, R2B_PERCENTILE ) )

        # Fill in the groundTruth[cell#] array: True if cell is time-cell.
        groundTruth = np.zeros( len( tiScore ), dtype = int )
        trueCells = np.array( dat[ ptc[idx][0] ], dtype=int )[:,0]
        for cc in trueCells:
            groundTruth[cc-1] = 1   # Convert from 1-base to 0-base arrays

        # Go through and count true pos, false pos, true neg, false neg
        # for each of 5 classification methods in the Truth table.
        for tt, rr, gg in zip( tiScore, r2bScore, groundTruth ):
            truth[0][int(tt.sigMean)][gg] += 1      # Transient score
            truth[1][int(tt.sigBootstrap)][gg] += 1 # TI score

            # Mau full form: Transient AND TI.
            truth[2][int(tt.sigMean and tt.sigBootstrap)][gg] += 1

            # R2B threshold form
            truth[3][int(rr.sigMean)][gg] += 1

            # R2B bootstrap.
            truth[4][int(rr.sigBootstrap)][gg] += 1

        
        # Print it all out.
        if idx % 20 == 0:
            print( "idx | tn fn fp tp | tn fn fp tp | tn fn fp tp | tn fn fp tp | tn fn fp tp" )
        print( "{:3d}".format( idx ), end = "" )
        for tr in truth:
            x = np.array(tr)
            x.shape = (4,)
            print( " |{:3d}{:3d}{:3d}{:3d}".format(x[0], x[1], x[2], x[3]), end = "" )
        print( flush = True )

def main():
    parser = argparse.ArgumentParser( description = "Perform r2b time cell analysis on given matlab file" )
    parser.add_argument( "datafile",  type = str, help = "Required. File name to load, in matlab format" )
    args = parser.parse_args()

    dat = h5py.File( args.datafile, 'r' )
    printDatasetInfo( dat )

if __name__ == '__main__':
    main()
