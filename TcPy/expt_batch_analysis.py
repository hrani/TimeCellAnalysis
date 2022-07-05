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
 * File:            run_batch_analysis.py
 * Description:     Do timeCell analysis on Matlab files, in batch mode.
 *                  Generates 3 output files in csv format:
 *                  ti.csv, r2b.csv, groundTruth.csv.
 *                  References: Mau et al. 2018 Current Biology
 *                              Modi et al. 2014 eLife.
 * Author:          Upinder S. Bhalla
 * E-mail:          bhalla@ncbs.res.in
 * Copyright (c) Upinder S. Bhalla
 ********************************************************************/
 '''

import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
import time
from scipy.io import loadmat
import tc       # This is the timeCell analysis code module.

R2B_THRESH = 3.0
R2B_PERCENTILE = 99.5
PEQ_THRESH = 0.002

# class AnalysisParams() and class TiAnalysisParams are described in the
# README.md
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

def scoreString( datasetIdx, cellIdx, x ):
    return "{},{},{:.4f},{:.4f},{:.4f},{:1d},{:1d},{:.4f},{}\n".format( datasetIdx, cellIdx, x.meanScore, x.baseScore, x.percentileScore, int(x.sigMean), int(x.sigBootstrap), x.fracTrialsFired, x.meanPkIdx )

def peqScoreString( datasetIdx, cellIdx, x ):
    if x.meanScore > 1e-5:
        noise = x.sdev/x.meanScore
    else:
        noise = 0.0
    return "{},{},{:1d},{:.6f},{:.4f},{:.4f},{:.4f},{:.4f}\n".format( datasetIdx, cellIdx, int(x.baseScore > PEQ_THRESH), x.baseScore, noise, x.eventWidthMean, x.imprecision, x.fracTrialsFired )

def convertAP( args ):
    ap = tc.AnalysisParams()        # Use defaults for AnalysisParams
    ap.csOnsetFrame = args.stimulusFrames[0]
    ap.usOnsetFrame = args.stimulusFrames[1]
    ap.circPad = args.circPad
    print( "AP BIN = ", args.binFrames )
    ap.binFrames = args.binFrames
    ap.numShuffle = args.numShuffle
    ap.epsilon = args.epsilon

    return ap

def printV5Info( dat, ap, r2bFile, tiFile, peqFile  ):
    pep = tc.PeqAnalysisParams()    # Use defaults for PeqAnalysisParams
    tip = tc.TiAnalysisParams()     # Use defaults for TIAnalysisParams
    tip.frameDt = 1.0/ 12.5         # Reassign default frameDt

    arg = np.swapaxes( dat, 0, 2 )
    tiScore = np.array(tc.tiScore( arg, ap, tip ) )
    r2bScore = np.array( tc.r2bScore( arg, ap, R2B_THRESH, R2B_PERCENTILE ) )
    peqScore = np.array( tc.peqScore( arg, ap, pep ) )
    idx = 0

    for cellIdx, (tt, rr, pp) in enumerate( zip( tiScore, r2bScore, peqScore ) ):
        tiFile.write( scoreString( idx, cellIdx, tt ) )
        r2bFile.write( scoreString( idx, cellIdx, rr ) )
        peqFile.write( peqScoreString( idx, cellIdx, pp ) )


def printHDF5Info( dat, location, ap, r2bFile, tiFile, peqFile  ):
    pep = tc.PeqAnalysisParams()    # Use defaults for PeqAnalysisParams
    tip = tc.TiAnalysisParams()     # Use defaults for TIAnalysisParams
    tip.frameDt = 1.0/ 12.5         # Reassign default frameDt

    sd0 = dat[location]             # dataset.

    # Go through all entries in synthetic dataset. Each corresponds to
    # a recording session with different conditions of noise, background...
    for idx, ss in enumerate( sd0 ):
        # These are the calls to the analysis routines. Return is an
        # array of CellScores, see above
        tiScore = np.array(tc.tiScore( dat[ss[0]], ap, tip ) )
        r2bScore = np.array( tc.r2bScore( dat[ss[0]], ap, R2B_THRESH, R2B_PERCENTILE ) )
        peqScore = np.array( tc.peqScore( dat[ss[0]], ap, pep ) )

        for cellIdx, (tt, rr, pp) in enumerate( zip( tiScore, r2bScore, peqScore ) ):
            tiFile.write( scoreString( idx, cellIdx, tt ) )
            r2bFile.write( scoreString( idx, cellIdx, rr ) )
            peqFile.write( peqScoreString( idx, cellIdx, pp ) )
        print( "Done dataset ", idx, flush=True)

def printDatasetInfo( dat, location, isHDF5, ap ):
    t0 = time.time()

    r2bFile = open("r2b.csv", "w")
    tiFile = open("ti.csv", "w")
    peqFile = open("peq.csv", "w")

    if isHDF5:
        printHDF5Info( dat, location, ap, r2bFile, tiFile, peqFile )
    else:
        printV5Info( dat[location], ap, r2bFile, tiFile, peqFile )
        
    # Print it all out.
    tiFile.write( "\n" ) # Put a spacer in case we rerun and append to file.
    r2bFile.write( "\n" )
    peqFile.write( "\n" )
    tiFile.close()
    r2bFile.close()
    peqFile.close()

def main():
    '''
    Perform time cell analysis on given matlab file, generate output csv files.\n
    File contents for ti,csv and r2b.csv are:\n
    datasetIdx, cellIdx, meanScore, baseScore, percentileScore, sigMean, sigBootstrap, fracTrialsFired, meanPkIdx\n
    Note that the fracTrialsFired is not computed for the r2b method.\n
    '''
    parser = argparse.ArgumentParser( description = main.__doc__ )
    parser.add_argument( "datafile",  type = str, help = "Required. File name to load, either in matlab >= 7.3 aka HDF5 format or in Matlab v5." )
    parser.add_argument( "-d", "--dataLocation",  type = str, help = "Optional. String to look up data in HDF file.", default = "/sdo_batch/syntheticDATA" )
    parser.add_argument( "-s", "--stimulusFrames",  type = int, nargs=2, help = "Optional. CS_frame US_frame. Default = [75,190]", default = [75,190], metavar = ("CS_frame", "US_frame" ) )
    parser.add_argument( "--circPad",  type = int, help = "Optional. Number of frames for padding circular shuffle, default=20.", default = 20 )
    parser.add_argument( "-b", "--binFrames",  type = int, help = "Optional. Bin time-series by this number of frames, default=3.", default = 3 )
    parser.add_argument( "--numShuffle",  type = int, nargs=1, help = "Optional. Number of times to shuffle dataset for bootstrap, default=1000", default = 1000 )
    parser.add_argument( "--epsilon",  type = float, nargs=1, help = "Optional. Espilon to use as smallest allowed value for mean df by f.", default = 1e-6 )
    args = parser.parse_args()

    ap = convertAP( args )

    if h5py.is_hdf5( args.datafile ):
        dat = h5py.File( args.datafile, 'r' )
        printDatasetInfo( dat, args.dataLocation, True, ap )
    else:
        dat = loadmat( args.datafile )
        printDatasetInfo( dat, args.dataLocation, False, ap )

if __name__ == '__main__':
    main()
