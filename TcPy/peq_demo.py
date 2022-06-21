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
 * File:            peq_demo.py
 * Description:     Example of use of timeCell python module on Matlab files
 *                  This extracts a range of parameters about the signal
 *                  including noise, jitter, hit-trial ratio and event width
 *                  This demo takes as input a MATLAB file in version 7.3.
 *                  It expects data in the form 
 *                  data[DATASET][CELL][TRIAL][FRAME]
 *                  It does a poor job of cell classification, but reports
 *                  this along with the extracted parameters.
 *                  It  reports the classification of cells as time/non time
 * Author:          Upinder S. Bhalla
 * E-mail:          bhalla@ncbs.res.in
 * Copyright (c) Upinder S. Bhalla  2022
 ********************************************************************/
 '''

import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
import tc       # This is the timeCell analysis code module.

DATA_LOCATION = "/sdo_batch/syntheticDATA"   # Synthetic dataset.
PEQ_THRESH = 0.002

def analyzeDatasets( dat ):
    ap = tc.AnalysisParams()        # Use defaults for AnalysisParams
    pep = tc.PeqAnalysisParams()    # Use defaults for PeqAnalysisParams
    sd0 = dat[DATA_LOCATION]        # Datasets.

    # Go through all entries in synthetic dataset. Each corresponds to
    # a recording session with different conditions of noise, background...
    for idx, ss in enumerate( sd0 ):
        # r2bScore returns an array of class CellScore.
        peqScore = np.array( tc.peqScore( dat[ss[0]], ap, pep) )

        # Print classification of first 30 cells for first dataset
        if idx == 0:   
            print( "Stats for first 30 cells for Dataset 0" )
            print( "CellIdx  sigPEQ      noise     eventwidth   imprecision    hit trial ratio   ")
            for cellIdx, pp in enumerate( peqScore ):
                print( "{:6d}{:8d}{:12.2f}{:12.2f}{:14.2f}{:12.2f}".format( cellIdx, (pp.baseScore > PEQ_THRESH), pp.sdev/pp.meanScore, pp.eventWidthMean, pp.imprecision, pp.fracTrialsFired ) )
                if cellIdx >= 30:
                    break
            print( "Number of time cells classified for each dataset")
            print( "Dataset    # time cells" )
        # Print number of classified cells for all the datasets.
        numSig = 0
        for pp in peqScore:
            numSig += (pp.baseScore > PEQ_THRESH)

        print( "{:4d}{:12d}".format( idx, numSig), flush = True )

def main():
    parser = argparse.ArgumentParser( description = "Perform r2b time cell analysis on given matlab file" )
    parser.add_argument( "datafile",  type = str, help = "Required. File name to load, in matlab format" )
    args = parser.parse_args()

    dat = h5py.File( args.datafile, 'r' )
    analyzeDatasets( dat )

if __name__ == '__main__':
    main()
