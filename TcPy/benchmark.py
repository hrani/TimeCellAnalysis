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
 * File:            benchmark.py
 * Description:     Benchmark timeCell python module on Matlab files
 * Author:          Upinder S. Bhalla
 * E-mail:          bhalla@ncbs.res.in
 * Copyright (C) Upinder S. Bhalla
 ********************************************************************/
 '''

import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
import tracemalloc
import time
import tc       # This is the timeCell analysis code module.

R2B_THRESH = 3.0
R2B_PERCENTILE = 99.5
NUM_ITER = 10

def doBenchmark( dat, method ):
    ap = tc.AnalysisParams()        # Use defaults for AnalysisParams
    tip = tc.TiAnalysisParams()     # Use defaults for TIAnalysisParams
    tip.frameDt = 1.0/ 12.5         # Reassign default frameDt
    pep = tc.PeqAnalysisParams()    # Use defaults for Parameter Eqn.
    sd0 = dat["/sdo_batch/syntheticDATA"]   # Synthetic dataset.
    t0 = time.time()

    # Go through all entries in synthetic dataset. Each corresponds to
    # a recording session with different conditions of noise, background...
    if method == "ti":
        for idx in range( NUM_ITER ):
            tiScore = np.array(tc.tiScore( dat[sd0[idx][0]], ap, tip ) )
            print( idx, end = " ", flush = True )

    elif method == "r2b":
        for idx in range( NUM_ITER ):
            r2bScore = np.array( tc.r2bScore( dat[sd0[idx][0]], ap, R2B_THRESH, R2B_PERCENTILE ) )
            print( idx, end = " ", flush = True )
    elif method == "peq":
        for idx in range( NUM_ITER ):
            peqScore = np.array( tc.peqScore( dat[sd0[idx][0]], ap, pep ) )
            print( idx, end = " ", flush = True )
    print()
    return time.time() - t0

def main():
    tracemalloc.start()
    parser = argparse.ArgumentParser( description = "Perform r2b time cell analysis on given matlab file" )
    parser.add_argument( "datafile",  type = str, help = "Required. File name to load, in matlab format" )
    parser.add_argument( "method",  type = str, help = "Required. Options: 'ti' or 'r2b' or 'peq'" )
    args = parser.parse_args()

    dat = h5py.File( args.datafile, 'r' )
    runtime = doBenchmark( dat, args.method )
    print( "Memory Use: ", tracemalloc.get_traced_memory()[1], " bytes" )
    print( "Time per dataset = ", runtime / NUM_ITER, " sec" )
    tracemalloc.stop()

if __name__ == '__main__':
    main()
