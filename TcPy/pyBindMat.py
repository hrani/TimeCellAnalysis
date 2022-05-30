import numpy as np
import tc

r2bThresh = 3.0
r2bPercentile = 99.5
NUM_CELLS = 135

'''
### The data structures below are now obtained from the C++ code.

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


'''

def convScoreToList( cs ):
    ret = [ cs.meanScore, cs.baseScore, cs.percentileScore, cs.sigMean, cs.sigBootstrap, cs.fracTrialsFired, cs.meanPkIdx ]
    return ret

def runR2Banalysis( DATA, csOnset, usOnset, circPad, circShuffleFrames, binFrames, numShuffle, r2bThresh, r2bPercentile ):
    ap = tc.AnalysisParams()
    ap.csOnsetFrame = int( csOnset )
    ap.usOnsetFrame = int( usOnset )
    ap.circPad = int( circPad )
    ap.circShuffleFrames = int( circShuffleFrames )
    ap.binFrames = int( binFrames )
    ap.numShuffle = int( numShuffle )
    # This assumes that the axis ordering is DATA[Cell][trial][frame]
    dat = np.swapaxes( np.array( DATA ), 0, 2 )
    # r2bScore assumes that the axis ordering is DATA[frame][trial][Cell]
    cs = tc.r2bScore( dat, ap, r2bThresh, r2bPercentile )
    ret = []
    for idx, rr in enumerate( cs ):
        ret.append( convScoreToList( rr ) )
    return ret

def runTIanalysis( DATA, csOnset, usOnset, circPad, circShuffleFrames, binFrames, numShuffle,  transientThresh, tiPercentile, fracTrialsFiredThresh, frameDt ):
    ap = tc.AnalysisParams()
    ap.csOnsetFrame = int( csOnset )
    ap.usOnsetFrame = int( usOnset )
    ap.circPad = int( circPad )
    ap.circShuffleFrames = int( circShuffleFrames )
    ap.binFrames = int( binFrames )
    ap.numShuffle = int( numShuffle )

    tip = tc.TiAnalysisParams()
    tip.transientThresh = transientThresh
    tip.tiPercentile = tiPercentile
    tip.fracTrialsFiredThresh = fracTrialsFiredThresh
    tip.frameDt = frameDt

    dat = np.swapaxes( np.array( DATA ), 0, 2 )
    cs = tc.tiScore( dat, ap, tip )
    ret = []
    for idx, rr in enumerate( cs ):
        ret.append( convScoreToList( rr ) )
    return ret

def main():
    parser = argparse.ArgumentParser( description = "Perform r2b time cell analysis on given matlab file" )
    parser.add_argument( "datafile",  type = str, help = "Required. File name to load, in matlab format" )
    args = parser.parse_args()

    dat = h5py.File( args.datafile, 'r' )
    printDatasetInfo( dat )

'''
if __name__ == '__main__':
    main()
'''
