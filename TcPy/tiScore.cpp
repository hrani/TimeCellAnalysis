#include <vector>
#include <random>
#include <numeric>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Python.h>

using namespace std;

namespace py = pybind11;
#include "tcHeader.h"
PYBIND11_MAKE_OPAQUE(std::vector<double>);

/// Returns the mean dfbf[ frame# ] for the specified shuffle[ trial# ]
/// using data[ trial# ][ frame# ]. Operates on data for a single cell.
vector< double > aveOfTrials( 
				const vector< vector< double > >& data, 
				const vector< unsigned int >& shuff, 
				unsigned int numFrames )
{
	unsigned int numTrials = data.size();
	vector< double > aveBin( numFrames, 0.0 );
	assert( numTrials == data.size() && numTrials > 0  );
	for ( unsigned int ii = 0; ii < numTrials; ++ii ) {
		vector< double >::iterator aptr = aveBin.begin();
		auto dd = data[ii].begin();
		for ( unsigned int bb = shuff[ii]; bb < shuff[ii] + numFrames; bb++, aptr++) {
			*aptr += *( dd + bb%numFrames );
		}
	}
	for( vector< double >::iterator aptr = aveBin.begin(); aptr != aveBin.end(); aptr++ ) {
		*aptr /= numFrames;
	}
	return aveBin;
}

/// returns binned data for a single trial.
vector< double > binner( const vector< double >& orig, unsigned int binFrames )
{
	if (binFrames <= 1 || binFrames >= orig.size()/2 ) return orig;
	unsigned int numRet = orig.size() / binFrames;
	vector< double > ret( numRet, 0.0 );
	for ( unsigned int idx = 0; idx < numRet; idx++ ) {
		for ( unsigned int bf = 0; bf < binFrames; bf++ )
			ret[idx] += orig[ idx * binFrames + bf ];
	}
	return ret;
}

/// returns binned data for a single cell, all trials.
vector< vector< double > > binnedCellData( const vector< vector< double > >& data, unsigned int numFrames, unsigned int binFrames )
{
	unsigned int numTrials = data.size();
	vector< vector< double > > ret( numTrials );
	for( unsigned int trialIdx = 0; trialIdx < numTrials; trialIdx++ ) {
		ret[ trialIdx ] = binner( data[ trialIdx ], binFrames );
	}
	return ret;
}

vector< unsigned int> findPeaks( const vector< double >& data ) {
	double pk = data[0];
	unsigned int pkIdx = 0;
	unsigned int otherPkIdx = 0;
	unsigned int numBins = data.size();
	for( unsigned int ii = 1; ii < numBins; ++ii ) {
		if ( pk < data[ii] ) {
			pk = data[ii];
			pkIdx = ii;
		}
	}
	if (pkIdx == 0) {
		otherPkIdx = 1;
	} else if (pkIdx == numBins - 1) {
		otherPkIdx = numBins - 2;
	} else if ( data[pkIdx-1] > data[pkIdx+1] ) {
		otherPkIdx = pkIdx-1;
	} else {
		otherPkIdx = pkIdx+1;
	}
	vector< unsigned int > ret( {pkIdx, otherPkIdx} );
	return ret;
}

/// Return the average tuning curve of the cell, and pass back an arg with
/// the bin# of the peak iff significant. Else return ~0U.
void tuningCurve( const vector< vector< double >>& data, const vector< vector< unsigned int > >& shuff, unsigned int numFrames, const AnalysisParams& ap, CellScore& cs )
{
	// Get the binned data
	unsigned int numTrials = data.size();
	vector< vector< double >> binnedCellData_ = binnedCellData( data, numFrames, ap.binFrames );
	unsigned int numBins = binnedCellData_[0].size();
	vector< unsigned int > zeros( numTrials, 0 );
	// Get the reference averaged binned trial data.
	cs.meanTrace = aveOfTrials( binnedCellData_, zeros, numBins );

	// Find the 2 consecutive biggest. Assume peak is one of them.
	vector< unsigned int > pks = findPeaks( cs.meanTrace );
	cs.meanPkIdx = pks[0];

	// Now I could be efficient and just average over the two bins
	// here, as per shuffling. Or I could use the simple aveOfTrials.
	unsigned int numOK = 0;
	for ( unsigned int ii = 0; ii < ap.numShuffle; ii++ ) {
		vector< double > temp = aveOfTrials( binnedCellData_, shuff[ii], numBins );
		numOK += ( cs.meanTrace[pks[0]] > temp[pks[0]] ) && 
				( cs.meanTrace[pks[1]] > temp[pks[1]] );
	}
	cs.meanScore = cs.meanTrace[cs.meanPkIdx];
	cs.sigMean = (  (numOK * 100) > (99 * ap.numShuffle) );
}

/// args: data[ frame# ][ trial# ][ cell# ]
/// Pass-back arg: ret[ cell# ][ trial# ][ frame# ]
/// Fills in only the data in the window around the stimuli, from 
/// CS_ONSET_FRAME - CIRC_PAD to US_ONSET_FRAME + CIRC_PAD.
void reorderData( const double* data, unsigned int numCells, unsigned int numTrials, unsigned int numFrames, vector< vector< vector< double >>>& ret, const AnalysisParams& ap ) 
{
	ret.clear();
	ret.resize( numCells );
	unsigned int ncxnt = numCells * numTrials;
	unsigned int nf = ap.circShuffleFrames;
	for ( unsigned int cell = 0; cell < numCells; cell++ ) {
		vector< vector< double > >& rc = ret[cell];
		rc.resize( numTrials );
		for( unsigned int tt = 0; tt < numTrials; tt++ ) {
			vector< double >& rct = rc[ tt ];
			rct.resize( nf,  0.0 );
			for( unsigned int ff = 0; ff < nf; ++ff ) {
				rct[ff] = data[ cell + tt * numCells + (ff + ap.csOnsetFrame - ap.circPad) * ncxnt];
			}
		}
	}
}

// Returns threshold for sig peak for data[trial#][frame#]. The
// data has already been sub-selected for a given cell.
double findCellThresh(
		const vector< vector< double > >& data, unsigned int numFrames, double transientThresh )
{
	double sum = 0.0;
	double sq = 0.0;
	for( auto trialD = data.begin(); trialD != data.end(); ++trialD ) {
		for( auto frameD = trialD->begin(); frameD != trialD->end(); ++frameD ) {
			sum += *frameD;
			sq += *frameD * *frameD;
		}
	}
	double numSamples = numFrames * data.size(); //numFrames * numTrials
	double mean = sum/numSamples;
	return mean + transientThresh * sqrt( sq / numSamples - mean*mean );
}

/// Returns time stamps of transients trial# ][ transient# ]. 
// given the reordered dfbf data[ trial# ][ frame# ]
// A given trial may have zero or more transients.
vector< vector< unsigned int > > findTransients(
		const vector< vector< double > >& data, double cellThresh, unsigned int numFrames )
{
	unsigned int numTrials = data.size();
	vector< vector< unsigned int > > ret ( numTrials );
	vector< vector< unsigned int > >::iterator trialT = ret.begin();
	for( auto trialD = data.begin(); trialD != data.end(); ++trialD, ++trialT ) {
		trialT->clear();
		double lastFrame = 0.0;
		bool refractory = 0;
		for( unsigned int ii = 0; ii < numFrames; ++ii ) {
			double frameD = (*trialD)[ii];
			// Don't permit another transient till signal goes < cellThresh
			if (frameD > cellThresh && frameD > lastFrame && !refractory) {
				trialT->push_back( ii );
				refractory = 1;
			} else {
				refractory = (frameD > cellThresh);
			}
			lastFrame = frameD;
		}
	}
	return ret;
}

// Returns transient rates for the current
double findAveTransientRate( const vector< vector< unsigned int > >& transients, double trialDuration, unsigned int numTrials ) 
{
	double ret = 0.0;
	for( auto trial: transients ) {
		for ( auto bin: trial ) {
			ret += bin;
		}
	}
	return ret/ ( trialDuration * numTrials );
}

/// Returns temporal information for the current cell.
double findTemporalInformation( const vector< vector< unsigned int > >& transients, double aveTransientRate, const vector< unsigned int >& shuff, unsigned int numFrames )
{
	double ret = 0.0;
	unsigned int numTrials = transients.size();

	vector< double > binAve( numFrames, 0.0 );
	for( unsigned int trialIdx = 0; trialIdx < numTrials; ++trialIdx ) {
		for ( unsigned int bin: transients[trialIdx] ) {
			unsigned int b = ( bin + shuff[trialIdx] ) % numFrames;
			binAve[b] += 1.0;
		}
	}
	for ( auto bb: binAve ) {
		if ( bb > 0.0 ) {
			ret += bb * log2( bb / aveTransientRate );
		}
	}
	if ( aveTransientRate > 0.0 )
		return ret / ( aveTransientRate * numFrames );
	return 0.0;
}

// Returns the specified percentile value of shuffled TIs, for current cell
double percentileTI( const vector< vector< unsigned int > >& transients, double aveTransientRate, const AnalysisParams& ap, double percentile )
{
	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> shuffler(0, ap.circShuffleFrames );
	unsigned int numTrials = transients.size();
	assert( numTrials > 0 );

	// allocate tiRank[iter#], for sorting.
	vector< double > tiRank( ap.numShuffle ); 

	for ( unsigned int iter = 0; iter < ap.numShuffle; ++iter ) {
		vector< unsigned int > shuff( numTrials );
		for ( unsigned int tt = 0; tt < numTrials; ++tt ) {
			shuff[tt] = shuffler( rng );
		}
		tiRank[iter] = findTemporalInformation( transients, aveTransientRate, shuff, ap.circShuffleFrames  );
	}

	unsigned int percIdx = (percentile * ap.numShuffle ) / 100;
	sort( tiRank.begin(), tiRank.end() );
	return tiRank[percIdx];
}

// Fraction of trials in which there was >= 1 transient. frac[ cell# ]
double fracTrialsFired( const vector< vector< unsigned int > >& transients )
{
	assert( transients.size() > 0 );
	double ret = 0.0;
	for( auto tt : transients ) {
		ret += (tt.size() > 0 );
	}
	return ret / transients.size();
}

/// Computes TI score for a given neuron. zero means it fails a criterion.
// Should really return a lot of stuff: the mean bin vector, the
// baseTI, the PTI, the frac fired, the percentile bin height.
// Suggest return a 2-D vector or a small struct.
CellScore cellTIscore( const vector< vector< double > >& data, const AnalysisParams& ap, const TiAnalysisParams& tip )
{
	CellScore cs;
	unsigned int numTrials = data.size();
	assert( numTrials > 0 );
	unsigned int numFrames = data[0].size();
	vector< unsigned int > nonShuff( numTrials, 0 );
	vector< double > aveOfTrials_ = aveOfTrials( data, nonShuff, ap.circShuffleFrames );

	double cellThresh = findCellThresh( data, ap.circShuffleFrames, tip.transientThresh );
	vector< vector< unsigned int > > transients =
			findTransients( data, cellThresh, ap.circShuffleFrames  );

	double trialDuration = tip.frameDt * ap.circShuffleFrames;
	double aveTransientRate = findAveTransientRate( transients, trialDuration, numTrials );
	cs.baseScore = findTemporalInformation( transients, aveTransientRate, nonShuff, ap.circShuffleFrames  );
	cs.percentileScore = percentileTI( transients, aveTransientRate, ap, tip.tiPercentile );
	cs.fracTrialsFired = fracTrialsFired( transients );

	// Fill up the RNGs
	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> shuffler(0, numFrames / ap.binFrames );
	vector< vector< unsigned int > > binShuff( ap.numShuffle );
	for ( unsigned int ii = 0; ii < ap.numShuffle; ii++ ) {
		binShuff[ii].resize( numTrials );
		for ( unsigned int jj = 0; jj < numTrials; jj++ ) {
			binShuff[ii][jj] = shuffler( rng );
		}
	}
	tuningCurve( data, binShuff, numFrames, ap, cs );

	cs.sigBootstrap = ( cs.baseScore > cs.percentileScore ) && (cs.fracTrialsFired > tip.fracTrialsFiredThresh );

	return cs;
}


// Returns the reliability index.
vector< CellScore > tiScore( py::array_t<double> xs, const AnalysisParams& ap, const TiAnalysisParams& tip )
{
	py::buffer_info info = xs.request();
	auto data = static_cast< double* >( info.ptr);
	unsigned int numCells = 1;
	// unsigned int numDat = info.itemsize;
	unsigned int numTrials = 1;
	unsigned int numFrames = 1;
	assert( info.ndim == 3 );
	numCells = info.shape[2];
	// numDat = info.shape[0] * info.shape[1];
	numTrials = info.shape[1];
	numFrames = info.shape[0];

	vector< vector< vector< double >>> reorderedData;
	reorderData( data, numCells, numTrials, numFrames, reorderedData, ap );

	vector< CellScore > ret( numCells );
	for( unsigned int cellIdx = 0; cellIdx < numCells; cellIdx++ ) {
		ret[cellIdx] = cellTIscore( reorderedData[cellIdx], ap, tip );
	}

	return ret;
}
