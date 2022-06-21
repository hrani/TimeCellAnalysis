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
#include <tcHeader.h>
PYBIND11_MAKE_OPAQUE(std::vector<double>);

#define EPSILON 1.0e-5

/// Returns the mean dfbf[ frame# ] for the specified shuffle[ trial# ]
/// using data[ trial# ][ frame# ]. Operates on data for a single cell.
static vector< double > aveOfTrials( 
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

static vector< unsigned int> findPeaks( const vector< double >& data ) {
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

/// args: data[ frame# ][ trial# ][ cell# ]
/// Pass-back arg: ret[ cell# ][ trial# ][ frame# ]
/// Fills in only the data in the window around the stimuli, from 
/// CS_ONSET_FRAME - CIRC_PAD to US_ONSET_FRAME + CIRC_PAD.
static void reorderData( const double* data, unsigned int numCells, unsigned int numTrials, unsigned int numFrames, vector< vector< vector< double >>>& ret, const AnalysisParams& ap ) 
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

// Returns mean and sdev across all frames and all trials for a given cell.
// data has already been sub-selected for a given cell.
pair< double, double > findCellStats(
		const vector< vector< double > >& data, unsigned int numFrames )
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
	return pair<double, double>( mean, sqrt( sq/numSamples - mean*mean ) );
}

/// Returns time stamps of transients trial# ][ transient# ]. 
// given the reordered dfbf data[ trial# ][ frame# ]
// A given trial may have zero or more transients.
static vector< vector< unsigned int > > findTransients(
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

/// Returns mean and SD of event widths for this cell,
// given the reordered dfbf data[ trial# ][ frame# ]
// Closely based on findTransients()
// A given trial may have zero or more events, it looks at all events.
pair< double, double > findEventWidthStats(
		const vector< vector< double > >& data, double cellThresh, unsigned int numFrames )
{
	double sum = 0.0;
	double sumSq = 0.0;
	double num = 0.0;

	for( auto trialD = data.begin(); trialD != data.end(); ++trialD ) {
		double lastFrame = 0.0;
		unsigned int startFrameOfEvent = 0;
		bool refractory = 0;
		for( unsigned int ii = 0; ii < numFrames; ++ii ) {
			double frameD = (*trialD)[ii];
			// Don't permit another transient till signal goes < cellThresh
			if (frameD > cellThresh && frameD > lastFrame && !refractory) {
				startFrameOfEvent = ii;
				refractory = 1;
			} else if ( refractory ) {
				refractory = (frameD > cellThresh);
				unsigned int dt = ii - startFrameOfEvent;
				if ( dt > 1 && !refractory ) { // End of event. get its stats.
					sum += dt;
					sumSq += dt*dt;
					num += 1.0;
				}
			}
			lastFrame = frameD;
		}
	}

	if ( num > 0.5 ) {
		double mean = sum/num;
		return pair< double, double >( mean, sqrt( sumSq/num - mean*mean) );
	}

	return pair< double, double >( -1.0, -1.0 );
}

// Returns hitTrialRation and the sdev of frames diff between meanPkIdx 
// and nearest transient.
pair< double, double > findHitRatioAndImprecision( unsigned int meanPkIdx, vector< vector< unsigned int > > transients, double hitWindow )
{
	double sumImp = 0.0;
	double sqImp = 0.0;
	double numImp = 0.0;
	double numHits = 0.0;
	for ( auto tt : transients ) {
		unsigned int numT = tt.size();
		if( numT > 0 ) {
			double minDt = 10000000;
			for (auto jj : tt ) {
				double dt = double( jj ) - meanPkIdx;
				if ( abs(minDt) > abs(dt) )
					minDt = dt;
			}
			sumImp += minDt;
			sqImp += minDt*minDt;
			numImp += 1.0;
			numHits += double( abs(minDt) < hitWindow );
		}
	}
	if ( numImp > 0.1 ) {
		double mean = sumImp / numImp;
		double sdev = sqrt( sqImp/numImp - mean*mean);
		return pair< double, double >( numHits/transients.size(), sdev );
	}
	return pair< double, double >( -1.0, -1.0 );
}

/// Computes TI score for a given neuron. zero means it fails a criterion.
// Should really return a lot of stuff: the mean bin vector, the
// baseTI, the PTI, the frac fired, the percentile bin height.
// Suggest return a 2-D vector or a small struct.
CellScore cellPeqScore( const vector< vector< double > >& data, const AnalysisParams& ap, const PeqAnalysisParams& pep )
{
	CellScore cs;
	cs.meanScore = 0.0; // Noise/signal
	cs.baseScore = 0.0; // This is the main return value: Q.
	cs.percentileScore = 0.0; // sdevEW/meanEW
	cs.eventWidthMean = 0.0;
	cs.eventWidthSdev = 0.0;
	cs.imprecision = 0.0;
	unsigned int numTrials = data.size();
	assert( numTrials > 0 );
	// unsigned int numFrames = data[0].size();
	vector< unsigned int > nonShuff( numTrials, 0 );
	cs.meanTrace = aveOfTrials( data, nonShuff, ap.circShuffleFrames );
	vector< unsigned int > pks = findPeaks( cs.meanTrace );
	cs.meanPkIdx = pks[0];

	pair<double, double> stats = findCellStats( data, ap.circShuffleFrames);
	double mean = cs.meanScore = stats.first;
	double sdev = cs.sdev =  stats.second;
	double cellThresh = mean + pep.transientThresh * sdev;

	vector< vector< unsigned int > > transients =
			findTransients( data, cellThresh, ap.circShuffleFrames  );

	// Imprecision is the mean # of frames by which the trial peak differs
	// from the mean peak.
	stats = findHitRatioAndImprecision( cs.meanPkIdx, transients, pep.hitWindow );
	cs.fracTrialsFired = stats.first; // Hit Trial Ratio
	double sdevImp = stats.second;
	double signal = -1.0;
	/*
	 * This doesn't work because mean can be close to zero.
	if ( mean > EPSILON )
		signal = sdev/mean;
	*/
	if ( cs.fracTrialsFired > 0.0 ) {
		// mean signal estimate as average of trials where there is a hit.
		signal = sdev * cs.fracTrialsFired / cs.meanTrace[cs.meanPkIdx];
	}

	// Now event Width
	stats = findEventWidthStats( data, cellThresh, ap.circShuffleFrames );
	double meanEW = stats.first;
	double sdevEW = stats.second;
	cs.eventWidthSdev = sdevEW;
	cs.eventWidthMean = meanEW;
	cs.imprecision = sdevImp;
	/// of these, the only that change with shuffling is HTR. So shuffling
	//is not a good way to analyze this. Just use the peq directly.
	
	double Q = 0.0;
	if ( !(sdevImp < 0.0 || meanEW < 0.0 || signal < 0.0 ) ) { 
		Q = cs.fracTrialsFired * exp( -pep.alpha*signal + pep.beta * sdevEW/meanEW + sdevImp/double(ap.circShuffleFrames) );
	}
	
	cs.baseScore = Q;

	return cs;
}


// Returns the reliability index.
vector< CellScore > peqScore( py::array_t<double> xs, const AnalysisParams& ap, const PeqAnalysisParams& pep )
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
		ret[cellIdx] = cellPeqScore( reorderedData[cellIdx], ap, pep );
	}

	return ret;
}
